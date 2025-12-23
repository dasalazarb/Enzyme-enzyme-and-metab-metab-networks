"""End-to-end network construction pipeline.

This module replaces the previous multi-step workflow that required executing
`Net_HMRA.py`, `metab_metab_network.py`, `enzyme_enzyme_network.py` and
`Output_list_metabs.py` separately.  It exposes a single CLI that reads one
input file and writes all expected artifacts in a single pass.

Usage (HMRA files):

    python network_pipeline.py --format hmra --input net_HMRA_test.txt --output-dir outputs

Usage (Recon files):

    python network_pipeline.py --format recon --input net_recon_test.txt --output-dir outputs

The code keeps the original output semantics while avoiding global state and
supporting Python 3.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


ARROW_REV = "<=>"
ARROW_IRREV = "=>"


def remove_coefficient(metabolite: str) -> str:
    """Strip an optional numeric stoichiometric coefficient from a metabolite string."""

    cleaned = metabolite.strip()
    return re.sub(r"^((\d+\.)?\d+\s+)?", "", cleaned)


def split_reaction(reaction: str) -> Tuple[List[str], List[str], bool]:
    """Split a reaction string into reactants, products and directionality.

    Reactions may contain either reversible ("<=>") or irreversible ("=>") arrows.
    Coefficients are removed and metabolites are returned as plain strings.
    """

    normalized = reaction.strip()
    reversible = ARROW_REV in normalized

    if reversible:
        left, right = normalized.split(ARROW_REV)
    elif ARROW_IRREV in normalized:
        left, right = normalized.split(ARROW_IRREV)
    else:
        raise ValueError(f"Reaction '{reaction}' does not contain a supported arrow")

    def _split_side(side: str) -> List[str]:
        if not side.strip():
            return []
        metabolites = re.split(r"\s*\+\s*", side.strip())
        return [remove_coefficient(met) for met in metabolites if met]

    return _split_side(left), _split_side(right), reversible


@dataclass
class ReactionRecord:
    identifier: str
    reaction: str
    ec_number: str
    gene: str
    compartment: str
    subsystem: str
    flux: str = ""
    reactants: List[str] = field(default_factory=list)
    products: List[str] = field(default_factory=list)
    reversible: bool = False

    @classmethod
    def from_row(cls, row: Sequence[str], dataset: str) -> "ReactionRecord":
        if dataset == "hmra":
            if len(row) < 6:
                raise ValueError("HMRA rows must contain at least 6 tab-separated fields")
            identifier, reaction, ec_number, gene, compartment, subsystem = row[:6]
            flux = row[6] if len(row) > 6 else ""
        else:  # recon
            if len(row) < 5:
                raise ValueError("Recon rows must contain at least 5 tab-separated fields")
            identifier, reaction, ec_number, gene = row[:4]
            compartment = row[4] if len(row) > 4 else ""
            subsystem = row[5] if len(row) > 5 else ""
            flux = ""

        reactants, products, reversible = split_reaction(reaction)

        return cls(
            identifier=identifier.strip(),
            reaction=reaction.strip(),
            ec_number=ec_number.strip(),
            gene=gene.strip(),
            compartment=compartment.strip(),
            subsystem=subsystem.strip(),
            flux=flux.strip(),
            reactants=reactants,
            products=products,
            reversible=reversible,
        )


def load_records(path: Path, dataset: str) -> List[ReactionRecord]:
    records: List[ReactionRecord] = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row or all(not cell.strip() for cell in row):
                continue
            try:
                records.append(ReactionRecord.from_row(row, dataset))
            except Exception as exc:  # pragma: no cover - defensive logging
                raise ValueError(f"Error parsing line {line_number} of {path}: {exc}") from exc

    if not records:
        raise ValueError(f"No records were found in {path}")
    return records


def metabolite_counts(records: Iterable[ReactionRecord]) -> Counter:
    counts: Counter = Counter()
    for record in records:
        counts.update(record.reactants)
        counts.update(record.products)
    return counts


def metabolite_network(records: Iterable[ReactionRecord]) -> List[str]:
    """Build metabolite-to-metabolite links mediated by each reaction."""

    edges: set[str] = set()
    for record in records:
        for reactant in record.reactants:
            for product in record.products:
                payload = f"{reactant}/{record.identifier}/{product}/{int(record.reversible)}/{record.flux}"
                edges.add(payload)
                if record.reversible:
                    reverse_payload = f"{product}/{record.identifier}/{reactant}/{int(record.reversible)}/{record.flux}"
                    edges.add(reverse_payload)
    return sorted(edges)


def enzyme_network(records: Iterable[ReactionRecord]) -> List[str]:
    """Build enzyme-enzyme edges via shared metabolites.

    The implementation mirrors the legacy "enzyme_enzyme_network.py" behaviour while
    removing the reliance on global variables and Python 2 constructs.
    """

    records_by_metabolite_react: Dict[str, List[ReactionRecord]] = defaultdict(list)
    records_by_metabolite_prod: Dict[str, List[ReactionRecord]] = defaultdict(list)

    for record in records:
        target_react = records_by_metabolite_react
        target_prod = records_by_metabolite_prod

        for metabolite in record.reactants:
            target_react[metabolite].append(record)
            if record.reversible:
                target_prod[metabolite].append(record)

        for metabolite in record.products:
            target_prod[metabolite].append(record)
            if record.reversible:
                target_react[metabolite].append(record)

    edges: set[str] = set()
    common_metabolites = set(records_by_metabolite_react) | set(records_by_metabolite_prod)

    for metabolite in sorted(common_metabolites):
        for react_record in records_by_metabolite_react.get(metabolite, []):
            for prod_record in records_by_metabolite_prod.get(metabolite, []):
                if react_record.identifier == prod_record.identifier:
                    continue
                if react_record.reversible and prod_record.reversible:
                    edges.add(
                        f"{react_record.identifier}/{metabolite}/{prod_record.identifier}/"
                        f"{int(react_record.reversible)}/{int(prod_record.reversible)}"
                    )
                    edges.add(
                        f"{prod_record.identifier}/{metabolite}/{react_record.identifier}/"
                        f"{int(prod_record.reversible)}/{int(react_record.reversible)}"
                    )
                else:
                    edges.add(
                        f"{prod_record.identifier}/{metabolite}/{react_record.identifier}/"
                        f"{int(prod_record.reversible)}/{int(react_record.reversible)}"
                    )
    return sorted(edges)


def write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(f"{line}\n")


def write_counts(path: Path, counts: Counter) -> None:
    ordered = sorted(counts.items(), key=lambda item: (item[1], item[0]))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for metabolite, count in ordered:
            handle.write(f"{metabolite}/{count}\n")


def run_pipeline(dataset: str, input_path: Path, output_dir: Path) -> None:
    records = load_records(input_path, dataset)
    counts = metabolite_counts(records)

    metabolite_edges = metabolite_network(records)
    enzyme_edges = enzyme_network(records)

    write_counts(output_dir / "metabolite_counts.txt", counts)
    write_lines(output_dir / "lista_metab_rxn_metab.txt", metabolite_edges)
    write_lines(output_dir / "react_metab_product_network.txt", enzyme_edges)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build enzyme and metabolite networks in one pass")
    parser.add_argument(
        "--format",
        choices=["hmra", "recon"],
        required=True,
        help="Input file format. Use 'hmra' for HMRA reconstructions and 'recon' for Recon files.",
    )
    parser.add_argument("--input", type=Path, required=True, help="Path to the tab-delimited input file")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="Directory where result files will be written",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_pipeline(args.format, args.input, args.output_dir)


if __name__ == "__main__":
    main()
