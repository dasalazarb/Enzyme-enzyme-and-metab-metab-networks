# Enzyme-enzyme-and-metab-metab-networks
This repository builds enzyme-enzyme and metabolite-metabolite networks for HMR
and Recon reconstructions.  Cytoscape is recommended for visualization.

## Quick start (one-step pipeline)

Use `network_pipeline.py` to run the full workflow end-to-end with Python 3:

```bash
python network_pipeline.py --format hmra --input net_HMRA_test.txt --output-dir outputs
```

or for Recon files:

```bash
python network_pipeline.py --format recon --input net_recon_test.txt --output-dir outputs
```

Outputs are written to the chosen directory:

* `metabolite_counts.txt` – metabolite occurrence counts.
* `lista_metab_rxn_metab.txt` – metabolite-metabolite edges mediated by reactions.
* `react_metab_product_network.txt` – enzyme-enzyme edges via shared metabolites.

## Input format (tab-delimited)

HMRA files (`--format hmra`):

```
name of rxn    rxn formula    ec number   gene    compartment    subsystem    flux(optional)
```

Recon files (`--format recon`):

```
name of rxn    rxn formula    ec number   gene    compartment(optional)    subsystem(optional)
```

Example HMRA input is provided in `net_HMRA_test.txt`.

## Legacy scripts

The original multi-step scripts (`Net_HMRA.py`, `Output_list_metabs.py`,
`metab_metab_network.py`, `enzyme_enzyme_network.py`, `Net_Recon.py`) are kept
for reference but the consolidated pipeline is recommended for new runs.
