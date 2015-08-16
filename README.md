# Enzyme-enzyme-and-metab-metab-networks
This code is useful for obtain enzyme-enzyme and metabolite-metabolite network for HMR and Recon.
Is useful to print the fluxes obtained by FBA.

You can use different atributtes of the nodes and edges. To visualize the network use cytoscape (recommended).

Enzyme-enzyme network: enzyme that produce a metabolite (node) - metabolite produce (edge) - enzyme that use this metabolite (node).
Metab-metab network: metabolite as reactive n reaction x (node) - enzyme of reaction x (edge) - metabolite as product in reaction x.

For reconstuction generated from HMRA:

Input format as txt separated by \tab:

        name of reacion     reaction formula    ec number   gene    compartment    subsystem.

For reconstuction generated from Recon:

Input format as txt separated by \tab:

        name of reacion     reaction formula    ec number   gene    compartment(optional)    subsystem.
