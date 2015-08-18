# Enzyme-enzyme-and-metab-metab-networks
This code is useful for obtain enzyme-enzyme and metabolite-metabolite network for HMR and Recon.
Is useful to print the fluxes obtained by FBA.

I use python language to write this scripts because I used other language but it was very delayed.

You can use different atributtes of the nodes and edges. To visualize the network use cytoscape (recommended).

Enzyme-enzyme network: enzyme that produce a metabolite (node) - metabolite produce (edge) - enzyme that use this metabolite (node).
Metab-metab network: metabolite as reactive n reaction x (node) - enzyme of reaction x (edge) - metabolite as product in reaction x.

For reconstuction generated from HMRA:

For test this script use the example called "net_HMRA_test.txt".

Input format as txt separated by \tab:

        name of rxn     rxn formula    ec number   gene    compartment    subsystem.

For reconstuction generated from Recon:

For test this script use the example called "net_recon_test.txt".

Input format as txt separated by \tab:

        name of rxn     rxn formula    ec number   gene    compartment(opt but with a '-' or 0 value)  subs.


I used the scripts as follows:

        - Net_HMRA.py
        
        - Output_list_metabs.py (optional)
        
        - metab_metab_network.py
        
        - enzyme_enzyme_network.py

Net_Recon.py conatins all the scripts but settled for Recon format.

End
