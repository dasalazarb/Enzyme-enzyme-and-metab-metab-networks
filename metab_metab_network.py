#metabolite-metabolite network
lista_metab_rxn_metab = []
for key in dict_ec_rxn_izq_der:
    for i in dict_ec_rxn_izq_der[key][0]:
        for j in dict_ec_rxn_izq_der[key][1]:
            if dict_irrev_rev[key] == 0:
                lista_metab_rxn_metab.append(i + "/" + key + "/" + j + "/" + str(dict_irrev_rev[key]) + "/" + dict_general[key][5])
            elif dict_irrev_rev[key] == 1:
                lista_metab_rxn_metab.append(i + "/" + key + "/" + j + "/" + str(dict_irrev_rev[key]) + "/" + dict_general[key][5])
                lista_metab_rxn_metab.append(j + "/" + key + "/" + i + "/" + str(dict_irrev_rev[key]) + "/" + dict_general[key][5])
lista_metab_rxn_metab = list(set(lista_metab_rxn_metab))

with open("lista_metab_rxn_metab.txt", "w") as f:
    for x in lista_metab_rxn_metab:
        f.write(x+"\n")
    f.close()
