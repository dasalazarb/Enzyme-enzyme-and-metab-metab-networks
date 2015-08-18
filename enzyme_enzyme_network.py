import re

#Location of metabs as reactive or product in rxn irrev and rev.
dict_reactives = {}
dict_products = {}
lista_react = []
lista_prod = []
for i in lista_final_metab:
    lista_react = []
    lista_prod = []
    for key, [react, products] in dict_ec_rxn_izq_der.iteritems():
        if dict_irrev_rev[key] == 1: #es reversible
            for j in react:
                verif_i = " " + i + " "
                verif_j = " " + j + " "
                if verif_i in verif_j:
                    lista_react.append(key)
                    lista_prod.append(key)
            for k in products:
                verif_i = " " + i + " "
                verif_k = " " + k + " "
                if verif_i in verif_k:
                    lista_react.append(key)
                    lista_prod.append(key)
        elif dict_irrev_rev[key] == 0: #es irreversible
            for v in react:
                verif_i = " " + i + " "
                verif_v = " " + v + " "
                if verif_i in verif_v:
                    lista_react.append(key)
            for m in products:
                verif_i = " " + i + " "
                verif_m = " " + m + " "
                if verif_i in verif_m:
                    lista_prod.append(key)
    dict_reactives[i] = lista_react
    dict_products[i] = lista_prod

#interjeccion: metab_rxn_metab

lista_react_metab_product = []
for key in dict_reactives:
    for i in dict_reactives[key]:
        for key_2 in dict_products:
            for j in dict_products[key_2]:
                if key == key_2:# and dict_reactives[key] != [] and dict_products[key_2] != []:
                    if i == j:
                        pass
                    else:
                        if dict_irrev_rev[i] == 1 and dict_irrev_rev[j] == 1:
                            lista_react_metab_product.append(i + "/" + key + "/" + j + "/" + str(dict_irrev_rev[i]) + "/" + str(dict_irrev_rev[j]))
                            lista_react_metab_product.append(j + "/" + key + "/" + i + "/" + str(dict_irrev_rev[j]) + "/" + str(dict_irrev_rev[i]))
                        else:
                            lista_react_metab_product.append(j + "/" + key + "/" + i + "/" + str(dict_irrev_rev[j]) + "/" + str(dict_irrev_rev[i]))

lista_react_metab_product = list(set(lista_react_metab_product))

with open("react_metab_product_network.txt", "w") as f:
    for x in lista_react_metab_product:
        f.write(x+"\n")
    f.close()
