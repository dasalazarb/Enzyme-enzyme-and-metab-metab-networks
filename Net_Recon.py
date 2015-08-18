#For recon al scripts in one.

import numpy as np
import re
#abrir archivo, separar por delimiter, emplear string.
hmr, rxn, ec, ge, co, sub = np.loadtxt("n_g_recon.txt", delimiter="\t", dtype="string", unpack=True)

def is_reversible_recon(expression):
    return "<=>" in expression

def is_irreversible_recon(expression):
    return " =>" in expression

def remove_coefficient(expression):
    return re.sub(r"^((\d+\.)?\d+\s)?","", expression)

def split_reaction(expression):
    return expression.replace(" => ", "/").replace(" <=> ", "/").split("/")
    
def split_metabolites(expression):
    return expression.replace(" + ", "/").split("/")

def split_metabolites_from_list(expression):
    return expression[0].replace(" + ", "/").split("/"), expression[1].replace(" + ", "/").split("/")
    
#General dictionary, donde 0->rxn - 1->ec - 2->ge - 3->co - 4->sub
dict_general = {}
for x in xrange(len(hmr)):
    dict_general[hmr[x]] = [rxn[x],ec[x],ge[x],co[x],sub[x]]

#Irrev rev dictionary.
irrev_rev = []
dict_irrev_rev = {}
for i in rxn:
    if is_reversible_recon(i):
        irrev_rev.append(1)
    elif is_irreversible_recon(i):
        irrev_rev.append(0)

for x in xrange(len(hmr)):
    dict_irrev_rev[hmr[x]] = irrev_rev[x]

#crear un diccionario: como key -> hmr_index como value -> rxn
l_hmr = {}
for x in xrange(len(hmr)):
    l_hmr[hmr[x]] = rxn[x].strip()

#dictionary which relates: E.C. - rxn_izq_der. 
#for rxns ('utp[e] <=> ') the solution is: [['utp[e]'], ['']]
l_ec_rxn_final = []
l_inter = []
dict_ec_rxn_izq_der = {}

for key in l_hmr:
    para_recon = " " + l_hmr[key] + " "
    l_ec_rxn = split_reaction(para_recon)
    try:
        l_inter = l_ec_rxn[0].replace(" + ", "/").split("/"), l_ec_rxn[1].replace(" + ", "/").split("/")
    except:
        l_inter = l_ec_rxn[0].replace(" + ", "/").split("/")
    for x in xrange(len(l_inter)):
        for y in xrange(len(l_inter[x])):
            l_inter[x][y] = l_inter[x][y].strip()
            l_inter[x][y] = remove_coefficient(l_inter[x][y])
    dict_ec_rxn_izq_der[key] = list(l_inter)
    #l_ec_rxn_final.append(l_inter)
    
#dictionary that has as values two lists, first list contains reactives and second contains products.
rxn_rev = list()
rxn_irrev = list()
rxn_irrev_final = dict()
for x in rxn:
    if "<=>" in x:
        rxn_rev = x.replace("<=>","/").split("/")
    elif "=>" in x:
        rxn_irrev = x.replace("=>","/").split("/")

#Metabolite list.
reaction = [x for x in rxn]
reaction = [x.replace(" + ", "/").replace("<=>", "/").replace("=>", "/") for x in reaction]
reaction = [x.split("/") for x in reaction]
l = list(reaction[0])
#a = []

#Counter of metabs.
for x in xrange(0, len(reaction)):
    for j in xrange(0,len(reaction[x])):
        reaction[x][j] = reaction[x][j].strip()
        reaction[x][j] = remove_coefficient(reaction[x][j])
    l+=reaction[x]

#dictcionary: hmr with metab as a list
d_hmr = {}
for y in xrange(len(hmr)):
    d_hmr[hmr[y]] = reaction[y]
d_recon_hmr = {}
for key in d_hmr:
    if '' in d_hmr[key]:
        d_hmr[key].sort(reverse=True)
        d_hmr[key].pop()
        d_recon_hmr[key] = d_hmr[key]
    else:
        d_recon_hmr[key] = d_hmr[key]
        
#dict of counts of metabs
counts = dict()
for metab in l:
    counts[metab] = counts.get(metab, 0) + 1

#Clean result.
l.sort()
lista_final_metab = list(set(l))

#######################################
#metab_metab_network

p = []
wi = []
valor_interact = 3
for i in lista_final_metab:
    for key in d_recon_hmr:
        for e in xrange(len(d_recon_hmr[key])):
            if i in d_recon_hmr[key]:
                if i == d_recon_hmr[key][e]:
                    pass
                else:
                    if counts[i] <= valor_interact:
                        try:
                            wi.append(d_recon_hmr[key][e] + "/" + "(000AAA000)" + key + "/" + i)
                            #p.append(j[e] + "/" + i)
                        except:
                            pass

l_final = []
l_final_2 = []
for i in wi:
    a = i.split("/")
    a.sort()
    if a[1] == a[2]:
        pass
    else:
        l_final.append(a)

for t in xrange(0,len(l_final)):
    if re.match(r"^" + "n_EX", l_final[t][0].replace("(000AAA000)","")) or re.match(r"^" + "g_EX", l_final[t][0].replace("(000AAA000)","")):
        l_final_2.append("EX_" + l_final[t][1] + "/" + l_final[t][0].replace("(000AAA000)","") + "/" + "EX_" + l_final[t][2] + "/" + 
        dict_general[l_final[t][0].replace("(000AAA000)","")][3] + "/"+ dict_general[l_final[t][0].replace("(000AAA000)","")][4])
    elif re.match(r"^" + "n_", l_final[t][0].replace("(000AAA000)","")):
        l_final_2.append("n_" + l_final[t][1] + "/" + l_final[t][0].replace("(000AAA000)","") + "/" + "n_" + l_final[t][2] + "/" + 
        dict_general[l_final[t][0].replace("(000AAA000)","")][3] + "/"+ dict_general[l_final[t][0].replace("(000AAA000)","")][4])
    elif re.match(r"^" + "g_", l_final[t][0].replace("(000AAA000)","")):
        l_final_2.append("g_" + l_final[t][1] + "/" + l_final[t][0].replace("(000AAA000)","") + "/" + "g_" + l_final[t][2] + "/" + 
        dict_general[l_final[t][0].replace("(000AAA000)","")][3] + "/"+ dict_general[l_final[t][0].replace("(000AAA000)","")][4])

l_final = list(set(l_final_2))

with open("metab_metab_network.txt", "w") as f:
    for x in l_final:
        f.write(x+"\n")
    f.close()

#######################################
#rxn_rxn-network

lista_hmr = []
valor_interact = 10
for key in dict_ec_rxn_izq_der:
    for i in dict_ec_rxn_izq_der[key][0]:
        for key_2 in dict_ec_rxn_izq_der:
            if counts[i] <= valor_interact:
                if i in dict_ec_rxn_izq_der[key_2][1]:
                    lista_hmr.append(key_2 + "/1/" + i + "/2/" + key)
                else:
                    pass
            else:
                pass
    for a in dict_ec_rxn_izq_der[key][1]:
        for key_2 in dict_ec_rxn_izq_der:
            if counts[i] <= valor_interact:
                if j in dict_ec_rxn_izq_der[key_2][0]:
                    lista_hmr.append(key + "/1/" + j + "/2/" + key_2)
                else:
                    pass
            else:
                pass

lista_hmr = list(set(lista_hmr))
with open("hmr_hmr_network.txt", "w") as f:
    for x in lista_hmr:
        f.write(x+"\n")
    f.close()
