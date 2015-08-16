#for HMRA´s networks
import numpy as np
import re
#Open the file
hmr, rxn, ec, gen, compart, subsys, flux = np.loadtxt("input_file.txt", delimiter="\t", dtype="string", unpack=True)
#metab_repetidos = np.loadtxt("lista_cof.txt", delimiter="\t", dtype="string", unpack=True) #This is optional for open-network
#general dictionary from input file, donde 0->rxn - 1->ec - 2->ge - 3->co - 4->sub
dict_general = {}
for x in xrange(len(hmr)):
    dict_general[hmr[x]] = [rxn[x],ec[x],gen[x],compart[x],subsys[x],flux[x]]

def is_reverisble(expression):
    return " <=> " in expression

def is_irreverisible(expression):
    return " => " in expression

def remove_coefficient(expression):
    return re.sub(r"^((\d+\.)?\d+\s)?","", expression)

def split_reaction(expression):
    return expression.replace(" => ", "/").replace(" <=> ", "/").split("/")
    
def split_metabolites(expression):
    return expression.replace(" + ", "/").split("/")

def split_metabolites_from_list(expression):
    return expression[0].replace(" + ", "/").split("/"), expression[1].replace(" + ", "/").split("/")
    
#Irrev rev dictionary.
irrev_rev = []
dict_irrev_rev = {}
for i in rxn:
    if is_irreverisible(i):
        irrev_rev.append(0)
    elif is_reverisble(i):
        irrev_rev.append(1)
        
for x in xrange(len(hmr)):
    dict_irrev_rev[hmr[x]] = irrev_rev[x]

#Createa dictionary only of: key -> hmr_index and value -> rxn
l_hmr = {}
for x in xrange(len(hmr)):
    l_hmr[hmr[x]] = rxn[x]

#dictionary:E.C. - rxn_izq_der
l_ec_rxn_final = []
l_inter = []
dict_ec_rxn_izq_der = {}

for key in l_hmr:
    l_ec_rxn = split_reaction(l_hmr[key])
    try:
        l_inter = l_ec_rxn[0].replace(" + ", "/").split("/"), l_ec_rxn[1].replace(" + ", "/").split("/")
    except:
        l_inter = l_ec_rxn[0].replace(" + ", "/").split("/")
    for x in xrange(len(l_inter)):
        for y in xrange(len(l_inter[x])):
            l_inter[x][y] = remove_coefficient(l_inter[x][y])
            #l_inter[x][y] = re.sub(r"^[0-9]+\D[0-9]+\s","", l_inter[x][y])
            #l_inter[x][y] = re.sub(r"^[0-9]+\s","", l_inter[x][y])
    dict_ec_rxn_izq_der[key] = list(l_inter)

#Metabolites list.
reaction = [x for x in rxn]
reaction = [x.replace(" + ", "/").replace(" => ", "/").replace(" <=> ", "/") for x in reaction]
reaction = [x.split("/") for x in reaction]
l = list(reaction[0])
a = []
#Counter of metabolites
for x in xrange(0, len(reaction)):
    for j in xrange(0,len(reaction[x])):
        reaction[x][j] = remove_coefficient(reaction[x][j])
        a.append(reaction[x])
    l+=reaction[x]

#Dictionary of hmr_ with metab in a list.
d_hmr = {}
for y in xrange(len(hmr)):
    d_hmr[hmr[y]] = reaction[y]

counts = dict()
for metab in l:
    counts[metab] = counts.get(metab, 0) + 1
#print counts

#organize and eliminate duplicates.
l.sort()
lista_final_metab = list(set(l))
