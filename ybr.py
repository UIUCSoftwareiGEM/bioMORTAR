import string
import re
from imptools.models import *
from imptools.stoich import *


def formwt(x):
    numelts = len(filter(lambda x: x in string.uppercase, x))
    coeffs =filter(lambda notempty: notempty != '',re.split('[^0-9]',x))
    coeffs = map(int,coeffs)
    return sum(coeffs) + numelts - len(coeffs)

def compwt(x,insigs=[]):
    if x in insigs:
        return 0
    elif int(x.formula[1:]<1000):
	return 0
    else:
        try:
            return formwt(x.formula)
        except:
            return formwt(Compound.objects.get(KeggID=x).formula)

def pathwt(rxn,cmps):
    (mat,seen) = matrix(rxn[:-1],cmps)
#    print mat,seen
#    print rxn, cmps
    weight = 0
    for line in enumerate(mat):
        weight += line[1][-1]*compwt(seen[line[0]])
    return weight


def optimize(rxnb,cmpsb):
	related = []
    aux = (list(),list())
	for rx in Reaction.objects.all():
		if set(cmpsb) & set(rx.products.all()) or set(cmpsb) & set(rx.reactants.all()):
			related.append(rx)
	for rxn in related:
		if pathwt(list(rxnb[:-1])+aux[0]+[rxn]+[1],list(cmpsb)+aux[1]+[rxn.products.all()[0]])>pathwt(rxnb,cmpsb):
			pass
		else:
            aux[0].append(rxn)
            #cmpsb.append(filter(lambda notin: notin not in cmpsb,rxn.products.all())[0])
			aux[1].append(filter(lambda notin: notin not in cmpsb,rxn.products.all())[0])
	return aux


