from igem.imptools.models import Reaction, Reactant, Product

def matrix_for_cobra(rxns, comps, forwards = []):
  matrix, seenCompounds = matrix(rxns, comps, forwards)
  seenCompoundNames = [map(lambda y: y.Name[:-1] if y.Name[-1]==';' else y.Name, x) for x in seenCompounds]
  return matrix, seenCompounds, seenCompoundNames
def matrix(rxs,comps,forwards = []):
    if forwards == []:
        forwards = [0]*len(rxs)
    seenCompounds=[]
    numberOfReactions=len(rxs)
    matrix=[]
    j=0
    #implement fuzzy guessing. I remember a comps string?
    for reaction in rxs:
        if forwards[rxs.index(reaction)] == 0:
            if not Reaction.objects.get(KeggID=reaction).reversibility:
                currentforward = 1
            elif comps[rxs.index(reaction)] in map(lambda name: name.reactname,
                                            Reaction.objects.get(KeggID=reaction).reactants.all()):
                currentforward = 1
            else:
                currentforward = -1
        else:
            currentforward = forwards[rxs.index(reaction)]
        reaction=Reaction.objects.get(KeggID=reaction)
        for reactant in reaction.reactants.all():
            try:
                index=seenCompounds.index(reactant.reactname)
                i=index
            except:
                seenCompounds.append(reactant.reactname)
                matrix.append([0]*(numberOfReactions+1))
                i=len(seenCompounds)-1
            matrix[i][j] -= currentforward*int(reactant.stoich)
            matrix[i][numberOfReactions] -=currentforward* int(reactant.stoich)
        for product in reaction.products.all():
            try:
                index=seenCompounds.index(product.prodname)
                i=index
            except:
                seenCompounds.append(product.prodname)
                matrix.append([0]*(numberOfReactions+1))
                i=len(seenCompounds)-1
            matrix[i][j] += currentforward*int(product.stoich)
            matrix[i][numberOfReactions] +=currentforward* int(product.stoich)
        j+=1
    return (matrix, seenCompounds)
