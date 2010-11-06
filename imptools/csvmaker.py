from igem.imptools.models import Organism, Compound, Reaction, Reactant, Product, Enzyme, Gene

##org=shortname or None
##given rxns a list of strings
##assuming an org
def csv1(writer,rxns, org):
    if not org=='None':
        org=Organism.objects.get(ShortName=org)
        for rx in rxns:
            rx=Reaction.objects.get(KeggID=rx)
            rxFlag=False
            for enz in rx.enzymes.all():
                enzFlag=False
                if enz.genes.filter(organism=org):
                    for gene in enz.genes.filter(organism=org):
                        if not rxFlag and not enzFlag:
                            row=[rx.KeggID , str(stoichlist(rx,True)) , str(stoichlist(rx,False)) , 'EC'+enz.ECnumber , gene.genename ]
                        elif not enzFlag:
                            row=['' , '' , '' , 'EC'+enz.ECnumber , gene.genename ]
                        else:
                            row = ['' , '', '', '', gene.genename] 
                        rxFlag=True
                        enzFlag=True
                        writer.writerow(row)
                else:
                    gene=str(len(enz.genes.all()))
                    if not rxFlag:
                        row=[rx.KeggID , str(stoichlist(rx,True)) , str(stoichlist(rx,False)) , 'EC'+enz.ECnumber , gene+' genes' ]
                    else:
                        row=['' , '' , '' , 'EC'+enz.ECnumber , gene+' genes']
                    rxFlag=True 
                    writer.writerow(row)
    else:
        for rx in rxns:
            rx=Reaction.objects.get(KeggID=rx)
            rxFlag=False
            for enz in rx.enzymes.all():
                gene=str(len(enz.genes.all()))
                if not rxFlag:
                    row=[rx.KeggID , str(stoichlist(rx,True)) , str(stoichlist(rx,False)) , 'EC'+enz.ECnumber , gene+' genes' ]
                else:
                    row=['' , '' , '' , 'EC'+enz.ECnumber , gene+' genes']
                rxFlag=True 
                writer.writerow(row)
                
def stoichlist(rx,react):
    if react:
        list=''
        for comp in rx.reactants.all():
            list+=(str(comp.stoich)+' '+comp.reactname.KeggID+' , ')
    else:
        list=''
        for comp in rx.products.all():
            list+=(str(comp.stoich)+' '+comp.prodname.KeggID+' , ')
    return list
    
    
    
def csv2(writer,rxns,mat,seen):
	rxns.insert(0,'')
	writer.writerow(rxns)
	i=0
	for row in mat:
		row.insert(0,seen[i])
		writer.writerow(row)
		i+=1
