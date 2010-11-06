from igem.imptools.models import Reaction

# To use this function make sure that rx and rxdir are initialized to the
# values of the first reaction in the list. 
def set_revers():
    file=open('reaction_mapformula.lst')
    rx='R00005' 
    rxdir='r'
    count=0
    for line in file:
      	if line.startswith(rx):
         	arrowloc=line.find('=')
			if (line[arrowloc-1]=='<' and line[arrowloc+1]=='>'):
          		arrow='b'	
			elif (line[arrowloc-1]=='<'):
          		arrow='l'	
			else:
         		arrow='r'
     		if arrow!=rxdir:
           		rxdir='b'
	    else:
			if rxdir=='b':
		    	try:
                  	dbrx=Reaction.objects.get(KeggID=rx)
					dbrx.reversibility=True
					dbrx.save()
					print 'this reaction is reversible',
					print rx
		    	except:
                  	print 'I couldnt change the db for this one ',
					print rx 
			elif rxdir=='l':
              	count+=1
               	switch(rx)
              	print 'this one has been switched ',
		    	print rx
		    	print count
			else:
              	print 'stayed right',
              	print rx 
		rx=line[:6]
		arrowloc=line.find('=')
		if (line[arrowloc-1]=='<' and line[arrowloc+1]=='>'):
           	arrow='b'
		elif (line[arrowloc-1]=='<'):
          	arrow='l'
		else:
		   	arrow='r'
		rxdir=arrow


def switch(rx):
    newprod=[]
    newreact=[]
    try:
    	rxn=Reaction.objects.get(KeggID=rx)
    	for elt in rxn.reactants.all():
        	newprod.append(elt)
    	for elt in rxn.products.all():
        	newreact.append(elt)
    	rxn.reactants=newreact
    	rxn.products=newprod
    except:
    	print 'couldnt switch',rx
