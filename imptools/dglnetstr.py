##new### 
def build_network():
    from igem.imptools import dongraphlib as dgl
    from igem.imptools.models import Reaction, Reactant, Product
    myGraph=dgl.network()
    for rxn in Reaction.objects.all():
        print rxn
        for react in rxn.reactants.all():
            myGraph.add_node(str(react.reactname))
            for prod in rxn.products.all():
		        myGraph.add_node(str(prod.prodname))
		        myGraph.add_edge(str(react.reactname),str(prod.prodname),[len(rxn.reactants.all())-1,len(rxn.products.all())-1,rxn.positive_atp(1),rxn.enzymespresent(),0,0],str(rxn)) #new
           #Divide the above by normalization factors, so rxn.enzymespresent()/avg[enzyme] 
            #print rxn,prod,react
        if rxn.reversibility==True:
            for react in rxn.products.all():
                for prod in rxn.reactants.all():
			        myGraph.add_edge(str(react.prodname),str(prod.reactname),[len(rxn.products.all())-1,len(rxn.reactants.all())-1,rxn.positive_atp(-1),rxn.enzymespresent(),0,0],str(rxn))
    finalgraph=finish_network(myGraph)    
    return finalgraph
    
def finish_network(G):
    from igem.imptools import dongraphlib as dgl
    import pickle
    rpair = pickle.load(open('rpair'))
    print 'Finishing'
    for pair in G.edges():
        for edge in G.edge_properties[pair]:
            edge[0][5]=G.node_order(pair[1])
    for pair in G.edges():
        G.edge_properties[pair] = filter(lambda x: x[1] in rpair.keys(), G.edge_properties[pair])
    for comps in G.edges():
        for rxns in G.edge_properties[comps]:
            if comps not in rpair[rxns[1]] and tuple(reversed(list(comps))) not in rpair[rxns[1]]:
                #print comps,rxns,len(G.edges())
                G.edge_properties[comps].remove(rxns) 
                if len(G.edge_properties[comps]) == 0 or len(G.edge_weights(*comps)) == 0:
                    print "edge:",comps, G.edge_properties[comps]
                    G.del_edge(*comps)
                    if comps in G.edges():
                        print "OH HELL NAW"
    for edge in G.edges():
        if G.edge_weights(*edge) == {}:
            G.del_edge(*edge)
    return G 

#weights [num reacts,num prods,atp(positive number), presence of enzymes(!bool), org(0), node_order]
#new weights [num reacts, num prods, atp(pos), enzymeinkegg, org, node_order, inrpair]

