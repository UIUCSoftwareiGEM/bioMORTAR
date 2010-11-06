def add_orgs():
    import time
    from igem.imptools.models import Reaction, Organism
    
    import pickle
    
    F = open("finalcrazydictionary1.txt","U")
    Di = pickle.load(F)
    F.close()
    z = 1
    #writefile=open("orgwriter.txt",'w')
    
    for K, V in Di.iteritems():
        if K=='ret':
            print time.ctime()
            try:
                org = Organism.objects.get(ShortName = K)
            except:
                print "Could not find this Organism." , K
            RxnList = []
            for item in V:
                try:
                    RxnList.append(Reaction.objects.get(KeggID = item))
                except:
                    print item, 'did not work'
            try:
                org.reaction_set = RxnList
            except:
                print 'Could not add this reaction to this organism.' , K
            org.save()
            del RxnList
            print 'finished' , K, z
        z = z+1
            
    

