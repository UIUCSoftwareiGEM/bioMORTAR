## dongraphlib
## taken extensively from python-graph
##which has this copyright:
# Copyright (c) 2007-2009 Pedro Matiello <pmatiello@gmail.com>
#                         Christian Muise <christian.muise@gmail.com>
#                         Johannes Reinhardt <jreinhardt@ist-dein-freund.de>
#                         Nathan Davis <davisn90210@gmail.com>
#                         Zsolt Haraszti <zsolt@drawwell.net>    
#import numpy ---without numpy
import heapq
import copy
from igem.imptools.models import Reaction, Organism

class network(object):
    def __init__(self):
        self.node_neighbors = {}     # Pairing: Node -> Neighbors
        self.edge_properties = {}    # Pairing: Edge -> (Label, Weight)
        self.node_incidence = {}     # Pairing: Node -> Incident nodes
    def __str__(self):
        return 'a dongraph network'
    def __len__(self):
        return len(self.node_neighbors)
    def __iter__(self):
        for each in self.node_neighbors.iterkeys():
            yield each
    def __getitem__(self, node):
        for each in self.node_neighbors[node]:
            yield each
    def nodes(self):
        return self.node_neighbors.keys()
    def neighbors(self, node,ignored=[]):
        neighs=self.node_neighbors[node][:]
        for ignore in ignored:
            try:
                neighs.remove(ignore)
            except:
                pass
        return neighs
        #return self.node_neighbors[node]
    def incidents(self, node,ignored=[]):
        inc=self.node_incidence[node][:]
        for ignore in ignored:
            try:
                inc.remove(ignore)
            except:
                pass
        return inc
        #return self.node_incidence[node]
    def edges(self):
        return self.edge_properties.keys()
    def has_node(self, node):
        return self.node_neighbors.has_key(node)
    def add_node(self, node, attrs = []):
        if (node not in self.node_neighbors):
            self.node_neighbors[node] = []
            self.node_incidence[node] = []
    def add_nodes(self, nodelist):
        for each in nodelist:
            self.add_node(each)
    def add_edge(self, u, v, wt =[0,0,0,0,0] , rx = ''): #wt is put in as a list
        if (v not in self.node_neighbors[u]):
            self.node_neighbors[u].append(v)
            self.node_incidence[v].append(u)
            self.edge_properties[(u, v)] = [(wt, rx)]
        else:
            self.edge_properties[(u,v)].append((wt,rx))
    def del_node(self, node):
        for each in list(self.incidents(node)):
            self.del_edge(each, node)
        for each in list(self.neighbors(node)):
            self.del_edge(node, each)
        del(self.node_neighbors[node])
        del(self.node_incidence[node])
    def del_edge(self, u, v):
        data=self.edge_properties[(u,v)]####
        self.node_neighbors[u].remove(v)
        self.node_incidence[v].remove(u)
        del(self.edge_properties[(u, v)])
        return data####
    def edge_weights(self, u, v):
        wtdata={}
        for edata in self.edge_properties[(u,v)]:
            wtdata[edata[1]]=edata[0]
        return wtdata
##i've chosen to ignore set edge weight for now
    def edge_labels(self, u, v):
        lbllist=[]
        for edata in self.edge_properties[(u,v)]: 
            lbllist.append(edata[1]) 
        return lbllist
    def has_edge(self, u, v):
        return self.edge_properties.has_key((u, v))
    def node_order(self, node):
        return len(self.neighbors(node))
    def node_degree(self, node):
        return len(self.node_incidence[node])
#node order and degree refer only to the number of neighboring compounds, multiedges are ignored

#i left out traversal but I may need it later
    def calc_weight(self,u,v,punishments,org): #org is new
        #print u,v
        best_weight=None
        if len(self.edge_weights(u,v)) == 0:
            print "EMPTY EDGE",u,v,punishments,org
        for kee,val in self.edge_weights(u,v).iteritems():   
            #array=numpy.array(copy.copy(val)) 
            array=copy.copy(val)
            if org and (punishments[4]>0):
                try:
                    Reaction.objects.get(KeggID=kee).organisms.get(ShortName=org)
                    array[4]=0
                except:
                    array[4]=1
            else:
                array[4]=1
            if best_weight: #this will compare multiple weights
                #cur_weight=numpy.dot(punishments,array)
               	cur_weight=dot(punishments,array)+1
                if cur_weight < best_weight[0]:
                    best_weight=(cur_weight,kee)
            else:
                #best_weight=(numpy.dot(punishments,array),kee) #punishments and weights are expected to be numpy arrays
                best_weight=(dot(punishments,array)+1,kee)
        return best_weight


def bidirectional_dijkstra(G, source, target,punishments=[0,0,0,0,0,0],org=False,ignored=[]):
    try:
    	G.node_neighbors[source]
    except:
    	return (False, 'Source compound has no reactions')
    try:
    	G.node_incidence[target]
    except:
    	return (False, 'Target compound has no reactions')
    if source == target: return (0, [source])
    if org=='None':
        org=False
    #Init:   Forward             Backward
    dists =  [{},                {}]# dictionary of final distances
    paths =  [{source:[source]}, {target:[target]}] # dictionary of paths 
    rxns= [{source:[]}, {target:[]}]
    fringe = [[],                []] #heap of (distance, node) tuples for extracting next node to expand
    seen =   [{source:0},        {target:0} ]#dictionary of distances to nodes seen 
    #initialize fringe heap
    heapq.heappush(fringe[0], (0, source)) 
    heapq.heappush(fringe[1], (0, target))
    #neighs for extracting correct neighbor information
    neighs = [G.neighbors, G.incidents]
    #variables to hold shortest discovered path
    #finaldist = 1e30000
    finalpath = []
    finalrxns=[]
    #puns=numpy.array(punishments)
    puns=punishments
    dir = 1
    while fringe[0] and fringe[1]:
        # choose direction 
        # dir == 0 is forward direction and dir == 1 is back
        dir = 1-dir
        # extract closest to expand
        (dist, v )= heapq.heappop(fringe[dir]) 
        if v in dists[dir]:
            # Shortest path to v has already been found 
            continue
        # update distance
        dists[dir][v] = dist #equal to seen[dir][v]
        if v in dists[1-dir]:
            # if we have scanned v in both directions we are done 
            # we have now discovered the shortest path
            return (finaldist,finalpath,finalrxns)
        for w in neighs[dir](v,ignored):
            if(dir==0): #forward
               # print "yay1"
                (wt,rx)=G.calc_weight(v, w,puns, org)
                vwLength = dists[dir][v] + 1+wt
            else: #back, must remember to change v,w->w,v
                (wt,rx)=G.calc_weight(w, v,puns, org)
                vwLength = dists[dir][v] + 1+wt
            
            if w in dists[dir]:
                if vwLength < dists[dir][w]:
                    raise ValueError,\
                        "Contradictory paths found: negative weights?"
            elif w not in seen[dir] or vwLength < seen[dir][w]:
                # relaxing        
                seen[dir][w] = vwLength
                heapq.heappush(fringe[dir], (vwLength,w)) 
                paths[dir][w] = paths[dir][v]+[w]
                rxns[dir][w] = rxns[dir][v] + [rx]
                if w in seen[0] and w in seen[1]:
                    #see if this path is better than than the already
                    #discovered shortest path
                    totaldist = seen[0][w] + seen[1][w] 
                    if finalpath == [] or finaldist > totaldist:
                        finaldist = totaldist
                        revpath = paths[1][w][:]
                        revrxns= rxns[1][w][:]
                        revpath.reverse()
                        revrxns.reverse()
                        finalrxns=rxns[0][w] + revrxns[:]
                        finalpath = paths[0][w] + revpath[1:]
    return (False,'No path found')




def three_paths(G, source, target,punishments=[0,0,0,0,0],org=False,ignored=[]):
    best= bidirectional_dijkstra(G, source, target,punishments,org,ignored)
    heap=[]
    if not best[0]:
        return (False,best[1])
    for x in xrange(0,len(best[2])):
        ux=best[1][x]
        vx=best[1][x+1]
        datax=G.del_edge(ux,vx)
        xpath=bidirectional_dijkstra(G, source, target,punishments,org,ignored)
        if not (xpath==False or xpath in heap):
            heapq.heappush(heap,xpath)
        for dx in datax:
            G.add_edge(ux,vx,dx[0],dx[1])
            ###########
    secondbest=heapq.heappop(heap)
    for y in xrange(0,len(secondbest[2])):
        uy=xpath[1][y]
        vy=xpath[1][y+1]
        datay=G.del_edge(uy,vy)
        ypath=bidirectional_dijkstra(G, source, target,punishments,org,ignored)
        if not (ypath==False or ypath in heap):
            heapq.heappush(heap,ypath)
        for dy in datay:
            G.add_edge(uy,vy,dy[0],dy[1])
    thirdbest=heapq.heappop(heap)
    return (best,secondbest,thirdbest)
    
    

            
def dot(arr1,arr2):
    result=0
    for i in xrange(0,len(arr1)):
        result+=arr1[i]*arr2[i]
    return result
