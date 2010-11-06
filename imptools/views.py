# Create your views here.
from django.http import HttpResponse
from django.core.context_processors import csrf
from django.shortcuts import render_to_response
import pickle
from igem.imptools.models import *
import igem.imptools.dongraphlib as dgl
from igem.imptools.xgmmlcreator import makeXgmml
from igem.imptools.stoich import matrix
#import csv
from igem.imptools import csvmaker
from igem.imptools import biobrick
from igem.ybr import *
import re
global G
G=pickle.load(open('/igem/network0607'))


def welcome(request):
    #return HttpResponse('HI')
    return render_to_response('igemwelcomepage.htm',{})

def start(request):
    c = {}
    all_orgs=Organism.objects.all()
    c.update(csrf(request))
    c.update({'all_orgs':all_orgs})
    return render_to_response('igemstartpage.html',c)

def allresults(request):
    c = {}
    c.update(csrf(request))
    # 'C[0-9]{4}
    wt1=float(request.POST['wt1']) #reactants
    wt2=float(request.POST['wt2']) #products
    wt3=float(request.POST['wt3']) #ATP
    wt4=float(request.POST['wt4']) #Enzyme
    wt5=float(request.POST['wt5']) #Node order
    wt6=float(request.POST['wt6']) #Org
    org=request.POST['orgchoice']
    rmlst=request.POST['rmlst'][:-1]
    starter=request.POST['typein']
    ender=request.POST['typeout']
#search by name here

    try:
        starter=CompoundName.objects.get(Name=starter).comp.KeggID
    except:
        try:
            starter= CompoundName.objects.get(Name=starter + ';').comp.KeggID
        except:
            if not re.match('C[0-9]{5}',starter):
                c.update({'message': 'Search by KeggID or check your spelling.'})
                return render_to_response('igemnopath.html',c)
            #regular expression, also igemnopath
    try:
        ender=CompoundName.objects.get(Name=ender).comp.KeggID
    except:
        try:
            ender = CompoundName.objects.get(Name=ender + ';').comp.KeggID
        except:
            if not re.match('C[0-9]{5}',starter):
                c.update({'message':'Search by KeggID or check your spelling'})
                return render_to_response('igemnopath.html',c)
    num=request.POST['num']
    if eval(num)==1:
        res=dgl.bidirectional_dijkstra(G,starter,ender,[wt1,wt2,wt3,wt4,wt5,wt6],org,eval(rmlst+']'))
        if not res[0]:
            c.update({'message':res[1]})
            return render_to_response('igemnopath.html',c)
        res[2].append(1)
        p1=zip(res[1],res[2], map(lambda Mname: Compound.objects.get(KeggID=Mname).MainName(),res[1] ))
        # return HttpResponse(res[0])
        c.update({'in':starter,'out':ender,'wt1':wt1,'wt2':wt2,'wt3':wt3,'wt4':wt4,'wt5':wt5,'wt6':wt6,'org':org,'removed':rmlst,'p1':p1,'num':num})
        return render_to_response('igemoneresult.html',c)


    else:
        res=dgl.three_paths(G,starter,ender,[wt1,wt2,wt3,wt4,wt5,wt6],org,eval(rmlst+']'))
        if not res[0]:
            c.update({'message':res[1]})
            return render_to_response('igemnopath.html',c)
        res[0][2].append(1)
        res[1][2].append(1)
        res[2][2].append(1)
        p1=zip(res[0][1],res[0][2])
        p2=zip(res[1][1],res[1][2])
        p3=zip(res[2][1],res[2][2])
        c.update({'in':starter,'out':ender,'wt1':wt1,'wt2':wt2,'wt3':wt3,'wt4':wt4,'wt5':wt5,'wt6':wt6,'org':org,'removed':rmlst,'p1':p1,'p2':p2,'p3':p3,'num':num})
        return render_to_response('igemallresults.html',c)


def compound(request,keggid):
    all_names=Compound.objects.get(KeggID=keggid).compoundname_set.all()
    return render_to_response('igemcompoundpage.html',{'all_names':all_names,'kegg':keggid})

def reaction(request,keggid):
    rx=Reaction.objects.get(KeggID=keggid)
    all_reacts=rx.reactants.all()
    all_prods=rx.products.all()
    rev=rx.reversibility
    return render_to_response('igemreactionpage.html',{'all_reacts':all_reacts,'all_prods':all_prods,'rev':rev,'kegg':keggid})

def detail(request):
    c = {}
    c.update(csrf(request))
    path=eval(request.POST['path'])
    wt1=float(request.POST['wt1'])
    wt2=float(request.POST['wt2'])
    wt3=float(request.POST['wt3'])
    wt4=float(request.POST['wt4'])
    wt5=float(request.POST['wt5']) #Node order
    wt6=float(request.POST['wt6']) #Org
    org=request.POST['orgchoice']
    rmlst=request.POST['rmlst'][:-1]
    starter=request.POST['typein']
    ender=request.POST['typeout']
    try:
        starter=CompoundName.objects.get(Name=starter).comp.KeggID
    except:
        try:
            starter= CompoundName.objects.get(Name=starter + ';').comp.KeggID
        except:
            if not re.match('C[0-9]{5}',starter):
                c.update({'message': 'Search by KeggID or check your spelling.'})
                return render_to_response('igemnopath.html',c)
            #regular expression, also igemnopath
    try:
        ender=CompoundName.objects.get(Name=ender).comp.KeggID
    except:
        try:
            ender = CompoundName.objects.get(Name=ender + ';').comp.KeggID
        except:
            if not re.match('C[0-9]{5}',starter):
                c.update({'message':'Search by KeggID or check your spelling'})
                return render_to_response('igemnopath.html',c)

    maincmps,rxns,name=zip(*path)
    (mat,seen)=matrix(rxns[:-1],maincmps)
    pair=zip(mat,seen)
    #from ybr import optimize
    #print optimize(rxns,maincmps)
    #return HttpResponse(path)
    aux = optimize(rxns,maincmps)
    c.update({'aux':aux})
    c.update({'path':path,'in':starter,'out':ender,'wt1':wt1,'wt2':wt2,'wt3':wt3,'wt4':wt4,'wt5':wt5,'wt6':wt6,'org':org,'removed':rmlst,'pair':pair})
    return render_to_response('igemdetailpage.html',c)

def xgmml(request):
    path=eval(request.POST['path'])
    maincmps,rxns,name=zip(*path)
    reactants=[]
    products=[]
    for rx in rxns[:-1]:
        thisrxr=[]
        for react in Reaction.objects.get(KeggID=rx).reactants.all():
            thisrxr.append(str(react))
        reactants.append(thisrxr)
        thisrxp=[]
        for prod in Reaction.objects.get(KeggID=rx).products.all():
            thisrxp.append(str(prod))
        products.append(thisrxp)
    #return render_to_response('igemtest.html',{'cmp':maincmps,'rxn':rxns[:-1],'re':reactants,'pr':products})
    response = HttpResponse(mimetype='text/xgmml')
    response['Content-Disposition'] = 'attachment; filename=imp.xgmml'
    file = makeXgmml(list(maincmps),list(rxns[:-1]),reactants,products)
    for line in file:
        response.write(line)
    return response

def enzyme(request):
    c = {}
    path=eval(request.POST['path'])
    org=request.POST['org']
    if not org=='None':
        org=Organism.objects.get(ShortName=org)
        all_orgs=[]
    else:
        org=False
        all_orgs=Organism.objects.all()
    maincmps,rxns,names=zip(*path)
    rxns=rxns[:-1]
    enzymes=[]
    genes=[]
    for rxn in rxns:
        thisenz=[]
        thisgene=[]
        rx=Reaction.objects.get(KeggID=rxn)
        for enz in rx.enzymes.all():
            thisenz.append(enz.ECnumber)
            if org:
                for gene in enz.genes.filter(organism=org):
                    thisgene.append(gene.genename)
        enzymes.append(thisenz)
        genes.append(thisgene)
    #return HttpResponse(enzymes)
    c.update({'enz':enzymes,'rxn':rxns,'genes':genes,'org':org,'all_orgs':all_orgs,'path':path})
    c.update(csrf(request))
    return render_to_response('igemenzyme.html', c )

def csv(request):
    import csv
    path=eval(request.POST['path'])
    org=request.POST['org']
    maincmps,rxns,name=zip(*path)
    rxns=list(rxns[:-1])
    response = HttpResponse(mimetype='text/csv')
    response['Content-Disposition'] = 'attachment; filename=imp.csv'
    writer = csv.writer(response)
    writer.writerow(['Reaction','Reactants','Products','Enzyme','Genes'])
    csvmaker.csv1(writer,rxns, org)
    writer.writerow(['','','','',''])
    writer.writerow(['','','','',''])
    writer.writerow(['Stoichiometric Matrix'])
    (mat,seen)=matrix(rxns)
    csvmaker.csv2(writer,rxns,mat,seen)
    return response

def biobricker(request,inputmethod):
    if inputmethod=="1":
    	gene=request.POST['gene']
    	seq=biobrick.get_sequence(gene)
    else:
    	seq=request.POST['userdefined'].upper()
    colors=biobrick.check_sequence(seq)
    if colors:
        new=biobrick.add_color(seq,colors)
    else:
        new=biobrick.add_breaks(seq)
    return render_to_response('igembiobrick.html',{'highlighted':new})


def prebiobrick(request):
    c = {}
    c.update(csrf(request))
    return render_to_response('igemgene.html',c)



def team(request):
    return render_to_response('igemteam.html',{})

def imp(request):
    return render_to_response('igemimp.html',{})

def demo(request):
    return render_to_response('igemdemo.html',{})

def feedback(request):
  return render_to_response('igemfeedback.html',{})
