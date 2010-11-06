# Create your views here.
from django.http import HttpResponse, HttpResponseRedirect
from django.template import RequestContext
from django.core.context_processors import csrf
from django.shortcuts import render_to_response
from django.utils import simplejson

from django.core.mail import send_mail

from PIL import Image

import pickle
import jsonpickle
import cairo
import rsvg
from os.path import exists

from igem.strain_designer.forms import FeedbackForm
import igem.parse_kegg_taxonomy

from igem.imptools.models import Organism, Compound, Reaction, Reactant, Product, CompoundName
from igem.strain_designer.models import Plasmid, Promoter, RBS, Terminator
import igem.imptools.dongraphlib as dgl

from reportlab.pdfgen import canvas
from datetime import date

from tempfile import TemporaryFile

from SOAPpy import WSDL

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Restriction import *
from Bio.SeqIO import parse as SeqIOparse
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
Entrez.email = 'bobak.hadidi@gmail.com'
Entrez.tool = 'BioMortar'

global G
G = pickle.load(open('/igem/network1022'))
import cPickle

# mod_wsgi is having difficulty loading /igem/parse_kegg_taxonomy.py
# for some reason, this is required to unpickle tax_root
#global tax_root
#tax_root = cPickle.load(open('/igem/taxonomy','r'))

global RESTRICTION_BATCH
RESTRICTION_BATCH = EcoRI + EcoRV + EcoRII + KpnI + BamHI + \
    HindIII + TaqI + NotI + SpeI + PstI + XbaI + SmaI
global FEATURE_COLORS
FEATURE_COLORS = {
    'promoter': colors.green,
    'rbs': colors.blue,
    'gene': colors.red,
    'terminator': colors.yellow,
    'EcoRI': colors.teal,
    'EcoRV': colors.slategrey,
    'EcoRII': colors.skyblue,
    'KpnI': colors.seashell,
    'BamHI': colors.salmon,
    'HindIII': colors.purple,
    'TaqI': colors.orange,
    'NotI': colors.navy,
    'SpeI': colors.mintcream,
    'PstI': colors.lavender,
    'XbaI': colors.fuchsia,
    'SmaI': colors.cyan
}

global org_short2full
org_short2full = cPickle.load(open('/igem/org_short2fullname','r'))
global plasmids
dplasmids = cPickle.load(open('/igem/pickled_plasmids'))

global SERV
SERV = WSDL.Proxy('http://soap.genome.jp/KEGG.wsdl')

def feedback(request):
  if request.method == 'POST':
    form = FeedbackForm(request.POST)
    if form.is_valid():
      subject = form.cleaned_data['subject']
      message = form.cleaned_data['message']
      sender = form.cleaned_data['sender']
      sender_email = form.cleaned_data['sender_email']
      recipients = ['UIUCSoftwareiGEM@gmail.com', 'bobak.hadidi@gmail.com']
      try:
        send_mail(subject+' from \''+sender+'\' at \''+sender_email+'\'', message, sender, recipients)
      except BadHeaderError:
        return HttpResponse('Invalid header found.')
      return HttpResponseRedirect('/strain_designer/thanks/')
  else:
    form = FeedbackForm()
  return render_to_response('straindes_feedback.html', RequestContext(request, {'form':form}))

def feedback_thanks(request):
  return render_to_response('straindes_feedback_thanks.html',{})

def add_features_to_seq_obj_and_complete_sequence(seq_obj, plasmid, prom, term, rbs):
  # assumes seq_obj.seq holds the gene sequence only (not the completed vector)
  add_plasmid_component_features(seq_obj, plasmid, prom, term, rbs)
  seq_obj.seq = Seq(new_create_sequence(seq_obj, plasmid, prom, term, rbs))
  add_restriction_site_features(seq_obj)

def add_plasmid_component_features(seq_obj, plasmid, prom, term, rbs):
  # assumes seq_obj.seq holds the gene sequence only (not the completed vector)
  prom_end = len(prom.sequence)
  rbs_end = prom_end + len(rbs.sequence)
  gene_end = rbs_end + len(seq_obj.seq)
  term_end = gene_end + len(term.sequence)
  seq_obj.features.append(SeqFeature(id=prom.name,\
      location=FeatureLocation(1+100, prom_end+150), type='promoter', strand=1))
  seq_obj.features.append(SeqFeature(id=rbs.name,\
      location=FeatureLocation(prom_end+1+250, rbs_end+300), type='rbs', strand=1))
  seq_obj.features.append(SeqFeature(id=seq_obj.id,\
      location=FeatureLocation(rbs_end+1, gene_end), type='gene', strand=1))
  seq_obj.features.append(SeqFeature(id=term.name,\
      location=FeatureLocation(gene_end+1+100, term_end+150), type='terminator', strand=1))

def add_restriction_site_features(seq_obj):
  for re,locs in RESTRICTION_BATCH.search(seq_obj.seq, linear=False).items():
    offset = re.elucidate().find('^') + 1
    for loc in locs:
      seq_obj.features.append(SeqFeature(id=str(re),\
          location=FeatureLocation(loc-offset, loc-offset+len(re.site)+50), type='restriction site', strand=-1))

def get_sequence(request, plasmid_id, prom_id, term_id, rbs_id, kegg_gene_id):
  response = HttpResponse(mimetype='text/plain')
  response['Content-Disposition'] = 'attachment; filename=\
    biomortar%s.%s.%s.%s.%s.fasta' % (plasmid_id, prom_id, term_id, rbs_id, kegg_gene_id)
  s = reconstruct_fasta_sequence(plasmid_id, prom_id, term_id, rbs_id, kegg_gene_id)
  response.write(s)
  return response

def reconstruct_fasta_sequence(plasmid_id, prom_id, term_id, rbs_id, kegg_gene_id):
  output = ''
  seq_obj = make_seq_obj_from_bget_result(kegg_gene_id)
  try:
    seq_obj = seq_obj.next()
  except StopIteration:
    seq_obj = None
    output += 'Could not retrieve %s from KEGG. ' % kegg_gene_id
  try:
    plasmid = Plasmid.objects.get(id=int(plasmid_id))
  except Plasmid.DoesNotExist:
    plasmid = None
    output += 'Could not find specified plasmid. '
  try:
    prom = Promoter.objects.get(id=int(prom_id))
  except Promoter.DoesNotExist:
    prom = None
    output += 'Could not find specified promoter. '
  try:
    term = Terminator.objects.get(id=int(term_id))
  except Terminator.DoesNotExist:
    term = None
    output += 'Could not find specified terminator. '
  try:
    rbs = RBS.objects.get(id=int(rbs_id))
  except RBS.DoesNotExist:
    rbs = None
    output += 'Could not find specified rbs. '
  if seq_obj and plasmid and prom and term and rbs:
    sequence = new_create_sequence(seq_obj, plasmid, prom, term, rbs)
    seq_obj.seq = Seq(sequence[0])
    seq_obj.id = ''
    seq_obj.description = 'Biobricks: %s %s %s %s; inset gene: %s' %\
        (plasmid.name, prom.name, term.name, rbs.name, seq_obj.description)
    output = seq_obj.format('fasta')
  return output

def serve_image(request, plasmid_id, prom_id, rbs_id, term_id, kegg_gene_id):
  response = HttpResponse(mimetype='image/png')
  bba_id_str = '.'.join([plasmid_id, prom_id, rbs_id, term_id])
  response['Content-Disposition'] = \
      'attachment; %s.%s.png' % (kegg_gene_id, bba_id_str)
  png = Image.open('/igem/media/plasmids/'+kegg_gene_id+'_'+bba_id_str+'.png')
  png.save(response, 'png')
  return response

def pdf_report(request):
  response = HttpResponse(mimetype='application/pdf')
  response['Content-Disposition'] = \
      'attachment; filename=biomortar%s.pdf' % str(date.today())
  p = canvas.Canvas(response)
  p.drawString(100, 100, "Hello world.")
  p.showPage()
  p.save()
  return response

def email_report(request):
  import subprocess
  c = {}
  if request.method == 'POST':
    inC = request.POST['inputConcen']
    bioC = request.POST['bioConcen']
    (rxns,comps,names) =zip(*request.POST['path'])
    org = request.POST['orgchoice']
    child = subprocess.Popen("strain_designer/cell_modeling_python.py",stdin=PIPE,cwd="/igem")
    child.stdin.write(org)
    child.stdin.write(rxns)
    child.stdin.write(comps)
    child.stdin.write(inC)
    child.stdin.write(bioC)
    child.stdin.flush()
    c.update({"message":"Your request is being processed"})
    c.update(csrf(request))
    return render_to_response("igemnopath.html",c)

def about(request):
  return render_to_response('aboutbiomort.html', {})

def index(request):
  return HttpResponseRedirect("/strain_designer/start")
  #return render_to_response('index.html', {})

def demo(request):
    return render_to_response("biomortdemo.html",{})

def start(request):
  c = {'all_orgs' : Organism.objects.all()}
  c.update(csrf(request))
  return render_to_response('start_strain_des.html', c)

def team(request):
	return render_to_response('biomortarteam.html',{})

def imptools_results(request):
  if(request.POST.has_key('update')):
    return update_results(request)
  c = {}
  c.update(csrf(request))
  promoters = Promoter.objects.all()
  rbs = RBS.objects.all()
  terminators = Terminator.objects.all()
  orgs = Organism.objects.all()
  c['all_orgs'] = orgs
  c['promoters'] = promoters
  c['terminators'] = terminators
  c['RBS'] = rbs
  plasmid_with_descript = [(key, dplasmids[key].description[:-20]) for key in dplasmids.keys()]
  c['plasmids'] = plasmid_with_descript
  wt1=float(request.POST['wt1']) #reactants
  wt2=float(request.POST['wt2']) #products
  wt3=float(request.POST['wt3']) #ATP
  wt4=float(request.POST['wt4']) #Enzyme
  wt5=float(request.POST['wt5']) #Node order
  wt6=float(request.POST['wt6']) #Org
  orgID=request.POST['orgchoice']
  rmlst=request.POST['rmlst'][:-1]
  starter= request.POST['typein']
  ender= request.POST['typeout']
  starter = convert_compound_name_to_id(starter)
  ender = convert_compound_name_to_id(ender)
  num=request.POST['num']
  org_obj = Organism.objects.get(id=orgID)
  if eval(num)==1:
    # res[0] = dist, res[1] = cmpds, res[2] = rxns
    res=dgl.bidirectional_dijkstra(G,starter,ender,[wt1,wt2,wt3,wt4,wt5,wt6],org_obj,eval(rmlst+']'))
    if not res[0]:
      c.update({'message':res[1]})
      return render_to_response('nopath_strain_des.html',c)
    res[2].append(1)
    p1=zip(res[1],res[2],[get_smallest_compound_name_from_kegg_id(comp) for comp in res[1]])
    json_p1 = jsonpickle.encode(p1)
    # return HttpResponse(res[0])
    c.update({'in':starter,'out':ender,'wt1':wt1,'wt2':wt2,'wt3':wt3,'wt4':wt4,'wt5':wt5,'wt6':wt6,'removed':rmlst,'p1':p1,'json_p1':json_p1,'num':num,'org':org_obj})
    return render_to_response('imptools_results_strain_des.html',c)
  else:
    res=dgl.three_paths(G,starter,ender,[wt1,wt2,wt3,wt4,wt5,wt6],org,eval(rmlst+']'))
    if not res[0]:
      c.update({'message':res[1]})
      return render_to_response('nopath_strain_des.html',c)
    res[0][2].append(1)
    res[1][2].append(1)
    res[2][2].append(1)
    p1=zip(res[0][1],res[0][2])
    p2=zip(res[1][1],res[1][2])
    p3=zip(res[2][1],res[2][2])
    c.update({'in':starter,'out':ender,'wt1':wt1,'wt2':wt2,'wt3':wt3,'wt4':wt4,'wt5':wt5,'wt6':wt6,'org':org,'removed':rmlst,'p1':p1,'p2':p2,'p3':p3,'num':num})
    return render_to_response('allresults_strain_des.html',c)

def plasmid_results(request):
  c = {}
  c.update(csrf(request))
  wt1=float(request.POST['wt1']) #reactants
  wt2=float(request.POST['wt2']) #products
  wt3=float(request.POST['wt3']) #ATP
  wt4=float(request.POST['wt4']) #Enzyme
  wt5=float(request.POST['wt5']) #Node order
  wt6=float(request.POST['wt6']) #Org
  plasmid=request.POST['plasmidchoice']
  orgID=request.POST['orgchoice']
  prom=request.POST['promoterchoice']
  rbs=request.POST['rbschoice']
  term=request.POST['terminatorchoice']
  rmlst=request.POST['rmlst'][:-1]
  starter=request.POST['typein']
  ender=request.POST['typeout']
  starter = convert_compound_name_to_id(starter)
  ender = convert_compound_name_to_id(ender)
  num=request.POST['num']
  path = jsonpickle.decode(request.POST['json_p1'])
  cmpds,rxns,names = zip(*path)
  cmpds = list(cmpds)
  rxns = list(rxns)[:-1]
  org_obj = Organism.objects.get(id=orgID)
  if ';' in prom:
    prom = prom[:prom.find(';')]
  if ';' in term:
    term = term[:term.find(';')]
  if ';' in rbs:
    rbs = rbs[:rbs.find(';')]
  if ';' in plasmid:
    plasmid = plasmid[:plasmid.find(';')]
  prom_obj = Promoter.objects.get(name__iexact=prom)
  rbs_obj = RBS.objects.get(name__iexact=rbs)
  term_obj = Terminator.objects.get(name__iexact=term)
  plasmid_obj = Plasmid.objects.get(name__iexact=plasmid)
  c.update({'in':starter,'out':ender,'wt1':wt1,'wt2':wt2,'wt3':wt3,'wt4':wt4,'wt5':wt5,'wt6':wt6,'org':org_obj,'removed':rmlst,'p1':path,'num':num,'plasmid':plasmid_obj,'promoter':prom_obj,'terminator':term_obj,'RBS':rbs_obj})
  enzymes = __get_enzymes__(rxns) # list of list of enzymes for each rxn
  seq_objs, enzyme_chosen_index = __new_get_genes__(enzymes, org_obj)
  for seq_obj in seq_objs:
    add_features_to_seq_obj_and_complete_sequence(seq_obj, plasmid_obj, prom_obj, term_obj, rbs_obj)
  bba_id_str = '.'.join([str(plasmid_obj.id), str(prom_obj.id), str(rbs_obj.id), str(term_obj.id)])
  image_paths = [draw_plasmid(seq_obj,bba_id_str) for seq_obj in seq_objs]
  enzymes = __get_enzyme_at_index__(enzymes, enzyme_chosen_index)
  kegg_gene_ids = [(x.id if x else None) for x in seq_objs]
  orig_orgs = __get_orig_orgs__(kegg_gene_ids)
  c['info'] = zip(rxns, enzymes, kegg_gene_ids, orig_orgs)#, sequences)
  c['image_paths'] = image_paths
  return render_to_response('new_plasmid_results_strain_des.html',c)

def __get_orig_orgs__(kegg_gene_ids):
  orgs = [(x[:x.find(':')] if x else None) for x in kegg_gene_ids]
  orig_orgs = []
  for org_shortname in orgs:
    if org_shortname:
      try:
        org_obj = Organism.objects.get(ShortName=org_shortname)
        orig_orgs.append(org_obj.FullName.strip())
      except:
        orig_objs.append(org)
    else:
      orig_orgs.append(None)
  return orig_orgs

def __get_enzyme_at_index__(enzymes, index):
  retval = []
  for e,i in zip(enzymes, index):
    if e == -1 or e == -2:
      retval.append(e)
    elif i is None:
      retval.append(e[0])
    else:
      retval.append(e[i])
  return retval

def __format_enzymes__(enzymes, index_from_gene):
  formatted_enzymes = []
  for a, i in enumerate(index_from_gene):
    if i == -1 or i == -2:
      formatted_enzymes.append('N/A')
    else:
      selected_enzyme = enzymes[a][i]
      formatted_enzymes.append(selected_enzyme)
  return formatted_enzymes

def __search_entrez__(ec_num, org=None):
  if org!=None:
    handle = Entrez.esearch(db='gene', term=org+'[ORGN] AND '+ec_num+'[ECNO]')
    record = Entrez.read(handle)
    if record['Count'] == '0':
      org == None
  if org == None:
    handle = Entrez.esearch(db='gene', term=ec_num+'[ECNO]')
    record = Entrez.read(handle)
    if record['Count'] == '0':
      return [('N/A', 'No gene information found for enzyme.', ['N/A'])]
  ids = record['IdList']
  handle = Entrez.efetch(db='gene', id=ids[0], rettype='docsum')
  all = handle.read()
  all = all[64:-20] # remove html tags
  org = all[all.find('[')+1:all.find(']')]
  loc = all.find('Official Symbol: ')+17
  symbol = all[loc:all.find(' ',loc)]
  all = all.split('\n')
  nonseq = ['No sequence information retrieved for %s.' % symbol]
  return [(all,nonseq,[org])]

def __new_get_genes__(enzymes, org_obj):
  org = str(org_obj.ShortName)
  seq_objs = []
  e_group_indices = []
  for e_group in enzymes:
    if e_group == -1:
      seq_objs.append(None)
      e_group_indices.append(None)
    elif e_group == -2:
      seq_objs.append(None)
      e_group_indices.append(None)
    else:
      seq_obj, index = get_gene_from_specific_org(e_group, org)
      if seq_obj:
        seq_objs.append(seq_obj)
        e_group_indices.append(index)
      else:
        seq_obj, index = get_gene_from_any_org(e_group)
        seq_objs.append(seq_obj)
        e_group_indices.append(index)
  return seq_objs, e_group_indices

def make_seq_obj_from_bget_result(query):
  if '(' in query:
    query = query[:query.find('(')]
  info = SERV.bget('-f -n n ' + query)
  tf = TemporaryFile()
  tf.write(info)
  tf.seek(0)
  seq_obj = SeqIOparse(tf, 'fasta')
  return seq_obj

def get_gene_from_any_org(e_group):
  seq_obj = None
  for index, e in enumerate(e_group):
    e = str(e)
    all_genes = SERV.bget('ec:'+e)
    seq_obj = parse_bget_all_genes(all_genes)
    if seq_obj:
      try:
        return seq_obj.next(), index
      except:
        continue
  return None, None

def parse_bget_all_genes(all_genes):
  seq_obj = None
  start = all_genes.find('GENES')
  end = all_genes.find('DBLINKS')
  if start != -1:
    if end == -1:
      all_genes = all_genes[start:].lower().split()[1:]
    else:
      all_genes = all_genes[start:end].lower().split()[1:]
    org = all_genes[0]
    gene = all_genes[1]
    query = org+gene
    seq_obj = make_seq_obj_from_bget_result(query)
  return seq_obj

def get_gene_from_specific_org(e_group, org):
  seq_obj = None
  for index, e in enumerate(e_group):
    e = str(e)
    org_genes = SERV.get_genes_by_enzyme('ec:'+e, org)
    if len(org_genes) != 0:
      seq_obj = make_seq_obj_from_bget_result(org_genes[0])
    if seq_obj:
      try:
        return seq_obj.next(), index
      except:
        continue
  return None,None

def __get_genes__(enzymes, org, orgfullname):
  #TODO: this code needs heavy refactoring/cleanup after said refactoring
  org = str(org)
  genes = []
  index_from_genes = []
  from SOAPpy import WSDL
  wsdl = 'http://soap.genome.jp/KEGG.wsdl'
  serv = WSDL.Proxy(wsdl)
  for e_group in enzymes:
    if e_group == -1: # could not find the rxn entry corresponding to this rxn
      genes.append([('N/A', 'No gene information found for enzyme.', ['N/A'])])
      index_from_genes.append(-1)
    elif e_group == -2: # could not find an enzyme corresponding to this rxn
      genes.append([('N/A', 'No gene information found for enzyme.', ['N/A'])])
      index_from_genes.append(-2)
    else:
      found_gene = False
      for index, e in enumerate(e_group): #look for genes from this org for each of the enzymes
        e = str(e)
        org_genes = serv.get_genes_by_enzyme('ec:'+e, org) #args are ec:'*', 'org'
        # enzyme_group is a list of a 3-tuple describing genes for this enzyme
        # first: background info of the gene
        # second: nucleotide sequence
        # third: a list containing the organism which the gene is from,
        #        followed by all organisms which have homologs for this gene
        # e.g. ('info', 'seq', ['escherichia coli', 'homo sapiens', 'mus musculus', ...])
        #enzyme_group = []
        if len(org_genes) != 0:
          # enzyme gene found in selected organism
          #for gene in org_genes:
          #  info = serv.bget('-f -n n ' + gene)
          #  enzyme_group.append(__process_gene_result__(info, orgfullname))
          #genes.append(enzyme_group)
          info = serv.bget('-f -n n ' + org_genes[0])
          genes.append([__process_gene_result__(info, orgfullname)])
          index_from_genes.append(index)
          found_gene = True
          break # do not process other enzymes in this e_group
      if not found_gene:
        for index,e in enumerate(e_group): #look for genes not from this org
          # enzyme gene not found in selected organism, find homolog
          e = str(e)
          all_genes = serv.bget('ec:'+e)
          loc = all_genes.find('GENES')
          if loc == -1:
            #enzyme_group.append(('N/A', 'No gene information found for enzyme.', ['N/A']))
            result = __search_entrez__(e)
            if result[0][0] != 'N/A':
              genes.append(result)
              index_from_genes.append(index)
              found_gene = True
              break
            else:
              continue
            #genes.append([('N/A', 'No gene information found for enzyme.', ['N/A'])])
          else:
            end = all_genes.find('DBLINKS')
            if end != -1:
              all_genes = [x.strip() for x in all_genes[loc+5:end-1].lower().split('\n')]
            else:
              all_genes = [x.strip() for x in all_genes[loc+5:].lower().split('\n')]
            all_genes = map(__format_bget_enzyme_genes__, all_genes)
            shortname = all_genes[0][0]
            info = serv.bget('-f -n n ' + shortname+':'+all_genes[0][1])
            homologs = __get_homologs_orgs__(all_genes[1:])
            if shortname in org_short2full.keys():
              fullname = org_short2full[shortname]
            else:
              fullname = 'Unknown'
            tupl = __process_gene_result__(info, fullname, homologs)
            #enzyme_group.append(tupl)
            genes.append([tupl])
            index_from_genes.append(index)
            found_gene = True
            break
        if not found_gene:
          genes.append([('N/A', 'No gene information found for enzyme.', ['N/A'])])
          index_from_genes.append(-2)
    #genes.append(enzyme_group)
  # ugly code - TODO: should create/change to a 'gene_response' object
  #a list of lists of tuples, each outer entry corresponds to an EC num
  # each inner list has (gene info, sequence, orgs) for each gene for that EC num
  # [ ec:[ gene: (info, seq, orgs) ] ]
  return genes, index_from_genes

def __get_homologs_orgs__(genes):
  retval = []
  for orglist in genes:
    if orglist[0] in org_short2full.keys():
      retval.append(org_short2full[orglist[0]].strip())
  return retval

def __format_bget_enzyme_genes__(string):
  string = string.split()
  retval = [string[0][:-1]] # remove ':' after org id
  for elem in string[1:]:
    loc = elem.find('(') # remove alternate names from gene id
    if loc != -1:
      retval.append(elem[:loc])
    else:
      retval.append(elem)
  return retval

'''
def __get_related_orgs__(org):
  if org == 'eco':
    n = tax_root.children[1].children[0].children[0].children[0].children
    return [x.name for x in n]
  else:
    cPickleearch = [tax_root]
    while(len(search)!=0):
      node = search.pop()
      names = [x.name for x in node.children]
      if org in names:
        sibling_groups = node.parent.children
        for group in sibling_groups:
          if group != node:
            names += __get_leaf_names__(group)
      else:
        search = node.children + search
    return []
'''

def __get_leaf_names__(orggroup):
  pass


def __process_gene_result__(result, orgfullname, homologs=None):
  result = result.split()
  info = [result[0][1:]]
  sequence = []
  one_way_flag = False
  for entry in result[1:]:
    if(one_way_flag):
      sequence.append(entry.strip().upper())
    elif(entry != '(N)'):
      entry = entry.strip()
      if(entry[-1] == ',' or entry[-1] == ';'):
        info.append(entry[:-1])
      else:
        info.append(entry)
    else:
      one_way_flag = True
  #from igem.imptools.biobrick import *
  #colors = check_sequence(sequence)
  #sequence = add_color(sequence, colors)
  if homologs == None:
    return (info, sequence, [orgfullname])
  else:
    return (info, sequence, [orgfullname]+homologs)

def __get_enzymes__(reactions):
  enzyme_list = []
  for rxn in reactions:
    r = Reaction.objects.filter(KeggID__exact=rxn)
    if len(r) == 0:
      enzyme_list.append(-1) # the reaction could be found
      continue
    else:
      r = r[0]
    enzymes = r.enzymes.all()
    if len(enzymes) == 0:
      enzyme_list.append(-2) # no enzymes could be found
      continue
    else:
      enzyme_list.append([e.ECnumber for e in enzymes])
  return enzyme_list

def new_create_sequence(seq_obj, plasmid_obj, prom_obj, term_obj, rbs_obj):
  vector = None
  if seq_obj:
    vector = ''
    vector += prom_obj.sequence.upper()
    vector += rbs_obj.sequence.upper() #short 5' UTR
    vector += seq_obj.seq.tostring().upper()
    vector += term_obj.sequence.upper() # short 3' UTR
    vector += plasmid_obj.sequence.upper()
  return vector

def create_sequence(backbone, enzymes, prom, term, rbs):
  all_sequences = []
  for enzyme, gene, rxn in enzymes:
    cur_vector = ''
    info,gene,orgs = gene[0]
    cur_vector += prom.sequence.upper()
    cur_vector += rbs.sequence.upper() #short 5' UTR
    cur_vector += (''.join(gene)).upper()
    cur_vector += term.sequence.upper() # short 3' UTR
    cur_vector += backbone.sequence.upper()
    #TODO: should add polyA if to express in eukaryotic...
    # buffer between sucessive genes
    all_sequences.append((enzyme, gene, rxn, cur_vector))
  return all_sequences

def create_old_sequence(backbone, enzymes, prom, term, rbs):
  start,end = __get_mcs_site__(backbone)
  sequence = backbone.seq.__str__()[0:start].upper()
  suffix = 'TTAGTTAGTTAG' # biobrick suffix (extra terminators)
  for enzyme, gene, rxn in enzymes:
    info,gene,orgs = gene[0]
    sequence += prom.sequence.upper()
    sequence += rbs.sequence.upper() #short 5' UTR
    sequence += (''.join(gene)).upper()
    sequence += term.sequence.upper() # short 3' UTR
    #TODO: should add polyA if to express in eukaryotic...
    # buffer between sucessive genes
    sequence += suffix
  # add rest of plasmid backbone
  sequence += backbone.seq.__str__()[start:].upper()
  return sequence

def __get_mcs_site__(backbone):
  for feature in backbone.features:
    notes = feature.qualifiers.values()
    for note in notes:
      if 'mcs' in note[0].lower() or 'cloning site' in note[0].lower():
        return int(feature.location.start.__str__()), int(feature.location.end.__str__())
    return 0,0

def __format_genbank__(sequence):
  rec = SeqRecord(Seq(sequence, generic_dna)).format('genbank')
  sequence = []
  one_way_flag = False
  for line in [x.strip() for x in rec.split('\n')]:
    if one_way_flag:
      if line == '//':
        return sequence
      sequence.append(__set_spacing__(line,5))
    elif line.startswith('1'):
      one_way_flag = True
      sequence.append(__set_spacing__(line,5))
  return sequence

def __set_spacing__(line, num):
  if line[num] == ' ':
    return line
  else:
    last = line.find(' ')
    return ('_'*(num-last))+line[:last]+line[last:]


def convert_compound_name_to_id(name):
  aka = name.find('(a.k.a. ')
  name = name[:aka-1] if aka != -1 else name
  c = [x.comp.KeggID for x in CompoundName.objects.filter(Name=name)]
  c = c + [x.comp.KeggID for x in CompoundName.objects.filter(Name=name+';')]
  c = list(set(c))
  return c[0] if c else name


##### ajax #####
def length_comparator_aka(a,b):
  x = a.find('(a.k.a. ')
  y = b.find('(a.k.a. ')
  return 1 if len(a[:x]) > len(b[:y]) else -1

def compound_name_lookup(request):
  results = []
  if request.method == "GET":
    if request.GET.has_key(u'term'):
      value = request.GET[u'term']
      # ignore queries shorter than length 3
      if len(value) > 2:
        results = get_pretty_compound_names(value)
  results.sort(length_comparator_aka)
  json = simplejson.dumps(results)
  return HttpResponse(json, mimetype='application/json')

def length_comparator(a,b):
  return 1 if len(a) > len(b) else -1

# this function should be in imptools.models.Compounds
def get_smallest_compound_name_from_kegg_id(kegg_id):
  comp = Compound.objects.get(KeggID=kegg_id)
  cnames = [x.Name for x in comp.get_names()]
  cnames = map(lambda x: x[:-1] if x[-1] == ';' else x, cnames)
  cnames.sort(length_comparator)
  return cnames[0]

def get_pretty_compound_names(query):
  results = []
  comps = Compound.objects.filter(compoundname__Name__istartswith=query)
  for comp in comps:
    compnames = [x.Name for x in comp.get_names()]
    compnames = map(lambda x: x[:-1] if x[-1] == ';' else x, compnames)
    compnames.sort(length_comparator)
    if len(compnames) > 1:
      descriptive_names = ' or '.join(compnames[1:])
      results.append('%s (a.k.a. %s)' % (compnames[0], descriptive_names))
    else:
      results.append(compnames[0])
  return results

def get_pretty_generic_biobrick_names(query_result_set):
  results = []
  for obj in query_result_set:
    string = obj.name + '; ' + obj.description
    results.append(string)
  return results

def promoter_name_lookup(request):
  results = []
  if request.method == "GET":
    if request.GET.has_key(u'term'):
      value = request.GET[u'term']
      # ignore queries shorter than length 3
      if len(value) > 2:
        query_result_set = Promoter.objects.filter(name__istartswith=value)
        results = get_pretty_generic_biobrick_names(query_result_set)
  results.sort()
  json = simplejson.dumps(results)
  return HttpResponse(json, mimetype='application/json')

def get_plasmid_query_result_set(query, res):
  query_result_set = Plasmid.objects.filter(name__istartswith=query)
  if 'A' in res:
    query_result_set = query_result_set.filter(biobrick_A_res=True)
  if 'C' in res:
    query_result_set = query_result_set.filter(biobrick_C_res=True)
  if 'K' in res:
    query_result_set = query_result_set.filter(biobrick_K_res=True)
  if 'T' in res:
    query_result_set = query_result_set.filter(biobrick_T_res=True)
  return query_result_set

def plasmid_name_lookup(request):
  results = []
  if request.method == "GET":
    res = ''
    if request.GET.has_key(u'res'):
      res = request.GET[u'res']
    if request.GET.has_key(u'term'):
      value = request.GET[u'term']
      # ignore queries shorter than length 3
      if len(value) > 2:
        query_result_set = get_plasmid_query_result_set(value, res)
        results = get_pretty_generic_biobrick_names(query_result_set)
  results.sort()
  json = simplejson.dumps(results)
  return HttpResponse(json, mimetype='application/json')

def terminator_name_lookup(request):
  results = []
  if request.method == "GET":
    if request.GET.has_key(u'term'):
      value = request.GET[u'term']
      # ignore queries shorter than length 3
      if len(value) > 2:
        query_result_set = Terminator.objects.filter(name__istartswith=value)
        results = get_pretty_generic_biobrick_names(query_result_set)
  results.sort()
  json = simplejson.dumps(results)
  return HttpResponse(json, mimetype='application/json')

def rbs_name_lookup(request):
  results = []
  if request.method == "GET":
    if request.GET.has_key(u'term'):
      value = request.GET[u'term']
      # ignore queries shorter than length 3
      if len(value) > 2:
        query_result_set = RBS.objects.filter(name__istartswith=value)
        results = get_pretty_generic_biobrick_names(query_result_set)
  results.sort()
  json = simplejson.dumps(results)
  return HttpResponse(json, mimetype='application/json')


#-outdated-----------------------
def update_results(request):
  c = {}
  c.update(csrf(request))
  starter = request.POST['typein']
  ender = request.POST['typeout']
  wt1 = request.POST['wt1']
  wt2 = request.POST['wt2']
  wt3 = request.POST['wt3']
  wt4 = request.POST['wt4']
  wt5 = request.POST['wt5']
  wt6 = request.POST['wt6']
  org = request.POST['orgchoice']
  rmlst = request.POST['rmlst']
  path = request.POST['path']
  num = request.POST['num']
  enzymes = request.POST['enzymes']
  c.update({'in':starter,'out':ender,'wt1':wt1,'wt2':wt2,'wt3':wt3,'wt4':wt4,'wt5':wt5,'wt6':wt6,'org':org,'removed':rmlst,'p1':path,'num':num})
  #c['enzymes'] = __find_homologous__(enzymes, org)
  c['enzymes'] = enzymes
  return render_to_response('oneresult_strain_des.html',c)

def __find_homologous__(info, org):
  #info = zip(enzymes, genes, res[2])
  for e,genes,r in info:
    # genes is a list of lists of a two-tuple
    for gene in genes:
      for info, seq in gene:
        if info == 'N/A':
          #TODO
          pass

def draw_plasmid(seq_obj, bba_id_str):
  if seq_obj is None:
    return ''
  gd = GenomeDiagram.Diagram(seq_obj.name + ' Biobrick Vector')
  gd_track = gd.new_track(2, name='Restriction Features')
  gd_track_g = gd.new_track(1, name='Plasmid Features')
  gd_feat_set = gd_track.new_set()
  gd_feat_set_g = gd_track_g.new_set()
  for feat in seq_obj.features:
    if feat.type == 'restriction site':
      color = FEATURE_COLORS[feat.id]
      gd_feat_set.add_feature(feat, name=feat.id, color=color, label=True, label_size=25)
    else:
      color = FEATURE_COLORS[feat.type]
      gd_feat_set_g.add_feature(feat, name=feat.id, sigil='ARROW', arrowshaft_height=0.4,\
          arrowhead_length=0.30, color=color, label=True, label_size=25)
  gd.draw(circular=True)
  fpath = '/igem/media/plasmids/'+seq_obj.name+'_'+bba_id_str+'.svg'
  gd.write(fpath, 'SVG')
  return __convert_svg_to_png__(fpath)

def __convert_svg_to_png__(fpath):
  # adapted from : http://guillaume.segu.in/blog/code/43/svg-to-png/
  file = fpath
  svg = rsvg.Handle (file = file)
  if file[-4:] == ".svg":
    file = file[:-4]
  output = "%s.png" % file
  #base = "%s%d.png"
  #i = 1
  #while exists(output):
  #  output = base % (file, i)
  #  i += 1
  width = svg.props.width
  height = svg.props.height
  surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
  cr = cairo.Context (surface)
  wscale = float (width) / svg.props.width
  hscale = float (height) / svg.props.height
  cr.scale (wscale, hscale)
  svg.render_cairo (cr)
  surface.write_to_png(output)
  return output


