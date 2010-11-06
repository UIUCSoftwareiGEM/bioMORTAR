#This file parses RPAIR dumps.
import urllib
#example = {R00023:[('C00004', 'C00005'),('C00003','C00475')]}

def simple_parse(fname='./rpair_proc'):
  f = open(fname)
  l = f.readlines()
  f.close()
  total = len(l)
  rpair = {}
  names = {}
  pair = ()
  i = 0
  line = l[i]
  while(i != len(l)):
    if i % 1000 < 20:
      print i, ' / ', total
    # read a compound pair
    # get first compound
    lst = line.split()
    first = lst[1]
    name = ' '.join(lst[2:])
    names[first] = name
    # get second compound
    i += 1
    line = l[i]
    lst = line.split()
    second = lst[0]
    name = ' '.join(lst[1:])
    names[second] = name
    if first < second:
      pair = first,second
    else:
      pair = second,first
    # read reactions
    i += 1
    line = l[i]
    # read first rxn line
    lst = line.split()
    for rxn in lst[1:]:
      if rxn in rpair.keys():
        rpair[rxn].append(pair)
      else:
        rpair[rxn] = [pair]
    i += 1
    if i >= total:
      break
    line = l[i]
    while(not line.startswith('COMPOUND')):
      lst = l[i].split()
      for rxn in lst:
        if rxn in rpair.keys():
          rpair[rxn].append(pair)
        else:
          rpair[rxn] = [pair]
      i += 1
      line = l[i]
      if i >= total:
        break
  print 'DONE'
  import cPickle as p
  p.dump(rpair, open('rpair_dict', 'w'))
  p.dump(names, open('compound_names', 'w'))
  return rpair, names

def process_raw(fname='./rpair_raw'):
  f = open(fname)
  l = f.readlines()
  f.close()
  f = open('./rpair_proc','w')
  state = 0
  for index,line in enumerate(l):
    if index % 10000 == 0:
      print index
    if state == 0:
      if line.startswith('COMPOUND'):
        # first compound
        if l[index+1][12] != 'C':
          continued = l[index+1].strip()
          f.write(line.strip() + continued + '\n')
          index += 1
        else:
          f.write(line)
        # second compound
        line = l[index+1]
        if l[index+2][12] != 'C' and l[index+2][0] == ' ':
          continued = l[index+2].strip()
          f.write(line.strip() + continued + '\n')
        else:
          f.write(line)
        state = 1
    elif state == 1:
      if line.startswith('REACTION'):
        #print '1: '+line
        # read first rxn line
        f.write(line)
        # read rest of rxn lines
        next = index+1
        while(True):
          if l[next].startswith('ENZYME') or l[next].startswith('ALIGN'):
            break
          else:
            f.write(l[next])
            next += 1
        state = 0          
  f.close()

def parse_single(file, d):

  entry = "" #RP00001
  name = "" #C00005_C00006
  compound = ('','') #('C00005','C00006')
  reaction = [] #['R00106','etc']
  line = file.readline()
  if len(line) == 0:
      return False
  entry = line[line.find('RP'):line.find('RP')+7]
  #print entry
  line = file.readline()
  name = line[line.find('C'):line.find('C')+13]
  #print name
  compound = (name[:6], name[7:])
  #print compound
  while line.find('REACTION') == -1: #doesn't read the first line of the REACTION table.
    line = file.readline()
    
  line = line[1:]
  while line.find('ALIGN') == -1:
    while line.find('R') != -1:
      reaction.append(line[line.find('R'): line.find('R') + 6])
      line = line[line.find('R') + 6:]
    line = file.readline()
  while  line.find('///') == -1:
    line= file.readline()
  #print reaction
  for react in reaction:
    #print react
    if d.has_key(react):
      if compound not in d[react]:
        d[react].append(compound)
    else:
      d[react] = [compound]
  print entry,name,len(reaction)
  return True
      

def fetchRPAIR(remote = 'ftp://ftp.genome.jp/pub/kegg/ligand/rpair/rpair',
               local = '/home/bhadidi/projects/python/igem/rpair'):
  """Works, but takes some time"""
  rem = urllib.urlopen(remote)
  loc = open(local,'w')
  if file.errors is not None:
    pass
  for line in rem:
    loc.write(line)
  loc.close()
  return local

def update_dict(d, path = '/home/bhadidi/projects/python/igem/rpair'):
  rfile = open(path, 'r')
  while parse_single(rfile,d):
    pass
  
    #print d
