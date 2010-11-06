import cPickle as p

def parse_kegg_enzyme():
  '''
  returns a dictionary of ec id ('3.2.2.1') -> org, geneID ('eco', 'b0002')
  '''
  try:
    f = open('enzyme_dict', 'r')
    return p.load(f)
  except:
    try:
      f = open('kegg_enzyme', 'r')
    except:
      print 'could not find pickled dict nor original kegg dump'
      return None
    eDict = {}
    # parse the file
    line = f.readline()
    while(True):
      if line == '':
        break
      elif line.startswith('ENTRY'):
        ec_num = line.strip().split()[2]
        line = f.readline()
        while(not(line.startswith('ENTRY') or line.startswith('GENES'))):
          line = f.readline()
        if line.startswith('ENTRY'):
          continue #found new ec_num
        else: #found genes
          line = line.
          #TODO
        
    # done
    f = open('enzyme_dict', 'w')
    p.dump(eDict, f)
    return eDict

