class node:
  def __init__(self, parent, name):
    self.parent = parent
    self.children = []
    self.name = name
  def __str__(self):
    return self.name
  def set_leaf_attribs(self, full):
    self.fullname = full
  def create_child(self, name):
    child = node(self, name)
    self.children.append(child)
    return child

def add_new_orgs():
  from imptools.models import Organism
  count = 0
  f = open('kegg_taxonomy','r')
  for line in f:
    if line[0] != '#':
      org = line.split('\t')
      short = org[1]
      full = org[3]
      org = Organism.objects.filter(ShortName=short)
      if not org.exists():
        new = Organism()
        new.ShortName = short
        new.FullName = full
        new.save()
        count += 1
  return count  

def parse_kegg_taxonomy():
  f = open('kegg_taxonomy','r')
  root = node(None, 'root')
  cur = root
  cur_depth = 0
  for line in f:
    depth = line.count('#')
    if depth == 0:
      # create leaf child of cur
      line = line.split('\t')
      child = cur.create_child(line[1])
      #full = line[3][:line[3].find(' ',line[3].find(' ')+1)]
      child.set_leaf_attribs(line[3].strip())
    else:
      diff = depth - cur_depth
      if diff == 0:
        # create sibling of cur
        cur = cur.parent
        name = line.strip().split()[1]
        cur = cur.create_child(name)
      elif diff == 1:
        # create child of cur
        cur_depth += 1
        name = line.strip().split()[1]
        cur = cur.create_child(name)
      elif diff < 0:
        # create ancestor of cur
        while(diff != 0):
          diff += 1
          cur = cur.parent
        # create sibling of cur
        cur_depth = depth
        cur = cur.parent
        name = line.strip().split()[1]
        cur = cur.create_child(name)
  return root


'''
def recurse(f, cur, depth):
  line = f.readline()
  if line == '':
    return
  elif line.startswith('#'):
    get_
    if(line[depth+1] == '#'):
      #create child node of cur
      name = line.strip().split()[1]
      child = cur.create_child(name)
      recurse(f, child, depth+1)      
      return
    elif(line[depth] == '#'):
      #create sibling of cur (child of parent of cur)
      name = line.strip().split()[1]
      cur.parent.create_child(name)
      recurse(f, child, depth-1)
      return
    else:
      #create sibling of parent (child of parent of parent of cur)
      name = line.strip().split()[1]
      cur.parent.parent.create_child
  else:
    line = line.split('\t')
    child = cur.create_child(line[1])
    full = line[3][:l[3].find(' ',l[3].find(' ')+1)] # a riddle for you
    child.set_leaf_attribs(full)
    return recurse(f, cur, depth)
'''
