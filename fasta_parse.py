#! /usr/bin/python

from strain_designer.models import *

if __name__ == "__main__":
  parse()

def parse():
  f = open('all_biobricks.fasta')
  f.seek(0,2)
  eof = f.tell()
  f.seek(0)
  f.readline()
  f.readline()
  plasmids = 0
  terminators = 0
  promoters = 0
  rbs = 0
  while f.tell() != eof:
    header = f.readline()
    assert header[0] == '>'
    sequence = get_sequence(f)
    description, header = get_description(header)
    try:
      name, status, id_num, typ = header[1:].split()
      if typ == 'Plasmid_Backbone':
        create_and_save_plasmid_model(name, status, id_num, description, sequence)
        plasmids += 1
      elif typ == 'Terminator':
        create_and_save_terminator_model(name, status, id_num, description, sequence)
        terminators += 1
      elif typ == 'Regulatory' and ('promoter' in description.lower()):
        create_and_save_promoter_model(name, status, id_num, description, sequence)
        promoters += 1
      elif typ == 'RBS':
        create_and_save_RBS_model(name, status, id_num, description, sequence)
        rbs += 1
      else:
        print 'Unknown biobrick type: %s' % typ
    except:
      print 'Biobrick in non-standard format skipped'
      continue
  print '# Plasmids parsed: %d' % plasmids
  print '# Terminators parsed: %d' % terminators 
  print '# Promoters parsed: %d' % promoters
  print '# RBS parsed: %d' % rbs

def get_description(header):
  descr_start = header.find('"') + 1
  descr_end = header.find('"', descr_start)
  return header[descr_start:descr_end], header[:descr_start-2]

def get_sequence(f):
  sequence = ''
  line = f.readline()
  while line != '\n':
    sequence += line.strip()
    line = f.readline()
  return sequence
  
def assign_standard_biobrick_attr(name, status, id_num, description, sequence, part):
  part.name = name
  part.biobrick_status = status
  part.biobrick_id_num = id_num
  part.description = description
  part.sequence = sequence
  part.isBioBrick = True

def create_and_save_plasmid_model(name, status, id_num, description, sequence):
  p = Plasmid()
  assign_standard_biobrick_attr(name, status, id_num, description, sequence, p)
  #p.save()

def create_and_save_terminator_model(name, status, id_num, description, sequence):
  t = Terminator()
  assign_standard_biobrick_attr(name, status, id_num, description, sequence, t)
  #t.save()

def create_and_save_promoter_model(name, status, id_num, description, sequence):
  p = Promoter()
  assign_standard_biobrick_attr(name, status, id_num, description, sequence, p)
  #p.save()

def create_and_save_RBS_model(name, status, id_num, description, sequence):
  r = RBS()
  assign_standard_biobrick_attr(name, status, id_num, description, sequence, r)
  r.save()

