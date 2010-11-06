from django.db import models
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import ExactPosition, FeatureLocation, SeqFeature
from Bio.Seq import Seq

class RBS(models.Model):
  sequence = models.CharField(max_length=2500)
  name = models.CharField(max_length=25, unique=True)
  isBioBrick = models.BooleanField()
  description = models.CharField(max_length=200)
  isConstitutive = models.BooleanField()
  isProkaryotic = models.BooleanField()
  specificOrganism = models.CharField(max_length=100)
  biobrick_id_num = models.PositiveIntegerField(null=True)
  biobrick_status = models.CharField(max_length=10, null=True)
  def __unicode__(self):
    return self.name +' '+ self.description

class Terminator(models.Model):
  sequence = models.CharField(max_length=2500)
  name = models.CharField(max_length=25, unique=True)
  isBioBrick = models.BooleanField()
  description = models.CharField(max_length=200)
  isForward = models.BooleanField()
  isBackward = models.BooleanField()
  isProkaryotic = models.BooleanField()
  specificOrganism = models.CharField(max_length=100)
  biobrick_id_num = models.PositiveIntegerField(null=True)
  biobrick_status = models.CharField(max_length=10, null=True)
  def __unicode__(self):
    return self.name +' '+ self.description

class Promoter(models.Model):
  sequence = models.CharField(max_length=2500)
  name = models.CharField(max_length=25, unique=True)
  isBioBrick = models.BooleanField()
  description = models.CharField(max_length=200)
  isConstitutive = models.BooleanField()
  isProkaryotic = models.BooleanField()
  specificOrganism = models.CharField(max_length=100)
  biobrick_id_num = models.PositiveIntegerField(null=True)
  biobrick_status = models.CharField(max_length=10, null=True)
  def __unicode__(self):
    return self.name +' '+ self.description

class RestrictionSites(models.Model):
  RE1 = models.CharField(max_length=25)
  RE2 = models.CharField(max_length=25)
  RE3 = models.CharField(max_length=25)
  RE4 = models.CharField(max_length=25)
  RE5 = models.CharField(max_length=25)
  def __unicode__(self):
    return self.name

class ResistanceMarker(models.Model):
  name = models.CharField(max_length=25)
  def __unicode__(self):
    return self.name

class Plasmid(models.Model):
  sequence = models.CharField(max_length=20000)
  name = models.CharField(max_length=25, unique=True)
  isBioBrick = models.BooleanField()
  copynum = models.CharField(max_length=25, null=True)
  description = models.CharField(max_length=200)
  biobrick_id_num = models.PositiveIntegerField(null=True)
  biobrick_status = models.CharField(max_length=10, null=True)
  biobrick_plasmid_oriNum = models.IntegerField(null=True)
  biobrick_plasmid_versNum = models.IntegerField(null=True)
  biobrick_A_res = models.NullBooleanField(null=True)   # ampicillin
  biobrick_C_res = models.NullBooleanField(null=True)   # chloramphenicol
  biobrick_E_res = models.NullBooleanField(null=True)   # erythromycin
  biobrick_G_res = models.NullBooleanField(null=True)   # gentamycin
  biobrick_K_res = models.NullBooleanField(null=True)   # kanamycin
  biobrick_N_res = models.NullBooleanField(null=True)   # neomycin
  biobrick_Na_res = models.NullBooleanField(null=True)  # nalidixic acid
  biobrick_R_res = models.NullBooleanField(null=True)   # rifampicin
  biobrick_S_res = models.NullBooleanField(null=True)   # spectinomycin
  biobrick_St_res = models.NullBooleanField(null=True)  # streptomycin
  biobrick_T_res = models.NullBooleanField(null=True)   # tetracycline
  biobrick_Tm_res = models.NullBooleanField(null=True)  # trimethoprim
  biobrick_Z_res = models.NullBooleanField(null=True)   # zeocin
  #oriStart = models.IntegerField()
  #oriEnd= models.IntegerField()
  resistance_marker = models.ManyToManyField(ResistanceMarker, through='resLocation')
  multiple_cloning_site = models.ManyToManyField(RestrictionSites, through='mcsLocation')
  def get_general_plasmid_purpose(self):
    return self.get_ori_purpose() + '; ' + self.get_versNum_purpose()
  def get_general_plasmid_description(self):
    return self.get_ori_description() + '; ' + self.get_versNum_description()
  def get_versNum_description(self):
    if (self.isBioBrick is False) or (self.biobrick_plasmid_versNum is None):
      return ''
    else:
      if self.biobrick_plasmid_versNum == 0:
        return 'BioBrick cloning site absent or incomplete'
      elif self.biobrick_plasmid_versNum == 1:
        return 'complete BioBrick cloning site (BCS)'
      elif self.biobrick_plasmid_versNum == 2:
        return '5\' terminator and BCS'
      elif self.biobrick_plasmid_versNum == 3:
        return '5\' terminator and BCS and 3\' terminator'
      elif self.biobrick_plasmid_versNum == 4:
        return 'pSB2K3-derived plasmid backbone free of many restriction sites'
      elif self.biobrick_plasmid_versNum == 5:
        return 'constructed from BioBrick base vector'
      elif self.biobrick_plasmid_versNum == 7:
        return 'BCS flanked by terminator BBa_B0015'
      elif self.biobrick_plasmid_versNum == 8:
        return 'BCS followed by eukaryotic terminator'
      elif self.biobrick_plasmid_versNum == 10:
        return 'Endy screening plasmid v1.0'
    return ''
  def get_versNum_purpose(self):
    if (self.isBioBrick is False) or (self.biobrick_plasmid_versNum is None):
      return ''
    else:
      if self.biobrick_plasmid_versNum == 0:
        return 'no BioBrick cloning site'
      elif self.biobrick_plasmid_versNum == 1:
        return 'assembly of BioBrick parts'
      elif self.biobrick_plasmid_versNum == 2:
        return 'transcriptional insulation of plasmid backbone upstream of cloned BioBrick part'
      elif self.biobrick_plasmid_versNum == 3:
        return 'transcriptional insulation of plasmid backbone downstream of cloned BioBrick part'
      elif self.biobrick_plasmid_versNum == 4:
        return 'Genome refactoring'
      elif self.biobrick_plasmid_versNum == 5:
        return 'standardized BioBrick vector design'
      elif self.biobrick_plasmid_versNum == 7:
        return 'transcriptional insulation of cloned BioBrick part'
      elif self.biobrick_plasmid_versNum == 8:
        return 'for cloning of eukaryotic composite parts'
      elif self.biobrick_plasmid_versNum == 10:
        return 'characterization of single input, single output transcriptional devices'
    return ''
  def get_ori_description(self):
    if (self.isBioBrick is False) or (self.biobrick_plasmid_oriNum is None):
      return ''
    else:
      if self.biobrick_plasmid_oriNum == 1:
        return 'modified pMB1 derived from pUC19'
      elif self.biobrick_plasmid_oriNum == 2:
        return 'F and P1 lytic derived from pSCANS-1-BNL'
      elif self.biobrick_plasmid_oriNum == 3:
        return 'p15A derived from pMR101'
      elif self.biobrick_plasmid_oriNum == 4:
        return 'rep101, repA derived from pSC101'
      elif self.biobrick_plasmid_oriNum == 5:
        return 'derived from F plasmid'
      elif self.biobrick_plasmid_oriNum == 6:
        return 'pMB1 derived from pBR322'
    return ''
  def get_ori_purpose(self):
    if (self.isBioBrick is False) or (self.biobrick_plasmid_oriNum is None):
      return ''
    else:
      if self.biobrick_plasmid_oriNum == 1:
        return 'Easy plasmid DNA purification'
      elif self.biobrick_plasmid_oriNum == 2:
        return 'Inducible copy number'
      elif self.biobrick_plasmid_oriNum == 3:
        return 'Multi-plasmid engineered systems'
      elif self.biobrick_plasmid_oriNum == 4:
        return 'Small cell to cell copy number variation'
      elif self.biobrick_plasmid_oriNum == 5:
        return 'Improved plasmid stability'
      elif self.biobrick_plasmid_oriNum == 6:
        return 'Multi-plasmid engineered systems'
    return ''
  def get_copy_num(self):
    if (self.isBioBrick is False) or (self.biobrick_plasmid_oriNum is None):
      return ''
    else:
      if self.biobrick_plasmid_oriNum == 1:
        return '500-700'
      elif self.biobrick_plasmid_oriNum == 2:
        return '1-2'
      elif self.biobrick_plasmid_oriNum == 3:
        return '10-12'
      elif self.biobrick_plasmid_oriNum == 4:
        return '4-6'
      elif self.biobrick_plasmid_oriNum == 5:
        return '1-2'
      elif self.biobrick_plasmid_oriNum == 6:
        return '15-20'
    return ''
  def get_antibiotic_resistances(self):
    resistances = []
    if (self.isBioBrick is False) or (self.biobrick_plasmid_oriNum is None):
      return ''
    else:
      if self.biobrick_A_res:
        resistances.append('ampicillin')
      if self.biobrick_C_res:
        resistances.append('chloramphenicol')
      if self.biobrick_E_res:
        resistances.append('erythromycin')
      if self.biobrick_G_res:
        resistances.append('gentamycin')
      if self.biobrick_K_res:
        resistances.append('kanamycin')
      if self.biobrick_N_res:
        resistances.append('neomycin')
      if self.biobrick_Na_res:
        resistances.append('nalidixic acid')
      if self.biobrick_R_res:
        resistances.append('rifampicin')
      if self.biobrick_S_res:
        resistances.append('spectinomycin')
      if self.biobrick_St_res:
        resistances.append('streptomycin')
      if self.biobrick_T_res:
        resistances.append('tetracycline')
      if self.biobrick_Tm_res:
        resistances.append('trimethoprim')
      if self.biobrick_Z_res:
        resistances.append('zeocin')
    return resistances
  def __unicode__(self):
    return self.description
    #length = len(self.description)
    #if length == 0:
    #  return self.seqrec.name + ': no description available' 
    #elif length > 35:
    #  return self.seqrec.name + ': ' + self.seqrec.description[:35] + '...'
    #else:
    #  return self.seqrec.name + ': ' + self.seqrec.description[:35] 
  class Meta:
    ordering = ['isBioBrick']

  #-TODO------------------------------
  def convert_to_seqrec(self):
    pass
  def __get_resistance_marker_features__(self, start, end, name, strand=None):
    startpos = ExactPosition(start)
    endpos = ExactPosition(end)
    loc = FeatureLocation(startpos, endpos)
    return SeqFeature(location=loc,type='res',strand=strand,id=name)
  def __get_restriction_site_features__(self, start, end, name, strand=None):
    startpos = ExactPosition(start)
    endpos = ExactPosition(end)
    loc = FeatureLocation(startpos, endpos)
    return SeqFeature(location=loc,type='mcs',strand=strand,id=name)
  #-----------------------------------

class mcsLocation(models.Model):
  plasmid = models.ForeignKey(Plasmid)
  restriction_sites = models.ForeignKey(RestrictionSites)
  startloc1 = models.IntegerField()
  endloc1 = models.IntegerField()
  startloc2 = models.IntegerField()
  endloc2 = models.IntegerField()
  startloc3 = models.IntegerField()
  endloc3 = models.IntegerField()
  startloc4 = models.IntegerField()
  endloc4 = models.IntegerField()
  startloc5 = models.IntegerField()
  endloc5 = models.IntegerField()

class resLocation(models.Model):
  plasmid = models.ForeignKey(Plasmid)
  resistance = models.ForeignKey(ResistanceMarker)
  startloc = models.IntegerField()
  endloc = models.IntegerField()

    
