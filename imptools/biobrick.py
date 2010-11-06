
def get_sequence(gene): 
    import urllib
    url='http://www.kegg.jp/dbget-bin/www_bget?-f+-n+n+'
    url+=gene
    f=urllib.urlopen(url).read()
    sequence=f[f.find('(N)')+4:f.find('</pre>')].replace('\n','')
    return sequence.upper()
    
    
def check_sequence(sequence):
    colors=[]
    EcoRIone=sequence.find('GAATTC')  ##yellow
    XbaIone=sequence.find('TCTAGA')     #Chartreuse
    SpecIone=sequence.find('ACTAGT')        #Aqua
    PstIone=sequence.find('CTGCAG')         ##DeepPink
    nothingFound=True
    while EcoRIone>=0:
        nothingFound=False
        #print "GAATTC found at", (EcoRIone,EcoRIone+5)
        colors.append((EcoRIone,'Yellow'))
        EcoRIone=sequence.find('GAATTC',EcoRIone+1)
    while XbaIone>=0:
        nothingFound=False
        #print "TCTAGA found at", (XbaIone,XbaIone+5)
        colors.append((XbaIone,'Chartreuse'))
        XbaIone=sequence.find('TCTAGA',XbaIone+1)
    while SpecIone>=0:
        nothingFound=False
        #print "ACTAGT found at", (SpecIone,SpecIone+5)
        colors.append((SpecIone,'Aqua'))
        SpecIone=sequence.find('ACTAGT',SpecIone+1)
    while PstIone>=0:
        nothingFound=False
        #print "CTGCAG found at", (PstIone,PstIone+5)
        colors.append((PstIone,'DeepPink'))
        PstIone=sequence.find('CTGCAG',PstIone+1)
    if nothingFound:
        return False
    else:
        return colors
        
def add_color(seq,colors):
    seq=add_breaks(seq)
    lowest=0
    new=''
    colors.sort()
    while True:
        try:
            (next,color)=colors.pop(0)
            new = new + seq[lowest:next+5*(next/50)]+'<FONT style="BACKGROUND-COLOR:'+color+'">'+seq[next+5*(next/50):next+6+5*((next+6)/50)]+'</FONT>'
            lowest=next+6+5*((next+6)/50)
        except:
            new+=seq[lowest:]
            break
    return __format_for_django__(new)
    
def add_breaks(seq):
    new=''
    x=50
    while x<len(seq):
        new+=seq[x-50:x]+'<br/>'
        x+=50
    new+=seq[x-50:]
    return new

def __format_for_django__(seq):
  retval = []
  seq = seq.split('<')
  for token in seq:
    if token.startswith('br'):
      retval.append('<br/>')
      retval.append(token[4:])
    elif token.__contains__('FONT'):
      #TODO
      pass
