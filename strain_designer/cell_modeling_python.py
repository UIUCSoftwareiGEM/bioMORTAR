####################################################################
# Convert from KEGG ID to BiGG
####################################################################

from imptools.stoich import matrix_for_cobra

####################################################################
# May have to change variable names!!!!!!!!!!!!!!!!!!
####################################################################
def cell_modeling_init(IMPhost_organism, rxns, compounds, initConcentrations, initBiomass):

####################################################################
# Import re
####################################################################
    import re
    matrix, kegg_names, names = matrix_for_cobra(rxns, compounds)
    bigg_names = []

####################################################################
# Check if the organism chosen matches database organisms
####################################################################
    if (re.search('coli',IMPhost_organism)):
        organism = 'Ecoli_iAF1260_flux1'
    elif (re.search('cerevisiae',IMPhost_organism)):
          organism = 'Ecoli_iAF1260_flux1'
    elif (re.search('aureus', IMPhost_organism)):
          organism = 'Saureus_iSB619'
    elif (re.search('tuberculosis',IMPhost_organism)):
          organism = 'Mtuberculosis_iNJ661'
    elif (re.search('barkeri',IMPhost_organism)):
          organism = 'Mbarkeri_iAF692_flux1'
    elif (re.search('pylori',IMPhost_organism)):
          organism = 'Hpylori_iTT341_min_medI'
    else:
          organism = ''

####################################################################
# Begin loop function to change reaction from KeggID to BiGG ID
####################################################################
    if (IMPhost_organism != ''):
        for i in range(len(matrix)):
            for line in open(str(IMPhost_organism)+'.xml',r):
                for name in names[i]:
                    if name in line:
                        if (re.search(r'id=\"\w+\"',line)):
                            idn = re.search(r'id=\"\w+\"',line).group(0)
                            idn = idn[4:-2]
                            if idn.endswith("_c"):
                                idn = idn[:-2] + "[c]"
                            if idn.endswith("_e"):
                                idn = idn[:-2] + "[e]"
                            idn = idn[2:]
                        else:
                            idn = name + ""

                        bigg_names.append(idn)

####################################################################
# Open MATLAB
####################################################################
    from mlabwrap import mlab;
    import numpy;
    # Start MATLAB engine session - by default it opens on current host
    # Does it need 'matrix' after organism and before bigg_names???
    modelinggraph = mlab.cell_modeling(organism,bigg_names,initConcentrations,initBiomass)

####################################################################

org = raw_input()
rxns = raw_input()
comps = raw_input()
initC = raw_input()
initB = raw_input()
cell_modeling_init(org,rxns,comps,initC,initB)

