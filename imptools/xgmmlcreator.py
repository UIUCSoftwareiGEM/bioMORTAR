def makeXgmml(maincmps,rxns,reactants,products):	
	#maincmps = ['C06990','C12837','C05618','C03572','C12838','C02222','C00846','C00042','C00104','C00031','C00469']
	#rxns = ['R06856','R06857','R06839','R06840','R06838','R02989','R02990','R00727','R00725','R04410']
	## from tools.app1.models import Reaction
	##>>> p = []
	##>>> r = []
	##>>> for rxn in rxns:
	##...     rx = Reaction.objects.get(KeggID = rxn)
	##...     reacts = rx.reactants.all()
	##...     re = []
	##...     for item in reacts:
	##...             item2 = str(item)
	##...             re.append(item2)
	##...     prods = rx.products.all()
	##...     pr = []
	##...     for item in prods:
	##...             item2 = str(item)
	##...             pr.append(item2)
	##...     r.append(re)
	##...     p.append(pr)
	#reactants = [['C00004', 'C00007', 'C00080', 'C06990'], ['C00003', 'C12837'], ['C00007', 'C05618'], ['C03572'], ['C00001', 'C12838'], ['C00006', 'C00846'], ['C00091', 'C00846'], ['C00010', 'C00042', 'C00081'], ['C00031', 'C00081'], ['C00469', 'C04164']]
	#products = [['C00003', 'C12837'], ['C00004', 'C00080', 'C05618'], ['C03572'], ['C01327', 'C12838'], ['C02222'], ['C00005', 'C00080', 'C02222'], ['C00042', 'C02232'], ['C00009', 'C00091', 'C00104'], ['C00092', 'C00104'], ['C00031', 'C06359']]
	reactants2 = []
	
	for lst in reactants:
		newlst = []
		for item in lst:
			try:
				maincmps.index(item)
				continue
			except:
				newlst.append(item)
		reactants2.append(newlst)
	products2 = []
	for lst in products:
		newlst = []
		for item in lst:
			try:
				maincmps.index(item)
			except:
				newlst.append(item)
		products2.append(newlst)
	crazy = []
	
	z = 0
	while z < len(reactants2):
		item2 = reactants2[z]
		crazy.append(item2)
		item = products2[z]
		crazy.append(item)
	
		z += 1
	q = 0
	extracmps = []
	while q < len(crazy)-1:
		extracmps.append(crazy[q:q+2])
		q += 2
	#########################################################################
	graph = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n', '<graph label="Network" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML"  directed="1">\n', '  <att name="documentVersion" value="1.1"/>\n', '  <att name="networkMetadata">\n', '    <rdf:RDF>\n', '      <rdf:Description rdf:about="http://www.cytoscape.org/">\n', '        <dc:type>Protein-Protein Interaction</dc:type>\n', '        <dc:description>N/A</dc:description>\n', '        <dc:identifier>N/A</dc:identifier>\n', '        <dc:title>gmlpractice1.gml</dc:title>\n', '        <dc:source>http://www.cytoscape.org/</dc:source>\n', '        <dc:format>Cytoscape-XGMML</dc:format>\n', '      </rdf:Description>\n', '    </rdf:RDF>\n', '  </att>\n', '  <att type="string" name="backgroundColor" value="#ccccff"/>\n', '  <att type="real" name="GRAPH_VIEW_ZOOM" value="0.35147127673782214"/>\n', '  <att type="real" name="GRAPH_VIEW_CENTER_X" value="1001.4205932617188"/>\n', '  <att type="real" name="GRAPH_VIEW_CENTER_Y" value="-100.0"/>\n']
	root_index = -1
	x = 0.0
	numnodes = 0
	for cmpd in maincmps:
		node = '  <node label="' + cmpd + '" id="' + str(root_index) + '">\n'
		canonicalName = '    <att type="string" name="canonicalName" value="' + cmpd + '"/>\n'
		actualname = '    <att type="string" name="Name" value="INSERTNAMEHERE"/>\n'
		graphics = '    <graphics type="ELLIPSE" h="40.0" w="40.0" x="' + str(x) + '" y="0.0" fill="#ff9999" width="1" outline="#666666" cy:nodeTransparency="1.0" cy:nodeLabelFont="Default-0-12" cy:borderLineType="solid"/>\n'
		newnode = [node, canonicalName, actualname, graphics, '  </node>\n']
		x = x + 200
		root_index = root_index-2
		graph.extend(newnode)
		numnodes += 1
	root_index = -2
	x = 100.0
	for rxn in rxns:
		node = '  <node label="' + rxn + '" id="' + str(root_index) + '">\n'
		canonicalName = '    <att type="string" name="canonicalName" value="' + rxn + '"/>\n'
		actualname = '    <att type="string" name="Enzyme" value="INSERTENZYMEHERE"/>\n'
		graphics = '    <graphics type="RECTANGLE" h="40.0" w="40.0" x="' + str(x) + '" y="0.0" fill="#CD3278" width="1" outline="#666666" cy:nodeTransparency="1.0" cy:nodeLabelFont="Default-0-12" cy:borderLineType="solid"/>\n'
		newnode = [node, canonicalName, actualname, graphics, '  </node>\n']
		x = x + 200
		root_index = root_index-2
		graph.extend(newnode)
		numnodes += 1
	
	numnodes *= -1
	root_index = numnodes-1
	numnodes *= -1
	x = 0.0
	move = 100
	checkcmplist = []
	d = {}
	for lsts in extracmps:
		lst = lsts[0]
		y = 0.0
		move *= -1
		if lst != []:
			for cmpd in lst:
				try:
					checkcmplist.index(cmpd)
					#print cmpd, 'a'
				except:
					y += move
					node = '  <node label="' + cmpd + '" id="' + str(root_index) + '">\n'
					canonicalName = '    <att type="string" name="canonicalName" value="' + cmpd + '"/>\n'
					actualname = '    <att type="string" name="Name" value="INSERTNAMEHERE"/>\n'
					graphics = '    <graphics type="ELLIPSE" h="40.0" w="40.0" x="' + str(x) + '" y="' + str(y) + '" fill="#ffe7ba" width="1" outline="#666666" cy:nodeTransparency="1.0" cy:nodeLabelFont="Default-0-12" cy:borderLineType="solid"/>\n'
					newnode = [node, canonicalName, actualname, graphics, '  </node>\n']
					d.update({cmpd: root_index})
					root_index = root_index-1
					graph.extend(newnode)
					checkcmplist.append(cmpd)
		lst = lsts[1]
		y = 0.0
		move *= -1
		if lst != []:
			for cmpd in lst:
				try:
					checkcmplist.index(cmpd)
					#print cmpd, 'b'
				except:
					y += move
					node = '  <node label="' + cmpd + '" id="' + str(root_index) + '">\n'
					canonicalName = '    <att type="string" name="canonicalName" value="' + cmpd + '"/>\n'
					actualname = '    <att type="string" name="Name" value="INSERTNAMEHERE"/>\n'
					graphics = '    <graphics type="ELLIPSE" h="40.0" w="40.0" x="' + str(x+200) + '" y="' + str(y) + '" fill="#ffe7ba" width="1" outline="#666666" cy:nodeTransparency="1.0" cy:nodeLabelFont="Default-0-12" cy:borderLineType="solid"/>\n'
					newnode = [node, canonicalName, actualname, graphics, '  </node>\n']
					d.update({cmpd: root_index})
					root_index = root_index-1
					graph.extend(newnode)
					checkcmplist.append(cmpd)
		x += 200 
	
	q = 0
	r = 0
	lbl = 0.1
	s = 0
	t = -1
	while r < len(rxns):
		lbl += 1
		s -= 1
		t -= 1
		edge = '  <edge label="' + maincmps[q] + ' (' + str(lbl) + ') ' + rxns[r] + '" source="' + str(s) + '" target="' + str(t) + '">\n'
		#print maincmps[q], rxns[r]
		canonicalname = '    <att type="string" name="canonicalName" value="' + maincmps[q] + ' (' + str(lbl) + ') ' + rxns[r] + '"/>\n'
		intvalue = '    <att type="string" name="interaction" value="' + str(lbl) + '"/>\n'
		graphics = '    <graphics width="1" fill="#000000" cy:sourceArrow="0" cy:targetArrow="3" cy:sourceArrowColor="#000000" cy:targetArrowColor="#000000" cy:edgeLabelFont="SanSerif-0-10" cy:edgeLineType="SOLID" cy:curved="STRAIGHT_LINES"/>\n'
		newedge = [edge, canonicalname, intvalue, graphics, '  </edge>\n']
		graph.extend(newedge)
		lbl+= 1
		s -= 1
		t -= 1
		edge = '  <edge label="' + rxns[r] + ' (' + str(lbl) + ') ' + maincmps[q+1] + '" source="' + str(s) + '" target="' + str(t) + '">\n'
		#print rxns[r], maincmps[q+1]
		canonicalname = '    <att type="string" name="canonicalName" value="' + rxns[r] + ' (' + str(lbl) + ') ' + maincmps[q+1] + '"/>\n'
		intvalue = '    <att type="string" name="interaction" value="' + str(lbl) + '"/>\n'
		graphics = '    <graphics width="1" fill="#000000" cy:sourceArrow="0" cy:targetArrow="3" cy:sourceArrowColor="#000000" cy:targetArrowColor="#000000" cy:edgeLabelFont="SanSerif-0-10" cy:edgeLineType="SOLID" cy:curved="STRAIGHT_LINES"/>\n'
		newedge = [edge, canonicalname, intvalue, graphics, '  </edge>\n']
		graph.extend(newedge)
		q += 1
		r += 1
	lbln = -1
	t = 0
	for lst in reactants2:
		lbld = 0.2
		t -= 2
		lbln += 2
		for item in lst:
			lbl = lbln + lbld
			lbld += 0.1
			s = d[item]
			r = reactants2.index(lst)
			edge = '  <edge label="' + item + ' (' + str(lbl) + ') ' + rxns[r] + '" source="' + str(s) + '" target="' + str(t) + '">\n'
			canonicalname = '    <att type="string" name="canonicalName" value="' + item + ' (' + str(lbl) + ') ' + rxns[r] + '"/>\n'
			intvalue = '    <att type="string" name="interaction" value="' + str(lbl) + '"/>\n'
			graphics = '    <graphics width="1" fill="#ff0000" cy:sourceArrow="0" cy:targetArrow="3" cy:sourceArrowColor="#000000" cy:targetArrowColor="#ff0000" cy:edgeLabelFont="SanSerif-0-10" cy:edgeLineType="SOLID" cy:curved="STRAIGHT_LINES"/>\n'
			newedge = [edge, canonicalname, intvalue, graphics, '  </edge>\n']
			graph.extend(newedge)
	
	lbln = 0
	s = 0
	for lst in products2:
		lbld = 0.2
		s -= 2
		lbln += 2
		for item in lst:
			lbl = lbln + lbld
			lbld += 0.1
			t = d[item]
			r = products2.index(lst)
			edge = '  <edge label="' + rxns[r] + ' (' + str(lbl) + ') ' + item + '" source="' + str(s) + '" target="' + str(t) + '">\n'
			canonicalname = '    <att type="string" name="canonicalName" value="' + rxns[r] + ' (' + str(lbl) + ') ' + item + '"/>\n'
			intvalue = '    <att type="string" name="interaction" value="' + str(lbl) + '"/>\n'
			graphics = '    <graphics width="1" fill="#0000ff" cy:sourceArrow="0" cy:targetArrow="3" cy:sourceArrowColor="#000000" cy:targetArrowColor="#0000ff" cy:edgeLabelFont="SanSerif-0-10" cy:edgeLineType="SOLID" cy:curved="STRAIGHT_LINES"/>\n'
			newedge = [edge, canonicalname, intvalue, graphics, '  </edge>\n']
			graph.extend(newedge)
	
		
	graph.append('</graph>\n')
	
	#wow = open("xgmmlpractice2.xgmml","w")
	#wow.writelines(graph)
	#wow.close()
	  
	return graph
