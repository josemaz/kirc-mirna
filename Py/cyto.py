# conda install -c conda-forge py2cytoscape
import glob, sys, os
import pandas as pd
import numpy as np
from py2cytoscape.data.cyrest_client import CyRestClient
import utils

print("Starting Ch3m Cytoscape...")

# cy = CyRestClient() - This default constructor creates connection to http://localhost:1234/v1
cy = CyRestClient(ip='127.0.0.1', port=1234)

# Cleanup: Delete all existing networks and tables in current Cytoscape session
cy.session.delete()

tipos = ["ctrl","stagei","stageii","stageiii","stageiv"]
nets = []
for t in tipos:
	fin = 'Output/Pearson/peaval-' + t + '-10K.sif'
	if os.path.exists(fin):
		utils.logprint(f'Opening: {fin}')
		net = cy.network.create_from(fin, collection="Pearson 10K")
		nets.append(net)
	else:
		utils.warprint(f'Error opening: {fin}')
		sys.exit(15)
	fin = 'Output/MI/expr-all-' + t + '-genmirna-10K.sif'
	if os.path.exists(fin):
		utils.logprint(f'Opening: {fin}')
		net = cy.network.create_from(fin, collection="MI 10K")
		nets.append(net)
	else:
		utils.warprint(f'Error opening: {fin}')
		sys.exit(15)
	


td = pd.read_csv('Output/annot-all-ctrl.tsv', sep='\t')
print(td)
for net in nets:
    net.update_node_table(df=td, data_key_col='gene')

# get the node table
node_table = nets[0].get_node_table()
colors = [ ','.join(map(str, np.random.randint(256, size=3))) for i in node_table.Chrom.unique()]
mapping = dict(zip(node_table.chromname.unique(),colors))
#mapping = {
#            '1': '202,75,78',
#            '2': '0,75,78'
#        }
print("Styling networks ...")
my_style = cy.style.create('GAL Style')
my_style.create_discrete_mapping(column='chromname', vp='NODE_FILL_COLOR', col_type='String', mappings=mapping )
for net in nets:
    cy.style.apply(my_style, net)







#######################################################################
# 
# for f in files:
#     collec = f.split('/')[-1]
#     collec = collec.split('.')[0]
#     collec = collec.split('-')[0] + '-' + collec.split('-')[2]
#     net = cy.network.create_from(f, collection=collec)
#     nets.append(net)
