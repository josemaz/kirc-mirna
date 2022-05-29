import py4cytoscape as p4c

df = pd.read_csv('data/tables/genes-only-stages.tsv',sep='\t')
net = df.loc[:,['gene', 'miR']]
net.columns = ["source","target"]


net  = pd.read_csv('data/tables/mi/ctrl-all-genmirna-100k.txt',sep='\t')
net = net.loc[net.Source.isin(df.gene),:]
p4c.create_network_from_data_frames(edges = net, title='From nodes')


phe = ["ctrl","stage1"]

for i in phe:
    print(i)
p4c.import_network_from_tabular_file(
        file='data/tables/mi/' + i + '-all-genmirna-100k.txt',
        first_row_as_column_names = True, column_type_list='s,t,i')

