import networkx as nx
import pandas as pd
from venn import venn


data = {
    'ctrl' : None, 
    'stagei' : None, 
    'stageii' : None, 
    'stageiii': None, 
    'stageiv': None
}


for key in data:
    fname = '../Output/MI/expr-all-' + key + '-1e4-genmirna.tsv'
    data[key] = pd.read_csv(fname, sep='\t')


redes = dict.fromkeys(data)
in_degrees = dict.fromkeys(data)
for key in redes:
    redes[key] = nx.from_pandas_edgelist(df=data[key], source='Source', target='Target', 
                                         edge_attr='MI', create_using=nx.DiGraph(),)
    lista = list(redes[key].in_degree())
    in_degrees[key] = sorted(lista, key=lambda tup: tup[1], reverse=True)[:10]
in_degrees_vals = {}
for key in in_degrees:
    in_degrees_vals[key] = [ x[1] for x in in_degrees[key] ]
df = pd.DataFrame.from_dict(in_degrees_vals)
# print(df.to_latex(index=False))
print(df)