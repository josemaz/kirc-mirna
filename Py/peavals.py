import pandas as pd
import utils
from re import match

tipos = ["ctrl","stagei","stageii","stageiii","stageiv"]
# tipos = ["ctrl","stagei"]

# Stage plot itration
for t in tipos:
	utils.logprint("Type: " + t)
	fin = "Output/expr-all-" + t + ".tsv"
	df = pd.read_csv(fin, sep = "\t")

	# Assuming genes are continues with Ensamble id
	genes = list(filter(lambda v: match('^ENS', v), df.gene))
	ngenes = len(genes)
	df.set_index('gene', inplace=True)

	

	data = df.iloc[:,1:].T.corr()

	mirna_gen = data.iloc[ngenes:,:ngenes]
	mirna_gen.index.name = None
	mirna_gen = mirna_gen.stack().reset_index()
	mirna_gen.columns = ['Source','Target','P']
	mirna_gen = mirna_gen.sort_values('P', ascending=False)

	mirna_gen.to_csv('Output/peaval-' + t + '.tsv', 
		index = False, header=True, sep='\t')

	utils.logprint(f'Type {t} Finished.')



	# fin = "Output/annot-all-" + t + ".tsv"
	# annot = pd.read_csv(fin, sep = "\t")
	# print(df.shape,annot.shape)
	# df = pd.merge(df, annot, on="gene", how="inner")


