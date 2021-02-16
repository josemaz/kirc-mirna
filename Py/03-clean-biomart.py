import pandas as pd
import re
from termcolor import colored, cprint

logprint = lambda x: cprint(x, 'red', attrs=["bold"])
# pd.set_option('display.max_columns', None)



#######################################################################
logprint("Loding and cleaning Biomart ...")
mart = pd.read_csv("Data/biomart-20210119.txt", sep = "\t")
colnames = {
	'Gene stable ID': 'GeneID', 
	'Karyotype band': 'KaryoBand',
	'Gene start (bp)': 'gStart',
	'Gene end (bp)': 'gEnd',
	'Gene % GC content': 'GC',
	'Gene name': 'gName',
	'Chromosome/scaffold name': 'Chrom',
	'miRBase ID': 'miRBaseID',
	'Transcript type': 'TransType',
}
mart = mart.rename(columns=colnames)
mart.Chrom = mart.Chrom.astype(str)
chrs = [str(i) for i in range(1,23)]
chrs.append('X')
chrs.append('Y')
mart = mart[mart['Chrom'].isin(chrs)]
# mart[(mart.Chrom == "19") & (mart.TransType == "miRNA" )]
tipo = ['protein_coding','miRNA']
mart = mart[mart['TransType'].isin(tipo)]
mart = mart.drop_duplicates(subset=['GeneID'])
# drops = mart[(mart['TransType']=="miRNA") & (mart['miRBaseID'].isna())].index
# mart = mart.drop(drops)
# def switch(row):
# 	if row['TransType'] == "miRNA":
# 		row['GeneID'] = row['gName']
# 	return row
# mart = mart.apply(lambda row: switch(row), axis=1)
def mods(row):
	if row['TransType'] == "miRNA":
		# Case 1
		if pd.isnull(row.miRBaseID): 
			row['gid'] = row['gName']
		# Case 2
		else: 
			div = row.miRBaseID.split("-")
			row['gid'] = div[1].upper() + div[2].upper()
	# Case 3
	else:
		row['gid'] = row['GeneID']
	return row
mart = mart.apply(lambda row: mods(row), axis=1)
mart = mart.drop_duplicates(subset = ["gid"]) # [20853 rows x 10 columns]



#######################################################################
logprint("Merging expression set and Biomart ...")
def clean(gen):
	if re.match('^hsa',gen):
		div = gen.split("-")
		return div[1].upper() + div[2].upper()
	else:
		return gen
ids = pd.read_csv("Output/ids.tsv", sep = "\t")
for i, group in ids.groupby("tipo"):
	logprint("Group: " + i)
	g = i.replace(' ', '')
	expr = pd.read_csv("Output/paste-miRNA-" + g + ".tsv", sep = "\t")
	expr['gid'] = expr['gene'].apply(lambda gen: clean(gen))
	# print(expr)	
	#! miRNAs hardcoding to merge with Biomart
	expr.loc[expr.gid == "MIR320B", "gid"] = "MIR320B1"
	expr.loc[expr.gid == "MIR320C", "gid"] = "MIR320C1"
	expr.loc[expr.gid == "MIR320D", "gid"] = "MIR320D1"

	merged = pd.merge(expr,mart, on='gid', how='left')
	print(merged[merged.TransType.isna( )][['gene','gid','miRBaseID','TransType']])
	drops = merged[merged.Chrom.isna()].index
	merged = merged.drop(drops)

	# print(merged[merged.TransType.isna( )][['gene','gid','miRBaseID']])

	cols = ['gene','gid'] + list(colnames.values())
	annot = merged[cols]
	cols = ['gid'] + list(colnames.values())
	expr = merged.drop(cols, axis = 1)

	fname = "Output/expr-all-" + g + ".tsv"
	expr.to_csv(fname, index = False, header=True, sep='\t')
	fname = "Output/annot-all-" + g + ".tsv"
	annot.to_csv(fname, index = False, header=True, sep='\t')













#! Filtros
# expr[expr.gene.str.contains('^hsa')]
# merged[merged.TransType == "miRNA"]
# merged[merged.TransType.isna()]
# merged[merged.TransType.isna()][['gene','gid','miRBaseID']]

# mart[mart.TransType=="miRNA"]














