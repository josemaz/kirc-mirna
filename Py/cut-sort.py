import pandas as pd
import sys, utils

tipos = ["ctrl","stagei","stageii","stageiii","stageiv"]

nl = int(sys.argv[1])
for t in tipos:
	fin = 'Output/peaval-'+ t +'.tsv'
	utils.logprint(fin)
	df = pd.read_csv(fin, sep = "\t")
	df['P'] = df['P'].abs()
	df = df.sort_values('P', ascending=False)
	fout = fin.split('.')[0] + '-' + \
		utils.sizeof_fmt(nl) + '.tsv'
	df.head(nl).to_csv(fout, 
		index = False, header=True, sep='\t')
	# print(df.tail(10))
