import pandas as pd
import utils, sys
import numpy as np
from scipy.stats import pearsonr
# import pingouin as pg
import time
from nancorrmp import NaNCorrMp
import statsmodels.stats.multitest as smt

param1 = "49"

tipos = ["ctrl","stagei","stageii","stageiii","stageiv"]
# tipos = ["ctrl","stagei"]
odir = f"Output/{param1}-Plots/"

# Stage plot itration
for t in tipos:
	utils.logprint("Tipo: " + t)
	# fin = f"Output/expr-all-{t}.tsv"
	fin = f"Output/expr-all-{t}-{param1}.tsv"
	utils.logprint(f"Expression File: {fin}")
	df = pd.read_csv(fin, sep = "\t")
	fin = "Output/annot-all-" + t + ".tsv"
	annot = pd.read_csv(fin, sep = "\t")
	print(df.shape,annot.shape)
	df = pd.merge(df, annot, on="gene", how="inner")
	df.Chrom.replace(['X', 'Y'], [23, 24], inplace=True)
	df.Chrom = df.Chrom.apply(pd.to_numeric)
	df = df.sort_values(["Chrom", "gStart"], ascending = True)
	# print(df.gStart.head(20))
	hlines = np.cumsum(df.groupby('Chrom').size().values)

	ncols = df.shape[1]-10
	dat = df.iloc[:,1:ncols].T

	utils.logprint("All genome correlation")
	corr, p_value = NaNCorrMp.calculate_with_p_value(dat)

	utils.logprint("Pvalues corrections")
	npflat = p_value.values.flatten() # 2d to 1d
	np_pvals = smt.multipletests(npflat, method="bonferroni")
	p_value = pd.DataFrame(
		np_pvals[1].reshape(p_value.shape[0],p_value.shape[1]))

	corr = corr[p_value < 0.05].replace(np.nan, 0)
			
	utils.plotcor(corr, "bonf_pval_01_all_" + t, odir, hlines)

	# break

	# # Chromosome plot itration
	# df.Chrom= df.Chrom.astype(str)
	# for gr_name, df_chr in df.groupby('Chrom'):
	# 	utils.logprint("Working in chr: " + gr_name)
	# 	d = df_chr.sort_values('gStart')
	# 	cor = d.iloc[:,1:ncols].T.corr()
	# 	print("Ploting chr: ", gr_name)
	# 	od = odir + t 
	# 	utils.plotcor(cor, "chr-" + str(gr_name), od)





########################################

	# OPCION 1
	# cor = df.iloc[:,1:ncols].T.corr()


	# OPCION 2
	# cor = pd.DataFrame([ [pearsonr(dat[c], dat[y])[0] for y in dat.columns] \
	# 	for c in dat.columns])
	# print(cor)
	
	# OPCION 3
	# cor = dat.rcorr(stars=False, decimals=4)
	# print(cor)
	# for c in dat.columns:		
	# 	ts = time.time()
	# 	for y in dat.columns:
	# 		pval = pearsonr(dat[c], dat[y])[1]
	# 	print(f"{c}: {time.time()-ts}")