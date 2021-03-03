import pandas as pd
import utils

tipos = ["ctrl","stagei","stageii","stageiii","stageiv"]
# tipos = ["ctrl","stagei"]

# Stage plot itration
for t in tipos:
	utils.logprint("Tipo: " + t)
	fin = "Output/expr-all-" + t + ".tsv"
	df = pd.read_csv(fin, sep = "\t")
	fin = "Output/annot-all-" + t + ".tsv"
	annot = pd.read_csv(fin, sep = "\t")
	print(df.shape,annot.shape)
	df = pd.merge(df, annot, on="gene", how="inner")
	df.Chrom.replace(['X', 'Y'], [23, 24], inplace=True)
	df.Chrom = df.Chrom.apply(pd.to_numeric)
	df = df.sort_values(["Chrom", "gStart"], ascending = True)
	# print(df)

	ncols = df.shape[1]-10
	cor = df.iloc[:,1:ncols].T.corr()

	utils.plotcor(cor, "all_" + t, "Output/Plots/")

	# Chromosome plot itration
	df.Chrom= df.Chrom.astype(str)
	for gr_name, df_chr in df.groupby('Chrom'):
		utils.logprint("Working in chr: " + gr_name)
		d = df_chr.sort_values('gStart')
		cor = d.iloc[:,1:ncols].T.corr()
		print("Ploting chr: ", gr_name)
		od = "Output/Plots/" + t 
		utils.plotcor(cor, "chr-" + str(gr_name), od)
