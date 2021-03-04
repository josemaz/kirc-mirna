import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from termcolor import colored, cprint

logprint = lambda x: cprint(x, 'red', attrs=["bold"])
warprint = lambda x: cprint(x, 'green', attrs=["bold"])

def sizeof_fmt(num, suffix=''):
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(num) < 1000.0:
            return "%.f%s%s" % (num, unit, suffix)
        num /= 1000.0
    return "%.1f %s%s" % (num, 'Y', suffix)

def plotcor(datos, lab, odir):

	Path(odir).mkdir(parents=True, exist_ok=True)

	print("Ploting ...")
	plt.figure(figsize=(10,10),dpi=300)
	masking = np.tril(datos)
	sns.heatmap(datos ,xticklabels=False, yticklabels=False, 
		mask=masking, vmin=-1., vmax=1., square=True,
		# cbar_kws={"shrink": 0.75}, cmap=RdBu_11.mpl_colormap)
		cbar_kws={"shrink": 0.75}, cmap=plt.get_cmap('seismic'))
	plt.title('Pearson correlation ' + lab)
	plt.ylabel('Gene start position')
	plt.xlabel('Gene start position')
	plt.tight_layout()
	plt.savefig(odir + '/' + lab + '.png')
	plt.clf()
	plt.close() # to clean memory