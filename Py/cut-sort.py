import pandas as pd
import glob, os, sys, utils
import click

@click.command()
@click.option('--head', default=10000, required=True, help='Header of interactions')
@click.option('--cor',  required=True, type=click.Choice(['Pearson', 'MI'], case_sensitive=False))

def prog(head, cor):
    """Simple program that greets NAME for a total of COUNT times."""    
    if cor == "Pearson":
    	pearson(head)
    else:
    	mi(head)

def pearson(nl):
	tipos = ["ctrl","stagei","stageii","stageiii","stageiv"]
	for t in tipos:
		fin = 'Output/Pearson/peaval-'+ t +'.tsv'
		utils.logprint(fin)
		df = pd.read_csv(fin, sep = "\t")
		df['Pabs'] = df['P'].abs()
		df = df.sort_values('Pabs', ascending=False)
		df = df[df.columns[[0,2,1]]]
		fout = fin.split('.')[0] + '-' + \
			utils.sizeof_fmt(nl) + '.sif'
		df.head(nl).to_csv(fout, 
			index = False, header=False, sep='\t')
		# print(df.tail(10))

def mi(nl):
	os.chdir('Output/MI')
	for f in glob.glob('*-genmirna.tsv'):
		utils.logprint(f)
		df = pd.read_csv(f, sep = "\t")
		df = df.sort_values('MI', ascending=False)
		df = df[df.columns[[0,2,1]]]	
		fout = f.split('.')[0] + '-' + \
			utils.sizeof_fmt(nl) + '.sif'
		df.head(nl).to_csv(fout, 
			index = False, header=False, sep='\t')


if __name__ == '__main__':
    prog()


