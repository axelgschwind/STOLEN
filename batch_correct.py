from combat.pycombat import pycombat
import pandas as pd
import matplotlib.pyplot as plt
from optparse import OptionParser
import sys
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import os
import numpy as np

def main(argv):
	if len(argv) == 1:
		argv.append("--help")
	usage = "usage: %prog -i counts.tsv -b metadata.tsv -o batch_corrected.tsv"
	desc = "Batch correction of gene expression. Input file must contain a matrix where the first column is labeled GENE and any other column represents a sample identifier. The metadata file needs a column RNA and COHORT."
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--in", action="store", dest="input", type="string", help="Input gene expression matrix. First column must be labeled GENE and any other column as a sample identifier.")
	parser.add_option("-o", "--out", action="store", dest="out", type="string", help="Output batch corrected TSV")
	parser.add_option("-b", "--batch", action="store", dest="batch", type="string", default="tpm", help="TSV file with batch information. Must contain a column RNA and a column COHORT.")
	parser.add_option("--plot", action="store_true", dest="plot", help="Plot a PCA analys with both corrected and uncorrected")
	(options, args) = parser.parse_args()
 
	#output directory
	odir = os.path.dirname(os.path.abspath(options.out))
 
	df_expression = pd.read_csv(options.input, delimiter="\t", index_col=0)
	df_expression.index.name = "GENE"

	dat = df_expression

	#Remove rows containing only zeroes (pyCombat cannot parse them)
	dat = df_expression.loc[(df_expression != 0).any(axis=1)]
	zero_values = df_expression.loc[~(df_expression != 0).any(axis=1)]

	metadata = pd.read_csv(options.batch, delimiter="\t")
	batch = []
	for col_name in dat.columns:
		metadata[metadata["RNA"] == col_name]
		batch.append(metadata.loc[metadata["RNA"] == col_name,"COHORT"].item())

	corrected = pycombat(dat, batch)
 
	corrected = pd.concat([corrected, zero_values])

	corrected.to_csv(options.out, sep="\t")


	if options.plot:
		pca = PCA()
		pca.fit(corrected.T)
		feature = pca.transform(corrected.T)
		#PLOT PCA of corrected data
		plt.figure(figsize=(12, 12))

		#encode colors for every cohort batch
		dictionary = {}
		for idx, entry in enumerate(set(batch)):
			dictionary[entry] = idx
		colors = [dictionary.get(e,e) for e in batch]
		
  

		plt.scatter(feature[:, 0], feature[:, 1], alpha=0.8, c=list(colors), cmap="jet")
		plt.grid()
		plt.xlabel("PC1")
		plt.ylabel("PC2")
		plt.savefig(odir + "/" + "pca_corrected.pdf")
	
		#PCA uncorrected data
		pca = PCA()
		df_expression = df_expression.transform(lambda x: np.log10(x+1))
		pca.fit(df_expression.T)
		feature = pca.transform(df_expression.T)

		plt.figure(figsize=(12, 12))
		plt.scatter(feature[:, 0], feature[:, 1], alpha=0.8, c=list(colors), cmap="jet")
		plt.grid()
		plt.xlabel("PC1")
		plt.ylabel("PC2")
		plt.savefig(odir + "/" + "pca_uncorrected.pdf")

if __name__ == "__main__":
	main(sys.argv)