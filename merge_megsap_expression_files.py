# Author: Axel Gschwind
# program: merge_megsap_expression_files.py
# purpose: Merges gene counts files from megSAP RNA pipeline into one tsv file

import sys
import pandas as pd
from optparse import OptionParser
import os.path


def main(argv):
	if len(argv) == 1:
		argv.append("--help")
	usage = "usage: %prog -o merged_counts_file.tsv INPUT_1.tsv INPUT_2.tsv INPUT_N.tsv"
	desc = "Merges gene counts files from megSAP RNA pipeline into one tsv file."
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-o", "--out", action="store", dest="out", type="string", help="Output merged TSV")
	parser.add_option("--column", action="store", dest="column", type="string", default="tpm", help="Column from input file to merge")
	(options, args) = parser.parse_args()

	merged_df = pd.DataFrame()
	for file in args:
		if not os.path.isfile(file):
			continue
		
		data = pd.read_csv(file, delimiter="\t", header=0, index_col=0, )
		data.index.name = "ID"

		#get file ID from megSAP file name
		id = os.path.basename(file).replace("_counts.tsv","")
		
		data.rename(columns={options.column: id}, inplace=True)
		data = data[id]
		merged_df = merged_df.merge(data, how="outer", left_index=True, right_index=True)

	merged_df.dropna(axis=0, inplace=True)

	merged_df.to_csv(options.out, sep="\t")

		

if __name__ == "__main__":
	main(sys.argv)