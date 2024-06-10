# Author: Axel Gschwind
# program: calc_expression_profile.py
# purpose: Calculate gene expression profiles

from optparse import OptionParser
import sys
import pandas as pd 
from scipy import stats

def main(argv):
	if len(argv) == 1:
		argv.append("--help")
	usage = "usage: %prog -i input_matrix.tsv -o output_list.tsv"
	desc = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--input", action="store", dest="input", type="string", help="Input gene expression matrix. First column Ensembl gene id. Any other column sample expression data.")
	parser.add_option("-o", "--output", action="store", dest="output", type="string", default="-", help="Output list")
	(options, args) = parser.parse_args()

	expr_data = pd.read_csv(options.input, delimiter="\t", index_col=0)

 
	signatures = {
		#IFN 6-genes, HLA-DRA, CXCL10, CXCL9, IDO1, STAT1, IFNG
		"GEP_IFNG6": ["ENSG00000204287", "ENSG00000169245", "ENSG00000138755", "ENSG00000131203", "ENSG00000115415", "ENSG00000111537"], 
		#preliminary 10 IFNG genes IFNG, STAT1, CCR5, CXCL9, CXCL10, CXCL11, IDO1, PRF1, GZMA, HLA-DRA
		"GEP_IFNG10": ["ENSG00000111537", "ENSG00000115415", "ENSG00000160791", "ENSG00000138755", "ENSG00000169245", "ENSG00000169248", "ENSG00000131203", "ENSG00000180644", "ENSG00000145649", "ENSG00000204287"], 
		#expanded immune gene signature: CD3D, IL2RG, IDO1, NKG7, CIITA, HLA-E, CD3E, CXCR6, CCL5, LAG3, GZMK, TAGAP, CD2, CXCL10, HLA-DRA, STAT1, CXCL13, GZMB
		"GEP_EXPANDED_IMMUNE": ["ENSG00000271503", "ENSG00000116824", "ENSG00000167286", "ENSG00000198851", "ENSG00000179583", "ENSG00000169245", "ENSG00000156234", "ENSG00000172215", "ENSG00000100453", "ENSG00000113088", 
					"ENSG00000204287", "ENSG00000204592", "ENSG00000131203", "ENSG00000147168", "ENSG00000089692", "ENSG00000105374", "ENSG00000115415", "ENSG00000164691"],
		#Immunsuppression signature
		"GEP_IMS": ["ENSG00000148848", "ENSG00000167601", "ENSG00000060982", "ENSG00000181374", "ENSG00000108691", "ENSG00000108700", "ENSG00000177575", "ENSG00000163359", 
					"ENSG00000078098", "ENSG00000136634", "ENSG00000122641", "ENSG00000187608", "ENSG00000162745", "ENSG00000113721", "ENSG00000088827", "ENSG00000159167", "ENSG00000233608", "ENSG00000038427"],
		#T cell inflamed signature gene expression profile (TIS-GEP)
		"GEP_TIS": ["ENSG00000271503", "ENSG00000139193", "ENSG00000120217", "ENSG00000103855", "ENSG00000153563", "ENSG00000174600", "ENSG00000138755", "ENSG00000172215", "ENSG00000196735", "ENSG00000196126",
				"ENSG00000204592", "ENSG00000131203", "ENSG00000089692", "ENSG00000105374", "ENSG00000197646", "ENSG00000205220", "ENSG00000115415", "ENSG00000181847"],
		"GEP_EGFR_MDM24": ["ENSG00000146648", "ENSG00000135679", "ENSG00000198625"],
		#CYT score, based on GZMA and PRF1, which are dramatically upregulated upon CD8+ T cell activation
		"CYT_SCORE": ["ENSG00000145649", "ENSG00000180644"],
		#HLA-A, HLA-B, HLA-C
		"GEP_HLA_I": ["ENSG00000206503", "ENSG00000234745", "ENSG00000204525"],
		#B2M
		"GEP_B2M": ["ENSG00000166710"],
  
		"GEP_MDM2": ["ENSG00000135679"],
		"GEP_MDM4": ["ENSG00000198625"],
		"GEP_EGFR": ["ENSG00000146648"]
	}
 
	#weights must be same order as gene IDs
	weights = {
		"GEP_TIS": [0.008346, 0.072293, 0.042853, -0.0239, 0.031021, 0.151253, 0.074135, 0.004313, 0.020091, 0.058806, 0.07175, 0.060679, 0.123895, 0.075524, 0.003734, 0.032999, 0.250229, 0.084767]
	}

 
	results = pd.DataFrame()
	for signature in signatures:
		data = expr_data.loc[signatures[signature]]

		#weighted mean
		if signature in weights:
			result = data.apply(lambda x: sum(x * weights[signature]) / sum(weights[signature]) , axis=0)
			results[signature] = result
		elif signature in ["CYT_SCORE"]:
			result = data.apply(stats.gmean, axis=0)
			results[signature] = result
		else:
			result = pd.DataFrame(data.mean(axis=0))
			results[signature] = result
 

	results.insert(0,"ID", results.index)
	if options.output == "-":
		print(results.to_csv(sep="\t", header=True, index=False))
	else:
		results.to_csv(options.output, sep="\t", header=True, index=False)


	

	

if __name__ == "__main__":
	main(sys.argv)