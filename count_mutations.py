# Author: Axel Gschwind
# program: count_mutations.py
# purpose: Merges gene counts files from megSAP RNA pipeline into one tsv file

import sys
from optparse import OptionParser
import pysam
import os.path
import pandas as pd

def main(argv):
	if len(argv) == 1:
			argv.append("--help")
	usage = "usage: %prog -i INPUT_PATHS.txt -o OUTPUT.tsv"
	desc = "Counts and summarizes alterations in pre-defined genes"
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--in", action="store", dest="input", type="string", help="Input file containing paths to VEP-annotated vcf files.")
	parser.add_option("--in_cnvs", action="store", dest="in_cnvs", type="string", help="Input file containing paths to VEP-annotated vcf files.")
	parser.add_option("-o", "--out", action="store", dest="out", type="string", help="Output TSV file containing mutation counts per genes.")
	parser.add_option("--genes", action="store", dest="genes", type="string", default="B2M,TAP1,TAP2,TAPBP,PDIA3,CALR,HLA-A,HLA-B,HLA-C,HLA_E,HLA-F,HLA-G", help="Comma-separated list of genes to be examined.")
	parser.add_option("--effects", action="store", dest="effects", type="string", default="frameshift_variant,splice_acceptor_variant,splice_donor_variant,start_lost,start_retained_variant,stop_gained,stop_lost,inframe_deletion,inframe_insertion,missense_variant,splice_region_variant", help="Comma-separated list of deleterious variant effects to be considered.")
	parser.add_option("--cadd_threshold", action="store", dest="cadd_threshold", type="int", default=20, help="Only consider variants with CADD higher than this.")
	parser.add_option("--cnv_effect", action="store", dest="cnv_effects", type="string", default="LOH,HET_DEL,HOM_DEL")
	(options, args) = parser.parse_args()

	genes_of_interest = str(options.genes).split(",")
	effects_of_interest = str(options.effects).split(",")

	data = {}
 
	cnv_paths = open(options.in_cnvs, "r")
	for path in cnv_paths:
		line = path.strip()
		if line == "":
			continue
		id = os.path.basename(line).replace("CNAs_","").replace(".txt","")
		data[id] = dict.fromkeys(genes_of_interest,[])
	cnv_paths.seek(0)

	snv_paths = open(options.input, "r") 
	for path in snv_paths:
		line = path.strip()
		if line == "":
			continue
		id = os.path.basename(line).replace(".vcf.gz", "")	
		if not id in data.keys():
			data[id] = dict.fromkeys(genes_of_interest,[])
	snv_paths.seek(0)


	
	for path in cnv_paths:
		line = path.strip()
		if line == "":
			continue
		id = os.path.basename(line).replace("CNAs_","").replace(".txt","")
		hits = dict.fromkeys(genes_of_interest,[])
		
		f = open(path.strip(), "r")
		for line in f:
			if line.startswith("#"):
				continue
			parts = line.strip().split("\t")
		
			cn_type = parts[4]

			if cn_type == "LOH" or cn_type == "DEL":
				genes = parts[15].split(",")

				for gene in genes:
					if gene in genes_of_interest and not cn_type in hits[gene]:
						content = hits[gene].copy()
						content.append(cn_type)
						hits[gene] = content
		data[id] = hits

	for path in snv_paths:
		line = path.strip()
		if line == "":
			continue
		id = os.path.basename(line).replace(".vcf.gz", "")
  
		vcf = pysam.VariantFile( line, "r" )
		
		hits = dict.fromkeys(genes_of_interest, 0)
		hits["ID"] = id
		for snv in vcf:
			vep_transcripts = snv.info["CSQ"]
			cadd_anno = pd.NA
			if "CADD_SNV" in snv.info.keys():
				cadd_anno = snv.info["CADD_SNV"]
			elif "CADD_INDEL" in snv.info.keys():
				cadd_anno = snv.info["CADD_INDEL"]

			if any(item in snv.filter.keys() for item in ["depth-nor", "depth-tum", "freq-nor", "freq-tum", "LowDepth", "LowEVS", "lt-3-reads", "special-chromosome"]):
				continue

			for transcript in vep_transcripts:
				parts = transcript.split("|")
				gene = parts[3]
				if not gene in genes_of_interest:
					continue
				effects = parts[1].split("&")

				for effect in effects:
					if not effect in effects_of_interest:
						continue
					if effect in data[id][gene]:
						continue
  
					if not pd.isna(cadd_anno) and cadd_anno < options.cadd_threshold:
						continue
					
					#Not all indels have a cadd anno. We only use VEP HIGH impact features in tis case
					if pd.isna(cadd_anno):
						if not effect in ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost" , "transcript_amplification", "feature_elongation", "feature_truncation"]:
							continue

					content = data[id][gene].copy()
					content.append(effect)
					data[id][gene] = content
					break

	snv_paths.close()

	df = pd.DataFrame.from_dict(data, orient="index")
	df = df.applymap(lambda x: ",".join(x))
	df.to_csv(options.out, sep="\t")


	
if __name__ == "__main__":
	main(sys.argv)