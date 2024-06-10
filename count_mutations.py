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
	parser.add_option("--genes", action="store", dest="genes", type="string", default="PDIA3,B2M,TAP1,TAP2,TAPBP,CALR,HLA-A,HLA-B,HLA-C,HLA-E,HLA-F,HLA-G,HLA-DMA,HLA-DMB,HLA-DOA,HLA-DOB,HLA-DPA1,HLA-DPB1,HLA-DQA1,HLA-DQA2,HLA-DQB1,HLA-DRA,HLA-DRB1,HLA-DRB5,JAK1,JAK2,JAK3,TYK2,STAT1,STAT2,STAT3,STAT4,STAT5A,STAT5B,STAT6,PTEN,IFNG,IFNGR1,IFNGR2,NLRC5,CIITA,BRAF,NRAS,ERBB2,EGFR,MDM2,MDM4,SERPINB3,SERPINB4,KRAS,CTNNB1,APC,AXIN1,HNF1A,CD274,KIT,TP53,STK11,IRF1,CDK4,SPOP,CMTM4,CMTM6,PBRM1,ARID1A,SMARCA4", help="Comma-separated list of genes to be examined.")
	parser.add_option("--effects", action="store", dest="effects", type="string", default="frameshift_variant,splice_acceptor_variant,splice_donor_variant,start_lost,start_retained_variant,stop_gained,stop_lost,inframe_deletion,inframe_insertion,missense_variant,splice_region_variant", help="Comma-separated list of deleterious variant effects to be considered.")	
	parser.add_option("--use_cosmic", action="store_true", dest="use_cosmic", help="Use COSMIC annotation.")
	parser.add_option("--cadd_threshold", action="store", dest="cadd_threshold", type="int", default=20, help="Only consider variants with CADD higher than this.")
	parser.add_option("--cnv_effect", action="store", dest="cnv_effects", type="string", default="LOH,HET_DEL,HOM_DEL")
	parser.add_option("--loglikelihood_threshold", action="store", dest="loglikelihood_threshold", type="int", default=100, help="CNV loglikelihood.")
	(options, args) = parser.parse_args()

	genes_of_interest = str(options.genes).split(",") + ["BRAF_V600"]
	effects_of_interest = str(options.effects).split(",")

	data = {}
 
	vep_high_impact = ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost" , "transcript_amplification", "feature_elongation", "feature_truncation"]
 
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
   
			log_likelihood = int(parts[9])
   
			if log_likelihood < options.loglikelihood_threshold:
				continue
		
			cn = int(parts[3])
			cn_type = parts[4]

			genes = parts[15].split(",")



			for gene in genes:
				if gene in genes_of_interest and not cn_type in hits[gene]:
					content = hits[gene].copy()
					if(cn == 0):
						content.append("HOM_DEL")
					elif(cn == 1):
						content.append("HET_DEL")
					elif(cn_type == "LOH"):
						content.append(cn_type)
					elif(cn >= 4):
						content.append("AMP")
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

			cmc_anno = pd.NA
			if "COSMIC_CMC" in snv.info.keys():
				parts = snv.info["COSMIC_CMC"].split("|")
				if parts[-1] != "Other":
					cmc_anno = int(parts[-1])

			if any(item in snv.filter.keys() for item in ["depth-nor", "depth-tum", "freq-nor", "freq-tum", "LowDepth", "LowEVS", "lt-3-reads", "special-chromosome"]):
				continue

			#BRAF V600E
			if snv.chrom == "chr7" and snv.pos == 140753336:
				data[id]["BRAF_V600"] = ["V600_positive"]



			for transcript in vep_transcripts:
				parts = transcript.split("|")
				gene = parts[3]
				if not gene in genes_of_interest:
					continue
				effects = parts[1].split("&")

				if not options.use_cosmic:
					for effect in effects:
						if not effect in effects_of_interest:
							continue
						if effect in data[id][gene]:
							continue

						if not pd.isna(cadd_anno) and int(cadd_anno) < options.cadd_threshold:
							continue
						#Not all indels have a cadd anno. We only use VEP HIGH impact features in tis case
						if pd.isna(cadd_anno):
							if not effect in vep_high_impact:
								continue

						content = data[id][gene].copy()
						content.append(effect)
						data[id][gene] = content
						break
				else:
					for effect in effects:
						if not effect in vep_high_impact:
							if pd.isna(cmc_anno):
								continue
						effect_tmp = effect
						if not pd.isna(cmc_anno):
							effect_tmp = str(cmc_anno)
						if effect_tmp in data[id][gene]:
							continue
						content = data[id][gene].copy()
						content.append(effect_tmp)
						data[id][gene] = content
						break

						

	snv_paths.close()

	df = pd.DataFrame.from_dict(data, orient="index")
	df = df.applymap(lambda x: ",".join(x))
	df["TUMOR"] = df.index
	df[["TUMOR","NORMAL"]] = df["TUMOR"].str.split("-", n=1, expand=True)
	df.reset_index(drop=True, inplace=True)

	# Shifting column "ABC" to the left
	df = pd.concat([df['NORMAL'], df.drop(columns=['NORMAL'])], axis=1)
	df = pd.concat([df['TUMOR'], df.drop(columns=['TUMOR'])], axis=1)
	
	#analyze pathways
	ctnnb1_pathway = []
	ifn_jak_stat_pathway = []
	b2m_pathway = []
	pten_pathway = []
	stk11_pathway = []
	tp53_pathway = []
	mdm2_pathway = []
	mdm4_pathway = []
	egfr_pathway = []
	kras_pathway = []
	
	#response pathways
	cd274_pathway = []
	pathway_serpinb3_4 = []
	pathway_chromatin_modifier = []
 
	for i,row in df.iterrows():
		#Beta catenin: CTNNB1 activated or negative regulators inactvated
		beta_catenin_mut = False
		if any(alteration in row["CTNNB1"] for alteration in ["1", "2", "3"]):
			beta_catenin_mut = True
		if any(alteration in row["APC"] for alteration in vep_high_impact + ["1", "2", "3", "HOM_DEL"]):
			beta_catenin_mut = True
		if any(alteration in row["AXIN1"] for alteration in vep_high_impact + ["1", "2", "3", "HOM_DEL"]):
			beta_catenin_mut = True
		if any(alteration in row["HNF1A"] for alteration in vep_high_impact + ["1", "2", "3", "HOM_DEL"]):
			beta_catenin_mut = True

		ctnnb1_pathway.append(beta_catenin_mut)
  
		#Interferon II to JAK-STAT
		jak_stat_mut = False
		if any(alteration in row["JAK1"] for alteration in vep_high_impact + ["1", "2", "3", "HOM_DEL"]):
			jak_stat_mut = True
		if any(alteration in row["JAK2"] for alteration in vep_high_impact + ["1", "2", "3", "HOM_DEL"]):
			jak_stat_mut = True
		ifn_jak_stat_pathway.append(jak_stat_mut)
  
		b2m_mut = False
		if any(alteration in row["B2M"] for alteration in vep_high_impact + ["1", "2", "3", "HET_DEL", "LOH"]): #het deletions only, hom dels. lead to fully unfunctional MHC complex: NK would exterminate
			b2m_mut = True
		b2m_pathway.append(b2m_mut)

		pten_mut = False
		if any(alteration in row["PTEN"] for alteration in vep_high_impact + ["1", "2", "3","HOM_DEL","HET_DEL"]):
			pten_mut = True
		pten_pathway.append(pten_mut)
  
		stk11_mut = False	
		if any(alteration in row["STK11"] for alteration in vep_high_impact + ["1", "2", "3","HOM_DEL","HET_DEL"]):
			stk11_mut = True
		stk11_pathway.append(stk11_mut)
  
		tp53_mut = False
		if any(alteration in row["TP53"] for alteration in vep_high_impact + ["1", "2", "3","HOM_DEL", "HET_DEL", "LOH"]):
			tp53_mut = True
		tp53_pathway.append(tp53_mut)
   
		mdm2_mut = False
		if any(alteration in row["MDM2"] for alteration in ["1", "2", "3","AMP"]):
			mdm2_mut = True
		mdm2_pathway.append(mdm2_mut)

		mdm4_mut = False
		if any(alteration in row["MDM4"] for alteration in ["1", "2", "3","AMP"]):
			mdm4_mut = True
		mdm4_pathway.append(mdm4_mut)

		egfr_mut = False
		if any(alteration in row["EGFR"] for alteration in ["1", "2", "3","AMP"]):
			egfr_mut = True
		egfr_pathway.append(egfr_mut)
  
		kras_mut = False
		if any(alteration in row["KRAS"] for alteration in ["1", "2", "3","AMP"]):
			kras_mut = True
		kras_pathway.append(kras_mut)
  
		cd274_mut = False
		if any(alteration in row["CD274"] for alteration in ["AMP"]):
			cd274_mut = True
		if any(alteration in row["SPOP"] for alteration in vep_high_impact +["1", "2", "3","DEL"]):
			cd274_mut = True
		if any(alteration in row["CMTM4"] for alteration in vep_high_impact + ["1", "2", "3","DEL"]):
			cd274_mut = True
		if any(alteration in row["CMTM6"] for alteration in vep_high_impact +["1", "2", "3","DEL"]):
			cd274_mut = True
	
		cd274_pathway.append(cd274_mut)
  
		serpinb3_4_mut = False
		if any(alteration in row["SERPINB3"] for alteration in vep_high_impact + ["1", "2", "3","DEL"]):
			serpinb3_4_mut = True
		if any(alteration in row["SERPINB4"] for alteration in vep_high_impact + ["1", "2", "3","DEL"]):
			serpinb3_4_mut = True
		pathway_serpinb3_4.append(serpinb3_4_mut)
  
		chrom_mut = False
		if any(alteration in row["PBRM1"] for alteration in vep_high_impact + ["1", "2", "3","DEL"]):
			chrom_mut = True
		if any(alteration in row["ARID1A"] for alteration in vep_high_impact + ["1", "2", "3","DEL"]):
			chrom_mut = True
		if any(alteration in row["SMARCA4"] for alteration in vep_high_impact + ["1", "2", "3","DEL"]):
			chrom_mut = True
		pathway_chromatin_modifier.append(chrom_mut)


	df["PATHWAY_MUT_BETA_CATENIN"] = ctnnb1_pathway
	df["PATHWAY_MUT_JAK"] = ifn_jak_stat_pathway
	df["PATHWAY_MUT_B2M"] = b2m_pathway
	df["PATHWAY_MUT_PTEN"] = pten_pathway
	df["PATHWAY_MUT_STK11"] = stk11_pathway
	df["PATHWAY_MUT_TP53"] = tp53_pathway
	df["PATHWAY_MUT_MDM2"] = mdm2_pathway
	df["PATHWAY_MUT_MDM4"] = mdm4_pathway
	df["PATHWAY_MUT_EGFR"] = egfr_pathway
	df["PATHWAY_MUT_KRAS"] = kras_pathway
	df["PATHWAY_MUT_CD274"] = cd274_pathway
	df["PATHWAY_MUT_SERPINB3_4"] = pathway_serpinb3_4
	df["PATHWAY_MUT_CHROMATIN_MODIFIER"] = pathway_chromatin_modifier
	df.to_csv(options.out, sep="\t", index=False)


	
if __name__ == "__main__":
	main(sys.argv)