# Author: Axel Gschwind
# program: prepare_strelka_for_pvacseq.py
# purpose: Add annotations in strelka2 VCF files for use with pVACseq

import pysam
import pysamstats
from optparse import OptionParser
import sys
import pandas as pd
variant_types = ["splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion",
		"missense_variant", "protein_altering_variant", "splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
		"incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant", "NMD_transcript_variant"]

def is_exonic(vep_info):
	global variant_types
	is_exonic = False

	for vep in vep_info:
		parts = vep.split("|")
		if parts[1] in variant_types:
			is_exonic = True
	return is_exonic



def main(argv):
	if len(argv) == 1:
		argv.append("--help")
	usage = "usage: %prog -i input.vcf.gz -o ready_for_pVACseq.vcf.gz"
	desc = "Prepares somatic strelka2 VCF files for use with pVACseq."
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--input", action="store", dest="input", type="string", help="Input somatic Strelka2 VCF")
	parser.add_option("-o", "--output", action="store", dest="output", type="string", default="-", help="Output VCF")
	parser.add_option("--rna_bam", action="store", dest="rna_bam", type="string", default="", help="RNA BAM file")
	parser.add_option("--expr_file", action="store", dest="expr_file", type="string", default="", help="RNA expression file")
	(options, args) = parser.parse_args()

	#Add Format fields for AF and AD for all samples
	vcf = pysam.VariantFile( options.input , "r")

	#Add genotype information dummy
	vcf.header.formats.add("GT", 1, "Float", "Genotype")
	vcf.header.formats.add("AF", 1, "Float", "Variant Allele Frequency.")
	vcf.header.formats.add("AD", 1, "Integer", "Alternate depth of the SNV.")

	
	if options.rna_bam != "":
		rnafile = pysam.AlignmentFile(options.rna_bam, "rb")
		vcf.header.formats.add("RDP", 1, "Integer", "RNA total read depth")
		vcf.header.formats.add("RAD", 1, "Integer", "RNA alternate allele read depth")
		vcf.header.formats.add("RAF", 1, "Float", "RNA variant allele frequency")
  
	if options.expr_file != "":
		expressionfile = pd.read_csv(options.expr_file, delimiter="\t", index_col=0)
		vcf.header.formats.add("GX", ".", "String", "Gene Expression")

	new_vcf = pysam.VariantFile(options.output, "w", header=vcf.header)
 
	for snv in vcf:
		for sample in snv.samples:
			#Add genotype dummy
			snv.samples[sample]["GT"] = (0,1)
			
			#Add tumor/normal AF and AD
			snv_fields = set(["AU", "CU", "GU", "TU"])
			indel_fields = set(["TIR", "TAR"])
			if snv_fields.issubset( set(snv.samples[sample].keys()) ) :
				(A1,A2) = snv.samples[sample]["AU"]
				(C1,C2) = snv.samples[sample]["CU"]
				(G1,G2) = snv.samples[sample]["GU"]
				(T1, T2) = snv.samples[sample]["TU"]
				counts = {"A" : A1, "C": C1, "G": G1, "T": T1}
				depth =  sum(counts.values())
    
				obs = snv.alts[0]
    
				#variant allele frequency
				af = counts[obs] / depth
				#alternate depth
				ad = counts[obs]
				snv.samples[sample]["AF"] = af
				snv.samples[sample]["AD"] = ad
			elif indel_fields.issubset( set(snv.samples[sample].keys()) ):
				#TIR: Reads strongly supporting indel allele for tiers 1,2
				TIR = snv.samples[sample]["TIR"][0]
				#TAR: Reads strongly supporting alternate allele for tiers 1,2
				TAR = snv.samples[sample]["TAR"][0]

				af = 0.
				if TIR+TAR != 0:
					af = TIR/(TIR+TAR)
				snv.samples[sample]["AF"] = af
				snv.samples[sample]["AD"] = TIR
			else:
				print("Could not parse STRELKA variant:", snv)
				exit(1)

			#Annoate read depth from RNA BAM
			if options.rna_bam != "":
				if not is_exonic(snv.info["CSQ"]):
						snv.samples[sample]["RDP"] = 0
						snv.samples[sample]["RAD"] = 0
						snv.samples[sample]["RAF"] = 0
				else:
					a = pysamstats.stat_variation(rnafile, fafile="/mnt/storage2/megSAP/data/genomes/GRCh38.fa", chrom=snv.chrom, start=snv.start, end=snv.stop, truncate=True, pad=True)
					for rec in a:
						if(len(snv.ref[0]) == 1 and len(snv.alts[0]) == 1):
							rna_dp = rec["matches"] + rec["mismatches"]
							rna_alt_dp = rec[snv.alts[0]]
							rna_af = 0.
							if rna_dp > 0: 
								rna_af = rna_alt_dp/rna_dp
							snv.samples[sample]["RDP"] = rna_dp
							snv.samples[sample]["RAD"] = rna_alt_dp
							snv.samples[sample]["RAF"] = rna_af
						else:
							rna_dp = rec["matches"] + rec["mismatches"]
							indel_count = rec["insertions"] +  rec["deletions"]
							rna_af = 0
							if rna_dp > 0:
								rna_af = indel_count / rna_dp

							snv.samples[sample]["RDP"] = rna_dp
							snv.samples[sample]["RAD"] = indel_count
							snv.samples[sample]["RAF"] = rna_af

			#Annotate gene expression data
			if options.expr_file != "":
				vep_info = snv.info["CSQ"]

				#annotate expression of first gene in list (pvactools can only parse one value per SNV)
				gene_expression = {}
				for vep in vep_info:
					parts = vep.split("|")
					gene_id = parts[4]
					if not gene_id in expressionfile.index:
						continue
					if gene_id in gene_expression.keys():
						continue
					gene_expression[gene_id] = gene_id + "|" +str(expressionfile.loc[gene_id]["tpm"].item())


				if len(gene_expression) > 0:
					snv.samples[sample]["GX"] = ",".join(gene_expression.values())
					#snv.samples[sample]["GX"] = list( gene_expression.values() )[0] #Take only one, arbitrary first value. See https://github.com/griffithlab/pVACtools/issues/991

		new_vcf.write(snv)
	vcf.close()
	new_vcf.close()

if __name__ == "__main__":
	main(sys.argv)