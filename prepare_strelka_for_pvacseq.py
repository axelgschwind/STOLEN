# Author: Axel Gschwind
# program: prepare_strelka_for_pvacseq.py
# purpose: Add annotations in strelka2 VCF files for use with pVACseq

import pysam
from optparse import OptionParser
import sys


def main(argv):
	if len(argv) == 1:
		argv.append("--help")
	usage = "usage: %prog -i input.vcf.gz -o ready_for_pVACseq.vcf.gz"
	desc = "Prepares somatic strelka2 VCF files for use with pVACseq."
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--input", action="store", dest="input", type="string", help="Input somatic Strelka2 VCF")
	parser.add_option("-o", "--output", action="store", dest="output", type="string", default="-", help="Output VCF")

	(options, args) = parser.parse_args()

	#Add Format fields for AF and AD for all samples
	vcf = pysam.VariantFile( options.input , "r")

	#Add genotype information dummy
	vcf.header.formats.add("GT", 1, "Float", "Genotype")
	vcf.header.formats.add("AF", 1, "Float", "Variant Allele Frequency.")
	vcf.header.formats.add("AD", 1, "Float", "Alternate depth of the SNV.")

	new_vcf = pysam.VariantFile(options.output, "w", header=vcf.header)
 
 
	for snv in vcf:
		for sample in snv.samples:
			#Add genotype dummy
			snv.samples[sample]["GT"] = (0,1)
			
			#Add AF and AD
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

		new_vcf.write(snv)
    
	vcf.close()
	new_vcf.close()

if __name__ == "__main__":
	main(sys.argv)