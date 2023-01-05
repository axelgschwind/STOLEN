# Author: Axel Gschwind
# program: remove_duplicate_read_names.py
# purpose: Remove duplicate read names from fastq.gz
import gzip
import os
import io
import sys
from optparse import OptionParser

def main(argv):
	parser = OptionParser(usage="usage: %prog -in FILE.fastq.gz -out", description="Removes duplicate read names from a fastq file. Only the first read is kept and subsequent reads with the same name are discarded.")
	parser.add_option("-i", "--in", dest="input", help="input Fastq")
	parser.add_option("-o", "--out", dest="output", help="output Fastq")

	(options, args) = parser.parse_args()

	if not os.path.exists(options.input):
		raise Exception("Input file %s does not exist!" % options.input)

	with io.TextIOWrapper(io.BufferedReader( gzip.open(options.input, "r"), buffer_size=262144)) as f:
		with gzip.open(options.output, "w") as out:
			ids = set()
			duplicate_read_name = False

			count_total = 0
			count_dismissed = 0

			duplicate_read_name = False
			line = "bla"
			while line != bytes("", "ascii"):
				line = f.readline().encode("ascii").strip()
				if line == bytes("", "ascii"):
					break
				count_total += 1
				if line in ids:
					duplicate_read_name = True
					count_dismissed += 1
				else:
					duplicate_read_name = False
					ids.add(line)

				#out block will contain 4 lines from input file
				out_block = [line, f.readline().encode("ascii").strip(), f.readline().encode("ascii").strip(), f.readline().encode("ascii").strip()]
				if not duplicate_read_name:
					out.writelines(line + bytes("\n", "ascii") for line in out_block)
			print("PARSED %i READS" % count_total)
			print("DISCARDED %i READS" % count_dismissed)

if __name__ == "__main__":
	main(sys.argv)
