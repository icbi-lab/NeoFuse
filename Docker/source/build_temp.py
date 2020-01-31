#!/usr/bin/python
import argparse, sys
import csv

def final_out(infileAssoc, outfile):
	fGene=[]
	files=[]
	gene1=[]
	gene2=[]
	ftype=[]
	confidence=[]
	stopCod=[]
	alleles=[]
	peptide=[]
	prediction=[]
	percentile=[]
	gene1=[]
	gene2=[]
	i=0

	with open(infileAssoc) as in_file:
		for line in in_file:
			fGene.append(line.split("#")[0])
			gene1.append(line.split("#")[0].split("_")[0])
			gene2.append(line.split("#")[0].split("_")[1])
			files.append(line.split("#")[1])
			ftype.append(line.split("#")[2])
			confidence.append(line.split("#")[3])
			stopCod.append(line.split("#")[4].replace("\n", ""))
	in_file.close()

	with open(outfile, "+w") as out_file:
		out_file.write("Fusion\tGene1\tGene2\tHLA_Type\tFusion_Peptide\tIC50\tRank\tEvent_Type\tStop_Codon\tConfidence\n")
		for file in files:
			with open(file) as csv_file:
				csv_reader = csv.reader(csv_file, delimiter=',')
				csv_file.readline() # skip the header line
				for row in csv_reader:
					alleles.append(row[0])
					peptide.append(row[1])
					prediction.append(row[2])
					percentile.append(row[5])
					out_file.write(fGene[i].replace("_", "-")+"\t"+gene1[i]+"\t"+gene2[i]+"\t"+row[0]+"\t"+row[1]+"\t"+row[2]+"\t"+row[5]+"\t"+ftype[i]+"\t"+stopCod[i]+"\t"+confidence[i]+"\n")
			i+=1
			csv_file.close()
	out_file.close()

if __name__ == "__main__":
	# Parse arguments
	parser = argparse.ArgumentParser(usage='build_temp.py [-h] -a {path/to/Associations.txt/input} -o {/path/to/output/}')
	parser.add_argument('-a','--AssociationsFile', help='Assosciations tmp file',required=True)
	parser.add_argument('-o','--outDir', help='Output dir', required=True)
	args = parser.parse_args()

	inFile = args.AssociationsFile
	outFile = args.outDir+"_final.tsv"

	final_out(inFile, outFile)
