#!/usr/bin/python
import argparse, sys
import csv

def desperation(inFile1, inFile2, outFile):
	gene=[]
	tpm=[]
	gene1tpm=[]
	gene2tpm=[]
	hlaTpm=[]
	tpmAvg=[]
	fusion=[]
	fpep=[]
	mIc=[]
	mRank=[]
	eveType=[]
	stopCod=[]
	arConf=[]

	with open(inFile1) as in_file:
		csv_reader = csv.reader(in_file, delimiter='\t')
		in_file.readline()
		for row in csv_reader:
			gene.append(row[0])
			if row[2] == "0":
				tpm.append("NA")
			else:
				tpm.append(row[2])
	in_file.close()

	with open(inFile2) as in_file:
		csv_reader = csv.reader(in_file, delimiter='\t')
		in_file.readline()
		for row in csv_reader:
			fusion.append(row[0])
			fpep.append(row[4])
			mIc.append(row[5])
			mRank.append(row[6])
			eveType.append(row[7])
			stopCod.append(row[8])
			arConf.append(row[9])
			if "," in row[1]:
				gene1tpm.append(row[1]+"\t"+"NA")
			elif "," in row[2]:
					gene2tpm.append(row[2]+"\t"+"NA")
			for i in range(0, len(tpm)):
				if gene[i] == row[1]:
					if tpm[i] == "NA":
						gene1tpm.append(row[1]+"\t"+tpm[i])
					elif tpm[i] != "NA":
						gene1tpm.append(row[1]+"\t"+"%.2f" % float(tpm[i]))
				if gene[i] == row[2]:
					if tpm[i] == "NA":
						gene2tpm.append(row[2]+"\t"+tpm[i])
					elif tpm[i] != "NA":
						gene2tpm.append(row[2]+"\t"+"%.2f" % float(tpm[i]))
				if gene[i] == row[3].split("*")[0]:
					if tpm[i] == "NA":
						hlaTpm.append(row[3]+"\t"+tpm[i])
					elif tpm[i] != "NA":
						hlaTpm.append(row[3]+"\t"+"%.2f" % float(tpm[i]))
	in_file.close()

	for i in  range(0, len(gene1tpm)):
		if "NA" in gene1tpm[i].split("\t")[1] or "NA" in gene2tpm[i].split("\t")[1]:
			tpmAvg.append("NA")
		else:
			tpmAvg.append("%.2f" % float((float(gene1tpm[i].split("\t")[1])+float(gene2tpm[i].split("\t")[1]))/2.0))

	with open(outFile, "+w") as out_file:
		out_file.write("Fusion\tGene1\tGene2\tHLA_Type\tFusion_Peptide\tIC50\tRank\tEvent_Type\tStop_Codon\tConfidence\tGene1_TPM\tGene2_TPM\tAvg_TPM\tHLA_TPM\n")
		for j in range(0, len(fusion)):
			out_file.write(fusion[j]+"\t"+gene1tpm[j].split("\t")[0]+"\t"+gene2tpm[j].split("\t")[0]+"\t"+hlaTpm[j].split("\t")[0]+"\t"+fpep[j]+"\t"+mIc[j]+"\t"+mRank[j]+"\t"+eveType[j]+"\t"+stopCod[j]+"\t"+arConf[j]+"\t"+gene1tpm[j].split("\t")[1]+"\t"+gene2tpm[j].split("\t")[1]+"\t"+tpmAvg[j]+"\t"+hlaTpm[j].split("\t")[1]+"\n")
		out_file.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(usage='final_out.py [-h] -t {/path/to/tmp/TPMinput/} -t {/path/to/tmp/filtered_input/} -o {/path/to/output/}')
	parser.add_argument('-t','--tpmInput', help='input file',required=True)
	parser.add_argument('-i','--Input', help='input file',required=True)
	parser.add_argument('-o','--Output', help='Output file',required=True)
	args = parser.parse_args()
	inFile1=args.tpmInput
	inFile2=args.Input
	outFile=args.Output

	desperation(inFile1, inFile2, outFile)
