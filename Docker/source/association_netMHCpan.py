#!/usr/bin/python
import argparse, sys
import csv, re

# Function to retrieve fusion gene names, peptide sequences and confidence level
def seek_fusePep(xenoInFile, fusionGene, peps, ftype, confidence, stopCod):
	with open(xenoInFile) as in_file:
		for line in in_file:
			gene1 = line.split("#", 1)[0].split("\t", 1)[0].split(" - ")[0]
			gene2 = line.split("#", 1)[0].split("\t", 1)[0].split(" - ")[1]
			fgene = f"{gene1}_{gene2}"
			fusionGene.append(fgene)
			peps.append(line.split("#", 1)[0].split("\t", 1)[1].replace("\t", " ").upper())
			ftype.append(line.split("#", 2)[1])
			confidence.append(line.split("#")[2])
			stopCod.append(line.split("#")[3].replace("\n", ""))
	in_file.close()
	return fusionGene, peps, ftype, confidence, stopCod

# Function to retrive HLAs predicted by OptiType
def seek_hla(optiInFile, outDir, hla):
	with open(optiInFile) as in_file:
		hl=[]
		for line in in_file:
			h = "HLA-%s" % line.replace("*", "")
			h=h.replace("\n", ",")
			hl.append(h)
		hla = "".join(hl)
	in_file.close()
	return hla

# Function to build the associations and Mhcflurry run intermidiate files
def tmp_out_pep(xenoInFile, tmpOutFile1, tmpOutFile2, outDir, fGenes, peptidesFile, hla, ftype, confidence, stopCod, cores, pepList, path):
	gene_file=[]
	counter=2 # used for parallelizing mhcflurry jobs
	fileID=1 # used in printing fusion gene names and corresponding peptide sequences at the MHCFlurry run file
	if "xeno_8.tsv" in xenoInFile:
		postfix1="_8.tsv"
		postfix2="_8.txt"
	elif "xeno_9.tsv" in xenoInFile:
		postfix1="_9.tsv"
		postfix2="_9.txt"
	elif "xeno_10.tsv" in xenoInFile:
		postfix1="_10.tsv"
		postfix2="_10.txt"
	elif "xeno_11.tsv" in xenoInFile:
		postfix1="_11.tsv"
		postfix2="_11.txt"
	with open(tmpOutFile1, "+w") as out_file:
		for i in range(0, len(fGenes)):
			out=outDir+"_"+fGenes[i].split(",")[0].split("(")[0]+"_"+str(fileID)+postfix1
			pepFile=outDir+"_"+fGenes[i].split(",")[0].split("(")[0]+"_"+str(fileID)+postfix2
			pepList.append(pepFile)
			gene_file.append(fGene[i]+"#"+out+"#"+ftype[i]+"#"+confidence[i]+"#"+stopCod[i]) #get the lines for the associations file
			if counter <= cores:
				out_file.write('''%s/Linux_x86_64/bin/netMHCpan -a %s -p %s -BA > %s &\n''' % (path, hla, pepFile, out))
			else:
				out_file.write('''%s/Linux_x86_64/bin/netMHCpan -a %s -p %s -BA > %s &\n''' % (path, hla, pepFile, out))
				out_file.write("wait\n")
				counter = 1
			counter += 1
			fileID += 1
	out_file.close()
	# write the associations file - to be used later for multiplexing the output files into a single final output
	with open(tmpOutFile2, "+w") as out_file2:
		for j in gene_file:
			out_file2.write("%s\n" % j)
	out_file2.close()

if __name__ == "__main__":
	# Parse arguments
	parser = argparse.ArgumentParser(usage='association_netMHCpan.py [-h] -x {/path/to/xeno/input/} -l {path/to/OptiType/input} -o {/path/to/output/} -p {path/to/netMHCpan}-c (NCores)')
	parser.add_argument('-x','--XenoFile', help='Xeno file',required=True)
	parser.add_argument('-l','--OptiFile', help='Optitype file',required=True)
	parser.add_argument('-o','--outDir', help='Output dir', required=True)
	parser.add_argument('-p','--path', help='natMHCpan dir', required=True)
	parser.add_argument('-c','--cores', help='Number of cores', type=int, required=False)
	args = parser.parse_args()
	fGene=[]
	peps=[]
	ftype=[]
	confidence=[]
	stopCod=[]
	pepList=[]
	hla=""
	xenoFile = args.XenoFile
	optiInFile = args.OptiFile
	outDir = args.outDir
	path = args.path
	cores = args.cores
	outFile1 = outDir+"_TEST_OUT.sh"
	outFile2 = outDir+"_ASSOCIATIONS_OUT.txt"

	seek_fusePep(xenoFile, fGene, peps, ftype, confidence, stopCod)
	hla = seek_hla(optiInFile, outDir, hla)

	tmp_out_pep(xenoFile, outFile1, outFile2, outDir, fGene, peps, hla, ftype, confidence, stopCod, cores, pepList, path)

	for i in range(0, len(pepList)):
		with open(pepList[i], "+w") as out_file:
			out_file.write(peps[i].replace(" ", "\n"))