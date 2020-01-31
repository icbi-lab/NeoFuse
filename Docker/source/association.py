#!/usr/bin/python
import argparse, sys
import csv

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
	supported = "BoLA-6*13:01, Eqca-1*01:01, H-2-Db, H-2-Dd, H-2-Kb, H-2-Kd, H-2-Kk, H-2-Ld, HLA-A*01:01, HLA-A*02:01, HLA-A*02:02, HLA-A*02:03, HLA-A*02:05, HLA-A*02:06, HLA-A*02:07, HLA-A*02:11, HLA-A*02:12, HLA-A*02:16, HLA-A*02:17, HLA-A*02:19, HLA-A*02:50, HLA-A*03:01, HLA-A*11:01, HLA-A*23:01, HLA-A*24:02, HLA-A*24:03, HLA-A*25:01, HLA-A*26:01, HLA-A*26:02, HLA-A*26:03, HLA-A*29:02, HLA-A*30:01, HLA-A*30:02, HLA-A*31:01, HLA-A*32:01, HLA-A*33:01, HLA-A*66:01, HLA-A*68:01, HLA-A*68:02, HLA-A*68:23, HLA-A*69:01, HLA-A*80:01, HLA-B*07:01, HLA-B*07:02, HLA-B*08:01, HLA-B*08:02, HLA-B*08:03, HLA-B*14:02, HLA-B*15:01, HLA-B*15:02, HLA-B*15:03, HLA-B*15:09, HLA-B*15:17, HLA-B*18:01, HLA-B*27:02, HLA-B*27:03, HLA-B*27:04, HLA-B*27:05, HLA-B*27:06, HLA-B*35:01, HLA-B*35:03, HLA-B*37:01, HLA-B*38:01, HLA-B*39:01, HLA-B*39:06, HLA-B*40:01, HLA-B*40:02, HLA-B*42:01, HLA-B*44:02, HLA-B*44:03, HLA-B*45:01, HLA-B*46:01, HLA-B*48:01, HLA-B*51:01, HLA-B*53:01, HLA-B*54:01, HLA-B*57:01, HLA-B*58:01, HLA-B*83:01, HLA-C*03:03, HLA-C*04:01, HLA-C*05:01, HLA-C*06:02, HLA-C*07:02, HLA-C*08:02, HLA-C*12:03, HLA-C*14:02, HLA-C*15:02, Mamu-A*01:01, Mamu-A*02:01, Mamu-A*02:0102, Mamu-A*07:01, Mamu-A*07:0103, Mamu-A*11:01, Mamu-A*22:01, Mamu-A*26:01, Mamu-B*01:01, Mamu-B*03:01, Mamu-B*08:01, Mamu-B*10:01, Mamu-B*17:01, Mamu-B*17:04, Mamu-B*39:01, Mamu-B*52:01, Mamu-B*83:01, Patr-A*01:01, Patr-A*04:01, Patr-A*07:01, Patr-A*09:01, Patr-B*01:01, Patr-B*13:01, Patr-B*24:01"
	with open(optiInFile) as in_file:
		hl=[]
		uns=[]
		for line in in_file:
			h = "HLA-%s" % line
			h=h.replace("\n", " ")
			if h.replace(" ", "") in supported: # check if HLA is supported by MHCFlurry, else append it to the unsupported list
				hl.append(h)
			else:
				uns.append(h)
		hla = " ".join(hl)
		unsupported = "\n".join(uns)
	in_file.close()
	out = outDir+"_unsupported.txt"
	with open(out, "+w") as out:
		out.write("The following HLA types were predicted by OptiType, but are not supported by MHCFlurry:\n")
		out.write(unsupported)
	out.close()
	return hla

# Function to build the associations and Mhcflurry run intermidiate files
def tmp_out_pep(xenoInFile, tmpOutFile1, tmpOutFile2, outDir, fGenes, peptides, hla, ftype, confidence, stopCod, cores):
	gene_file=[]
	counter=2 # used for parallelizing mhcflurry jobs
	fileID=1 # used in printing fusion gene names and corresponding peptide sequences at the MHCFlurry run file
	if "xeno_8.tsv" in xenoInFile:
		postfix="_8.tsv"
	elif "xeno_9.tsv" in xenoInFile:
		postfix="_9.tsv"
	elif "xeno_10.tsv" in xenoInFile:
		postfix="_10.tsv"
	elif "xeno_11.tsv" in xenoInFile:
		postfix="_11.tsv"
	with open(tmpOutFile1, "+w") as out_file:
		for i in range(0, len(fGenes)):
			out=outDir+"_"+fGenes[i].split(",")[0].split("(")[0]+"_"+str(fileID)+postfix
			if not hla: # check if there are supported by MHCFlurry HLAs
				pass
			else:
				gene_file.append(fGene[i]+"#"+out+"#"+ftype[i]+"#"+confidence[i]+"#"+stopCod[i]) #get the lines for the associations file
				if counter <= cores:
					out_file.write('''mhcflurry-predict --alleles %s --peptides %s --out %s &\n''' % (hla, peptides[i], out))
				else:
					out_file.write('''mhcflurry-predict --alleles %s --peptides %s --out %s &\n''' % (hla, peptides[i], out))
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
	parser = argparse.ArgumentParser(usage='association.py [-h] -x {/path/to/xeno/input/} -l {path/to/OptiType/input} -o {/path/to/output/} -c (NCores)')
	parser.add_argument('-x','--XenoFile', help='Xeno file',required=True)
	parser.add_argument('-l','--OptiFile', help='Optitype file',required=True)
	parser.add_argument('-o','--outDir', help='Output dir', required=True)
	parser.add_argument('-c','--cores', help='Number of cores', type=int, required=False)
	args = parser.parse_args()
	fGene=[]
	peps=[]
	ftype=[]
	confidence=[]
	stopCod=[]
	hla=""
	xenoFile = args.XenoFile
	optiInFile = args.OptiFile
	outDir = args.outDir
	cores = args.cores
	outFile1 = outDir+"_TEST_OUT.sh"
	outFile2 = outDir+"_ASSOCIATIONS_OUT.txt"

	seek_fusePep(xenoFile, fGene, peps, ftype, confidence, stopCod)
	hla = seek_hla(optiInFile, outDir, hla)

	tmp_out_pep(xenoFile, outFile1, outFile2, outDir, fGene, peps, hla, ftype, confidence, stopCod, cores)
