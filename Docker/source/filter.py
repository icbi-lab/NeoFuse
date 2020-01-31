#!/usr/bin/python
import argparse, sys
import csv

def purge(inFile, outFile, threshold, rank, confidence):
    inRows=[]
    outRows=[]
    with open(inFile) as in_file:
        csv_reader = csv.reader(in_file, delimiter='\t')
        in_file.readline()
        for row in csv_reader:
            inRows.append(row)
    in_file.close()
    
    for row in inRows:
        if float(row[5]) <= float(threshold) and float(row[5]) <= float(rank):
            for conf in confidence:
                if conf in row[9]:
                    outRows.append(row)
        else:
            pass

    with open(outFile, "+w") as out_file:
        out_file.write("Fusion\tGene1\tGene2\tHLA_Type\tFusion_Peptide\tIC50\tRank\tEvent_Type\tStop_Codon\tConfidence\n")
        for row in outRows:
            out_file.write("%s\n" % "\t".join(row))
    out_file.close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(usage='filter.py [-h] -i {/path/to/tmp/input/} -o {/path/to/output/} -t {affinity threshold} -T {rank threshold} -c {confidence level})')
	parser.add_argument('-i','--Input', help='input file',required=True)
	parser.add_argument('-o','--Output', help='Output file',required=True)
	parser.add_argument('-t','--AffinityThreshold', help='netMHCpan prediction threshold', required=True)
	parser.add_argument('-T','--RankThreshold', help='netMHCpan rank threshold', required=True)
	parser.add_argument('-c','--ConfidenceLevel', help='Arriba confidence level threshold', nargs='+', type=str,  required=False)
	args = parser.parse_args()
	inFile = args.Input
	outFile = args.Output
	thr = args.AffinityThreshold
	rank = args.RankThreshold
	conf = args.ConfidenceLevel

	purge(inFile, outFile, thr, rank, conf)
