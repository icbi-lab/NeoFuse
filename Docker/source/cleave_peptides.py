#!/usr/bin/python
import argparse, sys
import csv
import re


# Define the functions
## Rolling window
def rolling_window(pepSeq, window):
    return [pepSeq[i : i + window] for i in range(len(pepSeq) - (window - 1))]


## Rolling window for in-frame (+1 counter)
def rolling_window_infr(pepSeq, window):
    return [pepSeq[i : i + 1 + window] for i in range(len(pepSeq) - window)]


## Remove duplicates from list function
def remove_dups(myList):
    final_list = []
    for unique in myList:
        if unique not in final_list:
            final_list.append(unique)
    return final_list


## Function to get fusion gene, breakpoints, protein sequence, event type and confidence level from input file
def sequence_hound(
    inFile,
    fusionGene,
    bkpoint1,
    bkpoint2,
    split_read1,
    split_read2,
    disco_reads,
    closest_bkp1,
    closest_bkp2,
    pepSeq,
    ftype,
    confidence,
    stopCod,
):
    drows = []
    with open(inFile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter="\t")  # parser
        csv_file.readline()  # skip the header line
        for row in csv_reader:
            drows.append(row)

        rows = remove_dups(drows)

        for row in rows:
            # Check for early stop codons, sequences between the breakpoints and not translated transcripts
            if (
                re.match(".+[|][*]", row[28])
                or re.match(".*[*].*[|].*", row[28])
                or re.match(".*[|].*[|].*", row[28])
                or "." in row[28]
            ):
                pass
            elif str.startswith(row[28], "*") or "|" not in row[28]:
                pass
            else:
                if re.match(".+[|].+[*]", row[28]):
                    stopCod.append("yes")
                else:
                    stopCod.append("no")
                fusionGene.append(f"{row[0]} - {row[1]}")
                bkpoint1.append(row[4])
                bkpoint2.append(row[5])
                split_read1.append(row[9])
                split_read2.append(row[10])
                disco_reads.append(row[11])
                confidence.append(row[14])
                closest_bkp1.append(row[18])
                closest_bkp2.append(row[19])
                if row[15] == "in-frame":
                    ftype.append("Fusion")
                elif row[15] == "out-of-frame":
                    ftype.append("Fusion-out-of-frame")
                pepSeq.append(row[28])
    csv_file.close()
    return (
        fusionGene,
        bkpoint1,
        bkpoint2,
        split_read1,
        split_read2,
        disco_reads,
        closest_bkp1,
        closest_bkp2,
        pepSeq,
        ftype,
        confidence,
    )


## Sanitize sequences - cleave at specified length
def window_onslaught(
    outFile,
    pepSeq,
    fGene,
    bkpoint1,
    bkpoint2,
    split_read1,
    split_read2,
    disco_reads,
    closest_bkp1,
    closest_bkp2,
    winSize,
    ftype,
    confidence,
    stopCod,
):
    with open(outFile, "+w") as out_file:
        for i in range(0, len(fGene)):
            peptides = []
            # Remove stop codon "*" from peptide sequences
            seq = pepSeq[i].replace("*", "")
            # Keep calm and roll
            if (
                "Fusion" == ftype[i]
            ):  # check the event type and use the appropriate rolling window function
                peps = rolling_window_infr(seq, winSize)
                for pep in peps:
                    if "|" in pep:  # get only peptides that span the junction
                        # remove low confidence peptides and peptides starting/ending with the junction point
                        if "?" in pep or pep.startswith("|") or pep.endswith("|"):
                            pass
                        else:
                            peptides.append(pep.replace("|", ""))
                if not peptides:  # check if there are no peptides of length <winsize>
                    pass
                else:
                    out_file.write(
                        "%s\t%s#%s#%s#%s#%s#%s#%s#%s#%s#%s#%s\n"
                        % (
                            fGene[i],
                            ("\t".join(peptides)),
                            ftype[i],
                            confidence[i],
                            stopCod[i],
                            bkpoint1[i],
                            bkpoint2[i],
                            split_read1[i],
                            split_read2[i],
                            disco_reads[i],
                            closest_bkp1[i],
                            closest_bkp2[i],
                        )
                    )
            if "Fusion-out-of-frame" == ftype[i]:
                start = seq.index("|")  # get the position of the junction
                if (
                    winSize > start
                ):  # check if there are less amino acids than the window size
                    seq2 = seq.replace("|", "")
                    peps = rolling_window(seq2, winSize)
                    for pep in peps:
                        if "?" in pep:
                            pass
                        else:
                            peptides.append(pep)
                    if not peptides:
                        pass
                    else:
                        out_file.write(
                            "%s\t%s#%s#%s#%s#%s#%s#%s#%s#%s#%s#%s\n"
                            % (
                                fGene[i],
                                ("\t".join(peptides)),
                                ftype[i],
                                confidence[i],
                                stopCod[i],
                                bkpoint1[i],
                                bkpoint2[i],
                                split_read1[i],
                                split_read2[i],
                                disco_reads[i],
                                closest_bkp1[i],
                                closest_bkp2[i],
                            )
                        )
                if (
                    winSize <= start
                ):  # check if the amino acids are more or equal to the window size
                    splitat = start - (
                        winSize - 1
                    )  # Cut the part of gene 1 that doesn't span the junction
                    seq3 = seq[splitat:].replace("|", "")
                    peps = rolling_window(seq3, winSize)
                    for pep in peps:
                        if "?" in pep:
                            pass
                        else:
                            peptides.append(pep)
                    if not peptides:
                        pass
                    else:
                        out_file.write(
                            "%s\t%s#%s#%s#%s#%s#%s#%s#%s#%s#%s#%s\n"
                            % (
                                fGene[i],
                                ("\t".join(peptides)),
                                ftype[i],
                                confidence[i],
                                stopCod[i],
                                bkpoint1[i],
                                bkpoint2[i],
                                split_read1[i],
                                split_read2[i],
                                disco_reads[i],
                                closest_bkp1[i],
                                closest_bkp2[i],
                            )
                        )
    out_file.close()


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        usage="cleave_peptides.py [-h] -i {/path/to/input/} -o {/path/to/input/} -p {8,9,10,11,15,16,17,18,19,20,21,22,23,24,25}",
        description="pepCleave v1.0",
    )
    parser.add_argument("-i", "--input", help="Input arriba file name", required=True)
    parser.add_argument("-o", "--output", help="Output file name", required=True)
    parser.add_argument(
        "-p",
        "--peptides",
        help="Peptide Length (Default = 8 and 15)",
        nargs="+",
        type=int,
        choices=range(8, 26),
        required=True,
    )
    args = parser.parse_args()

    # Get values / Declare lists
    inFile = args.input
    outFile = args.output
    pepLen = args.peptides
    ftype = []
    fusionGene = []
    bkpoint1 = []
    bkpoint2 = []
    split_read1 = []
    split_read2 = []
    disco_reads = []
    closest_bkp1 = []
    closest_bkp2 = []
    confidence = []
    stopCod = []
    pep_seq = []
    rows = []

    # Commence the onslaught
    ## Release the hount - remove duplicates, fetch the sequences, gene names, breakpoints and confidence levels
    sequence_hound(
        inFile,
        fusionGene,
        bkpoint1,
        bkpoint2,
        split_read1,
        split_read2,
        disco_reads,
        closest_bkp1,
        closest_bkp2,
        pep_seq,
        ftype,
        confidence,
        stopCod,
    )
    ## Cleave the fusion protein sequences at <user_defined> length
    for a in pepLen:
        out = outFile
        out += "_xeno_%s.tsv" % a
        window_onslaught(
            out,
            pep_seq,
            fusionGene,
            bkpoint1,
            bkpoint2,
            split_read1,
            split_read2,
            disco_reads,
            closest_bkp1,
            closest_bkp2,
            a,
            ftype,
            confidence,
            stopCod,
        )
