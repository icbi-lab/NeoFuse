#!/usr/bin/python
"""
Requirements:
    * Python >= 3.6.2

Copyright (c) 2020 Georgios Fotakis <georgios.fotakis@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = (
    "0",
    "1",
)
__version__ = ".".join(__version_info__)
__version__ += "-dev" if not RELEASE else ""

import argparse, sys

def get_hlas(in_file, hlas=[]):
    with open(in_file, "r") as ifile:
        for line in ifile:
            hlas.append(line)
    
    return hlas

def filter_class_I(hlas=[]):
    class_I_genes = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-K", "HLA-L", "HLA-V"]
    class_I_types = ["A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "V"]
    class_II = []
    for hla in list(hlas):
        hla_gene = hla.split("*")[0]
        if hla_gene in class_I_genes:
            continue
        elif hla_gene not in class_I_genes and hla_gene in class_I_types:
            continue
        else:
            class_II.append(hla)

    return class_II

def filter_class_II(hlas=[]):
    class_I_genes = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-K", "HLA-L", "HLA-V"]
    class_I_types = ["A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "V"]
    hla_I = []
    for hla in list(hlas):
        hla_gene = hla.split("*")[0]
        if hla_gene in class_I_genes:
            hla_I.append(hla.split("-")[1])
        elif hla_gene not in class_I_genes and hla_gene in class_I_types:
            hla_I.append(hla)
        else:
            continue

    return hla_I

def out(out_file, hlas):
    hlas = ''.join(hlas)
    if hlas.endswith("\n"):
        pass
    else:
        hlas = hlas + "\n"
    with open(out_file, "w") as ofile:
        ofile.write(hlas)

if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description="Calculate CSiN")

    parser.add_argument(
        "--custom_list", 
        required=False, 
        help="HLAs custom list input"
    )
    parser.add_argument(
        "--output_dir",
        required=True, 
        help="Path to the output directory"
    )
    parser.add_argument(
        "--sample_name", 
        required=True, 
        help="Sample name"
    )
    args = parser.parse_args()
    # Parse arguments
    in_file = args.custom_list
    out_file_i = args.output_dir + args.sample_name + "_HLA_Optitype.txt"
    out_file_ii = args.output_dir + args.sample_name + "_custom_HLA_II.txt"


    hlas = []
    get_hlas(in_file, hlas)
    hla_i = filter_class_II(hlas)
    hla_ii = filter_class_I(hlas)
    out(out_file_i, hla_i)
    out(out_file_ii, hla_ii)
