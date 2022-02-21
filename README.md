# NeoFuse
NeoFuse is a user-friendly pipeline for the prediction of fusion neoantigens from tumor RNA-seq data.

NeoFuse takes single-sample FASTQ files of RNA-seq reads (single- or paired-end) as input and predicts putative fusion neoantigens through five main analytical modules based on state-of-the-art computational tools: 

* Genotyping of class-I Human Leukocyte Antigen (HLA) genes at 4-digit resolution using OptiType (Szolek et al., 2014).
* Prediction of fusion peptides using Arriba (https://github.com/suhrig/arriba), together with confidence scores reflecting the likelihood that a fusion is caused by a tumor-specific genomic rearrangement and is not due to technical artifacts.
* Prediction of the binding affinity of fusion peptides to HLA types, quantified as half maximal inhibitory concentration (IC50) and percentile rank, using MHCflurry (O’Donnell et al., 2018) or netMHCpan (Jurtz et al., 2017).
* Quantification of gene expression levels, as transcripts per million (TPM), using STAR (Dobin et al., 2013) and featureCounts (Liao et al., 2014).  
* Neoantigen prioritization based on IC50 binding affinity and confidence score, and annotation of each neoantigen with: IC50, percentile rank, confidence score, binding HLA type, expression of the fusion and HLA genes in TPM, and information about the presence of a premature stop codon that might cause nonsense mediated decay of the fusion transcript.

We advise using paired-end data to increase sensitivity and accuracy of gene fusion detection.

## Requirements
* At least 30 GB of RAM
* Docker  (version 19.03 or later) or
* Singularity (version 3.0 or later)

## 1. Installation

NeoFuse can be installed through the following four steps.

### 1.1. Install Docker or Singularity engine locally
[Instructions for Docker installation](https://docs.docker.com/install/)

[Instructions for Singularity installation](https://sylabs.io/docs/)

### 1.2. Clone the NeoFuse repo locally
```
$ git clone https://github.com/icbi-lab/NeoFuse.git
```

Or you can find the script freely available [here](https://icbi.i-med.ac.at/software/NeoFuse/downloads/NeoFuse-v1.1.1.zip) (deprecated).

Add NeoFuse to PATH:

```
$ export PATH=$PATH:~/path/to/NeoFuse
```

### 1.3. Pull/build the NeoFuse image
The NeoFuse image can be automatically generated using the NeoFuse script:

**Docker:**
```
$ NeoFuse -B --docker
```

**Singularity:**
```
$ NeoFuse -B --singularity
```

### 1.4. Download reference genomes and build STAR indexes
The NeoFuse script can be also used to generate the genomes and indexes required by the analysis:

**Docker:**
```
$ NeoFuse -R -o </path/to/output_folder> -n [cores] -V [genome version] --docker
```

**Singularity:**
```
$ NeoFuse -R -o </path/to/output_folder> -n [cores] -V [genome version] --singularity
```

**\<Arguments\>**

-o: Output directory

**[Options]**

-n: Number of cores (default: 1)

-V: Genome version, either “GRCh37” and “GRCh38” (default: GRCh38)

**Note:** this process may take more than 1 hour, depending on the internet connection and the processing power.

## 2. Usage

**Notes**

* On Mac OS X, you need to have Docker running to execute NeoFuse.

* On Mac OS X, you might need to increase CPUs, Memory and Swap in Docker settings (Settings > Preferences > Advanced). 

### 2.1. Analysis of single samples

NeoFuse can process single samples with the following command:

```
$ NeoFuse <arguments> [options] --singularity (or --docker) 
```

**\<Arguments\>**

**-1:** Path to read 1 FASTQ file (mandatory)

**-2:** Path to read 2 FASTQ fie (optional for single-end reads)

**-s:** Path to STAR index directory (mandatory)

**-g:** Path to reference genome FASTA file (mandatory)

**-a:** Path to annotation GTF file (mandatory)

**Note**: All input files passed as arguments must be unzipped.

**[Options]**

**-d:** Run ID (the name of the output files)

**-m:** Minimum peptide length (values: 8, 9, 10, or 11; default: 8)

**-M:** Maximum peptide length (values: 8, 9, 10, or 11; default: 8) *

**-n:** Number of cores (default: 1)

**-t:** IC50 binding affinity threshold (default: 500)

**-T:** Percentile rank threshold (default: Inf)

**-c:** Mimimum confidence score (values: H, M, or L; default: L) **

**-C:** Custom HLA I list (TXT file containing custom HLA I types)

**-k:** Keep STAR output (BAM files)

**-l:** int>=0:  maximum available RAM (bytes) for sorting BAM

**-K:** File containing known/recurrent fusions (see [Arriba manual](https://arriba.readthedocs.io/en/latest/input-files/#known-fusions) for more details)

**-f:** Comma separated list of Arriba filters to disable (do not use space separeted lists, see [Arriba manual](https://arriba.readthedocs.io/en/latest/command-line-options/) for more details)

**-v:** VCF or Tab-separated file with coordinates of structural variants found using whole-genome sequencing data (The file may be gzip-compressed, for more info refer to [Arriba manual](https://arriba.readthedocs.io/en/latest/input-files/#structural-variant-calls-from-wgs))

**-S:** Determines how far a genomic breakpoint may be away from a transcriptomic breakpoint to still consider it as a related event (used with -v parameter - Default = 100000)

**--singularity:** NeoFuse will use the Singularity image

**--docker:** NeoFuse will use the Docker image

\* NeoFuse will compute the binding affinity for all the possible lengths of peptides between the minimum and maximum input. For example if a user specifies '-m 8' and '-M 11', NeoFuse will compute the binding affinity for all peptides of length 8, 9, 10, and 11. To consider just one specific length, use **only** the '-m' argument.

** The [minimum Arriba confidence score](https://arriba.readthedocs.io/en/latest/interpretation-of-results/) can be set to: H (to return only high confidence fusions), M (for high and medium confidence fusions), or L (for high, medium, and low confidence fusions).

### 2.2. Analysis of multiple samples

For multiple-sample analysis, a TSV input file reporting the sample identifiers and path to input files has to be prepared. Format:

**Paired-end reads:**
```
#ID	Read1	Read2
Sample1 /path/to/Sample1_read_1.fastq	/path/to/Sample1_read_2.fastq
Sample2 /path/to/Sample2_read_1.fastq	/path/to/Sample2_read_2.fastq


```
**Single-end reads:**
```
#ID	Read1
Sample1	/path/to/Sample1_read_1.fastq
Sample2	/path/to/Sample2_read_1.fastq


```
**Notes:** The first line of the TSV should start with a hashtag. FASTQ files should be in the same directory. No whitespaces allowed.

Once the TSV file is created, the samples can be analyzed with the following command:


```
$ NeoFuse <arguments> [options] --singularity (or --docker)
```

**\<Arguments\>**

**-i:** Path to the input TSV file (mandatory)

**-s:** Path to STAR index directory (mandatory)

**-g:** Path to reference genome FASTA file (mandatory)

**-a:** Path to annotation GTF file (mandatory)

**Note**: All input files passed as arguments must be unzipped.

**[Options]**

**-m:** Minimum peptide length (values: 8, 9, 10, or 11; default: 8)

**-M:** Maximum peptide length (values: 8, 9, 10, or 11; default: 8) *

**-n:** Number of cores (default: 1)

**-t:** IC50 binding affinity threshold (default: 500)

**-T:** Percentile rank threshold (default: Inf)

**-c:** Mimimum confidence score (values: H, M, or L; default: L) **

**-C:** Custom HLA I list (TXT file containing custom HLA I types)

**-k:** Keep STAR output (BAM files)

**-l:** int>=0:  maximum available RAM (bytes) for sorting BAM

**-K:** File containing known/recurrent fusions (see [Arriba manual](https://arriba.readthedocs.io/en/latest/input-files/#known-fusions) for more details)

**-f:** Comma separated list of Arriba filters to disable (do not use space separeted lists, see [Arriba manual](https://arriba.readthedocs.io/en/latest/command-line-options/) for more details)

**-v:** VCF or Tab-separated file with coordinates of structural variants found using whole-genome sequencing data (The file may be gzip-compressed, for more info refer to [Arriba manual](https://arriba.readthedocs.io/en/latest/input-files/#structural-variant-calls-from-wgs))

**-S:** Determines how far a genomic breakpoint may be away from a transcriptomic breakpoint to still consider it as a related event (used with -v parameter - Default = 100000)

**--singularity:** NeoFuse will use the Singularity image

**--docker:** NeoFuse will use the Docker image

\* NeoFuse will compute the binding affinity for all the possible lengths of peptides between the minimum and maximum input. For example if a user specifies '-m 8' and '-M 11', NeoFuse will compute the binding affinity for all peptides of length 8, 9, 10, and 11. To consider just one specific length, use **only** the '-m' argument.

** The [mimimum Arriba confidence score](https://arriba.readthedocs.io/en/latest/interpretation-of-results/) can be set to: H (to return only high confidence fusions), M (for high and medium confidence fusions), or L (for high, medium, and low confidence fusions).

### 2.3. Binding affinity prediction with netMHCpan
Due to license compatibility issues, netMHCpan is fully integrated but **not distributed** as part of NeoFuse. 

If there is an existing local [installation](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan) of netMHCpan, peptide-HLA binding affinity (IC50 and rank) can be predicted with netMHCpan instead of MHCflurry using the following command:

```
$ NeoFuse <arguments> [options] -N [/path/to/netMHCpan_direcotry] --singularity (or --docker)
```

## 3. Results

### 3.1. Main output directory
NeoFuse will create an output directory with the following structure:

```
/NeoFuse/output/directory/
├── Sample1
│   ├── Arriba
│   ├── LOGS
│   ├── NeoFuse
│   ├── OptiType
│   └── TPM
├── Sample2
│   ├── Arriba
│   ├── LOGS
│   ├── NeoFuse
│   ├── OptiType
│   └── TPM
…
└── SampleN
    ├── Arriba
    ├── LOGS
    ├── NeoFuse
    ├── OptiType
    └── TPM
```

### 3.2. Output subdirectories
#### 3.2.1. Arriba

Sample.fusions.tsv file contains a list of gene fusions sorted from highest to lowest confidence. 

Sample.fusions.discarded.tsv contains all events that Arriba classified as artifacts or that are also observed in healthy tissues.


```
/Arriba
├── Sample1.fusions.discarded.tsv
└── Sample1.fusions.tsv
```

#### 3.2.2. LOGS
The standard output (sdout and stderr) for every tool used in the run is stored in the LOGS directory.
File names may differ depending on the tools, peptide length, etc.

```
/LOGS
├── Sample1_10_MHCFlurry.log
├── Sample1_11_MHCFlurry.log
├── Sample1_8_MHCFlurry.log
├── Sample1_9_MHCFlurry.log
├── Sample1.arriba.err
├── Sample1.arriba.log
├── Sample1.cleave.log
├── Sample1.counts_to_tpm.log
├── Sample1.featureCounts.log
├── Sample1.final.log
├── Sample1.Log.final.out
├── Sample1.Log.out
├── Sample1.Log.std.out
├── Sample1.optitype.log
├── Sample1.samtools.err
├── Sample1.samtools.log
├── Sample1.STAR.err
├── Sample1.STAR.log
└── Sample1.association.log
```
#### 3.2.3. OptiType
HLA_Optitype.txt contains the HLA types predicted by OptiType

```
/OptiType
├── Sample1_coverage_plot.pdf
└── Sample1_HLA_Optitype.txt
```

#### 3.2.4. TPM
Contains all TPM expression values for all the genes

```
/TPM
└── Sample1.tpm.txt
```

#### 3.2.5. NeoFuse
Contains the **final output** of the pipeline, which consists of three files:

```
/NeoFuse
├── Sample1_filtered.tsv
└── Sample1_unfiltered.tsv
```

Sample_unfiltered.tsv contains all the predicted fusion peptides and their annotations.

Sample_filtered.tsv contains a list of putative fusion neoantigens (selected considering the user-defined IC50/rank and confidence score thresholds).
This file reports for each putative neoantigen: confidence score, binding HLA type, expression of the fusion and HLA genes in TPM, and information about the presence of a premature stop codon that might cause nonsense mediated decay of the fusion transcript.
Example format:

```
Fusion	Gene1	Gene2	Breakpoint1	Breakpoint2	Split_Reads1	Split_Reads2	Discordant_Reads	Closest_Breakpoint1	Closest_Breakpoint2	HLA_Type	Fusion_Peptide	IC50	Rank	Event_Type	Stop_Codon	Confidence	Gene1_TPM	Gene2_TPM	Avg_TPM	HLA_TPM
SETD2-ELP6	SETD2	ELP6	chr3:47056821	chr3:47504448	11	19	30	chr3:47053118(3703)	chr3:47053118(3703)	HLA-B*35:02	TPPIVQGVSL	337.8201495211543	0.16212499999999996	Fusion	no	high	37.73	34.87	36.30	1296.32
SETD2-ELP6	SETD2	ELP6	chr3:47056821	chr3:47504448	11	19	30	chr3:47053118(3703)	chr3:47053118(3703)	HLA-B*35:08	TPPIVQGVSL	174.1196814867352	0.4007499999999997	Fusion	no	high	37.73	34.87	36.30	1296.32
NPAS3-ABHD4	NPAS3	ABHD4	chr14:33215587	chr14:22603390	2	2	32	.	.	HLA-B*35:08	TPCEGLQNKF	71.80709202149346	0.13862500000000008	Fusion	no	high	3.45	136.77	70.11	1296.32
```

## 4. References
Dobin,A. et al. (2013) STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29, 15–21.

Jurtz, V. et al. (2017) NetMHCpan-4.0: Improved Peptide–MHC Class I Interaction Predictions Integrating Eluted Ligand and Peptide Binding Affinity Data. J. Immunol., 199, 3360-3368.

Liao,Y. et al. (2014) featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30, 923–930.

O’Donnell,T.J. et al. (2018) MHCflurry: Open-Source Class I MHC Binding Affinity Prediction. Cell Syst, 7, 129–132.e4.

Szolek,A. et al. (2014) OptiType: precision HLA typing from next-generation sequencing data. Bioinformatics, 30, 3310–3316.
