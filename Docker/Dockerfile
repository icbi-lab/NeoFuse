# work from ubuntu 18.04 base image
FROM ubuntu:18.04

# set the environment variables
ENV optitype_version 1.3.2
ENV samtools_version 1.9
ENV bedtools_version 2.27.1
ENV PATH /usr/local/bin/:${PATH}

# run update and install necessary utils
RUN apt-get update -y && apt-get install -y \
	build-essential \
	curl \
	unzip \
	python-minimal \
	bzip2 \
	zlib1g-dev \
	libncurses5-dev \
	libncursesw5-dev \
	libnss-sss \
	libbz2-dev \
	liblzma-dev \
	libhdf5-dev \
	glpk-utils \
	python-pip \
	libpng-dev \
	libfreetype6-dev \
	libfreetype6 \
	pkg-config \
	vim \
	less \
	coinor-cbc \
	wget \
	gcc \
	g++ \
	make \
	zlib1g-dev \
	unzip \
	git \
	cmake \
	autoconf \
	pigz \
	python3-pip \
	coreutils \
	git-lfs \
	tcsh \
	gawk \
	locales \
	vim \
	procps

# Fix locale issues
RUN export LC_ALL=en_US.UTF-8 && \
	export LANG=en_US.UTF-8 && \
	locale-gen en_US.UTF-8

# Set working directory
WORKDIR /usr/local/bin/

#Install STAR
ARG STAR_VERSION=2.7.1a
RUN set -ex && \
	wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.zip && \
	unzip ${STAR_VERSION}.zip && \
	cd STAR-${STAR_VERSION}/source && \
	make STARstatic && \
	cp STAR /usr/local/bin && \
	cd /usr/local/bin && \
	'rm' -rf STAR-${STAR_VERSION} && \
	rm -rf ${STAR_VERSION}.zip

# Install YARA
RUN git clone https://github.com/seqan/seqan.git && \
	mkdir yara-build
WORKDIR /usr/local/bin/yara-build/
RUN cmake ../seqan -DSEQAN_BUILD_SYSTEM=APP:yara && \
	make all && \
	ln -s /usr/local/bin/yara-build/bin/yara* /usr/local/bin/
WORKDIR /usr/local/bin/

# Install samtools
ADD https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 /usr/local/bin/

RUN tar -xjf /usr/local/bin/samtools-${samtools_version}.tar.bz2 -C /usr/local/bin/ && \
	cd /usr/local/bin/samtools-${samtools_version}/ && make && \
	cd /usr/local/bin/samtools-${samtools_version}/ && make install

# Install required python modules
RUN python -m pip install --upgrade pip && \
	python3 -m pip install --upgrade pip && \
	python -m pip install 'NumPy==1.16.6' && \
	python -m pip install 'numexpr==2.7.0' && \
	python -m pip install 'six==1.16.0' && \
	python -m pip install 'tables==3.5.2' && \
	python -m pip install 'Pyomo' && \
	python -m pip install 'Pandas' && \
	python -m pip install 'Pysam' && \
	python -m pip install 'Matplotlib' && \
	python -m pip install 'Future' && \
	python3 -m pip install 'keras==2.3' && \
	python3 -m pip install 'setuptools==47.3.1' && \
	python3 -m pip install 'Pandas' && \
	python3 -m pip install 'six' && \
	python3 -m pip install 'pandas>=0.20.3' && \
	python3 -m pip install 'tensorflow>=2.2.0,<2.3.0' && \
	python3 -m pip install 'appdirs' && \
	python3 -m pip install 'scikit-learn' && \
	python3 -m pip install 'mhcgnomes' && \
	python3 -m pip install 'pyyaml' && \
	python3 -m pip install 'tqdm' && \
	python3 -m pip install 'np_utils'

# Install python modules
# RUN pip install 'NumPy==1.15.4' && \
# 	pip install 'Pyomo==4.2.10784' && \
# 	pip install 'Pandas==0.16.2' && \
# 	pip install 'Pysam==0.8.3' && \
# 	pip install 'Matplotlib==1.4.3' && \
# 	pip install 'Future==0.15.2' && \
# 	pip install 'tables==3.2.2'

# Install OptiType
RUN curl -SL https://github.com/FRED-2/OptiType/archive/v${optitype_version}.zip \
  >  v${optitype_version}.zip && \
	unzip v${optitype_version}.zip

## set up default configuration
RUN mv /usr/local/bin/OptiType-${optitype_version}/config.ini.example /usr/local/bin/OptiType-${optitype_version}/config.ini && \
	sed -i 's/\/path\/to\//\/usr\/local\/bin\//' /usr/local/bin/OptiType-${optitype_version}/config.ini && \
	sed -i 's/threads=16/threads=8/g' /usr/local/bin/OptiType-${optitype_version}/config.ini && \
	sed -i 's/glpk/cbc/' /usr/local/bin/OptiType-${optitype_version}/config.ini

## Create YARA indices from OptiType HLA reference
RUN mkdir /usr/local/bin/yara_idx/ && \
	yara_indexer /usr/local/bin/OptiType-1.3.2/data/hla_reference_rna.fasta -o /usr/local/bin/yara_idx/hla_reference_rna

# Install R and packages 
RUN export DEBIAN_FRONTEND=noninteractive && \
	apt-get update -y && \
	apt-get install -y --no-install-recommends r-base ca-certificates libcurl4-openssl-dev libxml2-dev && \
	Rscript -e 'install.packages("circlize", repos="http://cran.r-project.org"); source("https://bioconductor.org/biocLite.R"); biocLite(c("GenomicRanges", "GenomicAlignments"))'

# Install arriba 2.1.0 verion
RUN mkdir -p /usr/local/bin/arriba_v2
WORKDIR /usr/local/bin/arriba_v2
RUN git clone https://github.com/suhrig/arriba.git
WORKDIR /usr/local/bin/arriba_v2/arriba
RUN git checkout develop \
	&& make \
	&& ln -s /usr/local/bin/arriba_v2/arriba/arriba /usr/local/bin/arriba
WORKDIR /usr/local/bin/

# # Install arriba
# RUN wget -q -O - "https://github.com/suhrig/arriba/releases/download/v1.1.0/arriba_v1.1.0.tar.gz" | tar -xzf -

# ## make wrapper script for download_references.sh
# RUN echo '#!/bin/bash\n\
# cd /references\n\
# /arriba*/download_references.sh $1 $2 && \\\n\
# cp /arriba*/database/*${1%+*}* /references' > /usr/local/bin/download_references.sh && \
# 	chmod a+x /usr/local/bin/download_references.sh

# ## make wrapper script for run_arriba.sh
# RUN echo '#!/bin/bash\n\
# cd /output\n\
# /arriba*/run_arriba.sh /references/STAR_index_* /references/*.gtf /references/*.fa /references/blacklist_*.tsv.gz /read1.fastq.gz /read2.fastq.gz ${1-8}' > /usr/local/bin/arriba.sh && \
# chmod a+x /usr/local/bin/arriba.sh && \
# 	ln -s /usr/local/bin/arriba_v1.1.0/arriba /usr/local/bin/arriba

# ## make wrapper script for draw_fusions.R
# RUN echo '#!/bin/bash\n\
# Rscript /arriba*/draw_fusions.R --annotation=$(ls /references/*.gtf) --fusions=/fusions.tsv --output=/output/fusions.pdf --proteinDomains=$(ls /references/protein_domains_*.gff3) --alignments=/Aligned.sortedByCoord.out.bam --cytobands=$(ls /references/cytobands_*.tsv)' > /usr/local/bin/draw_fusions.sh && \
# 	chmod a+x /usr/local/bin/draw_fusions.sh

# Install featureCounts
ARG PACKAGE_VERSION=1.6.1
RUN wget -q https://downloads.sourceforge.net/project/subread/subread-${PACKAGE_VERSION}/subread-${PACKAGE_VERSION}-source.tar.gz && \
	tar -xzf subread-${PACKAGE_VERSION}-source.tar.gz && \
	cd subread-${PACKAGE_VERSION}-source/src && \
	make -f Makefile.Linux && \
	cd ../bin && \
	mv utilities/* . && \
	rmdir utilities && \
	ln -s /usr/local/bin/subread-1.6.1-source/bin/featureCounts /usr/local/bin/featureCounts

# Copy reference files and scripts
RUN mkdir -p /usr/local/bin/genomes/ && \
	mkdir -p /usr/local/bin/source/ && \
	mkdir -p /usr/local/bin/tmp/
COPY ./source/blacklist_hg38_GRCh38_2018-11-04.tsv.gz /usr/local/bin/genomes/
COPY ./source/ /usr/local/bin/source/
WORKDIR /usr/local/bin/source/
RUN chmod +x NeoFuse_single.sh && \
	chmod +x NeoFuse_multi.sh && \
	chmod +x NeoFuse.sh && \
	chmod +x NeoFuse_download.sh && \
	ln -s /usr/local/bin/source/NeoFuse_single.sh /usr/local/bin/NeoFuse_single && \
	ln -s /usr/local/bin/source/NeoFuse_multi.sh /usr/local/bin/NeoFuse_multi && \
	ln -s /usr/local/bin/source/NeoFuse.sh /usr/local/bin/NeoFuse && \
	ln -s /usr/local/bin/source/NeoFuse_download.sh /usr/local/bin/NeoFuse-download

# Clean up
RUN apt-get autoclean -y && apt-get clean -y && \
	rm /usr/local/bin/samtools-1.9.tar.bz2 && \
	rm /usr/local/bin/subread-1.6.1-source.tar.gz && \
	rm /usr/local/bin/v1.3.2.zip

WORKDIR /usr/local/bin/

# Install MHCFlurry
RUN pip install mhcflurry

ENV MHCFLURRY_DEFAULT_CLASS1_MODELS /home/neofuse/.local/share/mhcflurry/4/2.0.0
# ENV MHCFLURRY_DEFAULT_CLASS1_MODELS /home/neofuse/.local/share/mhcflurry/4/1.6.0
# ENV MHCFLURRY_DOWNLOADS_DIR /home/neofuse/.local/share/mhcflurry/4/1.4.0

# Add user to install MHCflurry models in "~/" dir
RUN useradd -ms /bin/bash neofuse
USER neofuse

# Download MHCflurry class I pan models
RUN mhcflurry-downloads fetch models_class1_pan

WORKDIR /home/
