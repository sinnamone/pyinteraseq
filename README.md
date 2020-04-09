# InteractomSeq

### A web server for the identification and profiling of domains and epitopes from Phage Display and Next Generation Sequencing data.

InteractomeSeq: a web server for the identification and profiling of domains and epitopes from Phage Display and Next Generation Sequencing data. InteractomeSeq is a webtool allowing either genomic or single gene domainome analysis of phage libraries generated and selected by following the interactome-sequencing approach.

InteractomeSeq data analysis workflow is composed of four sequential steps that, starting from raw sequencing reads, generate the list of putative domains with genomic annotations. In the first step, InteractomeSeq checks if the input files (raw reads, reference genome sequence, annotation list) are properly formatted. In the second step low-quality sequencing data are first trimmed and then discarded basing on minimal length requirements. In the third step, the remaining reads are aligned to the reference genome, a SAM file is generated and only reads with high quality score are processed and converted into a BAM file. After alignment, InteractomeSeq performs the step for domains detection and generates a list of CDS portions classified as putative domains/epitopes, the list is ranked taking in consideration the focus which is an index obtained from the ratio between maximal depth and coverage of a specific genic portion.

# 1)  Installation 

With this steps below users can mirror the excecution of the webserver https://interactomeseq.ba.itb.cnr.it/ in order to perform analysis from command line. 

### Install Miniconda

Miniconda is a Python distribution, package manager, and virtual environment solution. While QIIME 1 is Python 2 software, we recommend installing Miniconda with Python 3 (miniconda3), as many bioinformatics packages are now transitioning to Python 3. You can still install Python 2 software with miniconda3 by passing the python=2.7 flag when you create a new environment; otherwise the default Python version will be Python 3.

Begin by downloading Miniconda and following the associated installation instructions.

https://docs.conda.io/en/latest/miniconda.html

### Create your InteractomeSeq environment and install the dependences

InteractomeSeq consists of native Python 2.7 code and additionally wraps many external applications. These instructions describe how to perform a base installation of QIIME using Miniconda.


```bash
conda create -n interactomeseq python=2.7 matplotlib=1.5.3 trimmomatic=0.36 kallisto=0.46.0 cutadapt=1.12 blast r-base=3.4.1 pandas=0.21.0 samtools=0.1.19  cython biom-format biopython pybedtools bioconductor-edger ucsc-bedgraphtobigwig -c conda-forge -c bioconda -c r -c anaconda
```

### Activate your InteractomeSeq environment and test your installation


```bash
conda activate interactomeseq
```


```bash
conda list | less 
```

If you decide later that you don’t want the environment or its packages anymore, deactivate the environment and then run this command:


```bash
conda remove --name interactomeseq --all
```

### Installing additional dependencies accessible in your interactomeseq environment


```bash
pip install https://pypi.python.org/packages/source/m/matplotlib/matplotlib-1.4.3.tar.gz
```

Then install qiime1 version 1.91.1


```bash
conda install -c bioconda pysam
```


```bash
conda install -c conda-forge readline=6.2
```


```bash
pip install qiime
```

### Clone the Github repository


```bash
git clone https://github.com/sinnamone/InteractomeSeq
```
