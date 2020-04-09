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
which conda
```

    /Users/simone/miniconda3/bin/conda



```bash
conda create -n interactomeseq python=2.7 matplotlib=1.5.3 trimmomatic=0.36 kallisto=0.46.0 cutadapt=1.12 blast r-base=3.4.1 pandas=0.21.0 samtools=0.1.19  cython biom-format biopython pybedtools bioconductor-edger ucsc-bedgraphtobigwig -c conda-forge -c bioconda -c r -c anaconda -y
```

    Collecting package metadata (current_repodata.json): done
    Solving environment: failed with repodata from current_repodata.json, will retry with next repodata source.
    Collecting package metadata (repodata.json): done
    Solving environment: done
    
    ## Package Plan ##
    
      environment location: /Users/simone/.conda/envs/interactomeseq
    
      added / updated specs:
        - bioconductor-edger
        - biom-format
        - biopython
        - blast
        - cutadapt=1.12
        - cython
        - kallisto=0.46.0
        - matplotlib=1.5.3
        - pandas=0.21.0
        - pybedtools
        - python=2.7
        - r-base=3.4.1
        - samtools=0.1.19
        - trimmomatic=0.36
        - ucsc-bedgraphtobigwig
    
    
    The following packages will be downloaded:
    
        package                    |            build
        ---------------------------|-----------------
        blast-2.7.1                |      boost1.64_3        16.0 MB  bioconda
        boost-1.64.0               |           py27_4         307 KB  conda-forge
        boost-cpp-1.64.0           |                1        16.2 MB  conda-forge
        cairo-1.14.8               |                0         572 KB  anaconda
        curl-7.69.1                |       h2d98d24_0         131 KB  conda-forge
        fontconfig-2.12.1          |                3         176 KB  anaconda
        freetype-2.5.5             |                2         738 KB  anaconda
        gettext-0.19.8.1           |       h1f1d5ed_1         3.4 MB  conda-forge
        glib-2.50.2                |                1         2.9 MB  anaconda
        gnutls-3.5.19              |       h2a4e5f8_1         1.9 MB  conda-forge
        harfbuzz-0.9.39            |                2         302 KB  anaconda
        krb5-1.17.1                |       h1752a42_0         1.2 MB  conda-forge
        libcurl-7.69.1             |       hc0b9707_0         535 KB  conda-forge
        libiconv-1.14              |                4         1.3 MB  conda-forge
        libxml2-2.9.5              |                1         1.9 MB  conda-forge
        matplotlib-1.5.3           |      np111py27_1         6.0 MB  anaconda
        nettle-3.3                 |                0         869 KB  conda-forge
        pandas-0.21.0              |           py27_0        10.8 MB  conda-forge
        pango-1.40.3               |                1         448 KB  anaconda
        perl-5.22.0.1              |                0        14.6 MB  conda-forge
        perl-archive-tar-2.18      |                1          27 KB  bioconda
        perl-list-moreutils-0.15   |       pl5.22.0_0          37 KB  bioconda
        pyqt-5.6.0                 |py27hc26a216_1008         4.4 MB  conda-forge
        qt-5.6.3                   |       h1d42b2f_0        67.0 MB  anaconda
        samtools-0.1.19            |       hdd8ed8b_4         1.8 MB  bioconda
        sip-4.18.1                 |py27h0a44026_1000         243 KB  conda-forge
        sqlite-3.26.0              |       ha441bb4_0         2.3 MB  anaconda
        ucsc-bedgraphtobigwig-366  |       h1341992_0         392 KB  bioconda
        ------------------------------------------------------------
                                               Total:       156.5 MB
    
    The following NEW packages will be INSTALLED:
    
      _r-mutex           conda-forge/noarch::_r-mutex-1.0.1-anacondar_1
      bedtools           bioconda/osx-64::bedtools-2.29.2-h37cfd92_0
      bioconductor-edger bioconda/osx-64::bioconductor-edger-3.22.5-r341hfc679d8_0
      bioconductor-limma bioconda/osx-64::bioconductor-limma-3.36.5-r341h470a237_0
      biom-format        conda-forge/osx-64::biom-format-2.1.7-py27h917ab60_1002
      biopython          conda-forge/osx-64::biopython-1.76-py27h0b31af3_0
      blast              bioconda/osx-64::blast-2.7.1-boost1.64_3
      boost              conda-forge/osx-64::boost-1.64.0-py27_4
      boost-cpp          conda-forge/osx-64::boost-cpp-1.64.0-1
      bz2file            conda-forge/noarch::bz2file-0.98-py_0
      bzip2              conda-forge/osx-64::bzip2-1.0.8-h0b31af3_2
      ca-certificates    conda-forge/osx-64::ca-certificates-2020.4.5.1-hecc5488_0
      cairo              anaconda/osx-64::cairo-1.14.8-0
      certifi            conda-forge/osx-64::certifi-2019.11.28-py27h8c360ce_1
      click              conda-forge/noarch::click-7.1.1-pyh8c360ce_0
      curl               conda-forge/osx-64::curl-7.69.1-h2d98d24_0
      cutadapt           bioconda/osx-64::cutadapt-1.12-py27_1
      cycler             conda-forge/noarch::cycler-0.10.0-py_2
      cython             conda-forge/osx-64::cython-0.29.15-py27hde1d37b_1
      fontconfig         anaconda/osx-64::fontconfig-2.12.1-3
      freetype           anaconda/osx-64::freetype-2.5.5-2
      future             conda-forge/osx-64::future-0.18.2-py27h8c360ce_1
      gettext            conda-forge/osx-64::gettext-0.19.8.1-h1f1d5ed_1
      glib               anaconda/osx-64::glib-2.50.2-1
      gmp                conda-forge/osx-64::gmp-6.1.2-h0a44026_1000
      gnutls             conda-forge/osx-64::gnutls-3.5.19-h2a4e5f8_1
      gsl                conda-forge/osx-64::gsl-2.6-ha2d443c_0
      harfbuzz           anaconda/osx-64::harfbuzz-0.9.39-2
      hdf5               conda-forge/osx-64::hdf5-1.10.5-nompi_h0cbb7df_1103
      icu                conda-forge/osx-64::icu-58.2-h0a44026_1000
      jpeg               conda-forge/osx-64::jpeg-9c-h1de35cc_1001
      kallisto           bioconda/osx-64::kallisto-0.46.0-hc0c33ce_1
      krb5               conda-forge/osx-64::krb5-1.17.1-h1752a42_0
      libblas            conda-forge/osx-64::libblas-3.8.0-14_openblas
      libcblas           conda-forge/osx-64::libcblas-3.8.0-14_openblas
      libcurl            conda-forge/osx-64::libcurl-7.69.1-hc0b9707_0
      libcxx             conda-forge/osx-64::libcxx-10.0.0-0
      libdeflate         conda-forge/osx-64::libdeflate-1.5-h01d97ff_0
      libedit            conda-forge/osx-64::libedit-3.1.20170329-0
      libffi             conda-forge/osx-64::libffi-3.2.1-h4a8c4bd_1007
      libgcc             conda-forge/osx-64::libgcc-4.8.5-hdbeacc1_10
      libgfortran        conda-forge/osx-64::libgfortran-3.0.1-0
      libiconv           conda-forge/osx-64::libiconv-1.14-4
      liblapack          conda-forge/osx-64::liblapack-3.8.0-14_openblas
      libopenblas        conda-forge/osx-64::libopenblas-0.3.7-hd44dcd8_1
      libpng             conda-forge/osx-64::libpng-1.6.34-ha92aebf_2
      libssh2            conda-forge/osx-64::libssh2-1.8.2-hcdc9a53_2
      libtiff            conda-forge/osx-64::libtiff-4.0.9-h79f4b77_1002
      libxml2            conda-forge/osx-64::libxml2-2.9.5-1
      matplotlib         anaconda/osx-64::matplotlib-1.5.3-np111py27_1
      ncurses            conda-forge/osx-64::ncurses-5.9-10
      nettle             conda-forge/osx-64::nettle-3.3-0
      numpy              conda-forge/osx-64::numpy-1.11.3-py27hdf140aa_1207
      openjdk            conda-forge/osx-64::openjdk-11.0.1-hbbe82c9_1018
      openssl            conda-forge/osx-64::openssl-1.1.1f-h0b31af3_0
      pandas             conda-forge/osx-64::pandas-0.21.0-py27_0
      pango              anaconda/osx-64::pango-1.40.3-1
      pcre               conda-forge/osx-64::pcre-8.39-0
      perl               conda-forge/osx-64::perl-5.22.0.1-0
      perl-archive-tar   bioconda/osx-64::perl-archive-tar-2.18-1
      perl-list-moreuti~ bioconda/osx-64::perl-list-moreutils-0.15-pl5.22.0_0
      pip                conda-forge/noarch::pip-20.0.2-py_2
      pixman             conda-forge/osx-64::pixman-0.34.0-h1de35cc_1003
      pybedtools         bioconda/osx-64::pybedtools-0.8.1-py27h2dec4b4_0
      pyparsing          conda-forge/noarch::pyparsing-2.4.7-pyh9f0ad1d_0
      pyqi               conda-forge/noarch::pyqi-0.3.2-py_0
      pyqt               conda-forge/osx-64::pyqt-5.6.0-py27hc26a216_1008
      pysam              bioconda/osx-64::pysam-0.15.4-py27h4ace0ce_0
      python             conda-forge/osx-64::python-2.7.15-h932b40d_1008
      python-dateutil    conda-forge/noarch::python-dateutil-2.8.1-py_0
      python_abi         conda-forge/osx-64::python_abi-2.7-1_cp27m
      pytz               conda-forge/noarch::pytz-2019.3-py_0
      qt                 anaconda/osx-64::qt-5.6.3-h1d42b2f_0
      r-base             conda-forge/osx-64::r-base-3.4.1-4
      r-lattice          conda-forge/osx-64::r-lattice-0.20_35-r341hc070d10_0
      r-locfit           conda-forge/osx-64::r-locfit-1.5_9.1-r341h470a237_2
      r-rcpp             conda-forge/osx-64::r-rcpp-1.0.0-r341h9d2a408_0
      readline           conda-forge/osx-64::readline-7.0-0
      samtools           bioconda/osx-64::samtools-0.1.19-hdd8ed8b_4
      scipy              conda-forge/osx-64::scipy-1.2.1-py27hbd7caa9_1
      setuptools         conda-forge/osx-64::setuptools-44.0.0-py27_0
      sip                conda-forge/osx-64::sip-4.18.1-py27h0a44026_1000
      six                conda-forge/noarch::six-1.14.0-py_1
      sqlite             anaconda/osx-64::sqlite-3.26.0-ha441bb4_0
      tk                 conda-forge/osx-64::tk-8.6.10-hbbe82c9_0
      trimmomatic        bioconda/noarch::trimmomatic-0.36-6
      ucsc-bedgraphtobi~ bioconda/osx-64::ucsc-bedgraphtobigwig-366-h1341992_0
      wheel              conda-forge/noarch::wheel-0.34.2-py_1
      xopen              conda-forge/osx-64::xopen-0.8.4-py27h8c360ce_1
      xz                 conda-forge/osx-64::xz-5.2.5-h0b31af3_0
      zlib               conda-forge/osx-64::zlib-1.2.11-h0b31af3_1006
    
    
    
    Downloading and Extracting Packages
    gettext-0.19.8.1     | 3.4 MB    | ##################################### | 100% 
    perl-list-moreutils- | 37 KB     | ##################################### | 100% 
    sip-4.18.1           | 243 KB    | ##################################### | 100% 
    curl-7.69.1          | 131 KB    | ##################################### | 100% 
    harfbuzz-0.9.39      | 302 KB    | ##################################### | 100% 
    libxml2-2.9.5        | 1.9 MB    | ##################################### | 100% 
    pango-1.40.3         | 448 KB    | ##################################### | 100% 
    blast-2.7.1          | 16.0 MB   | ##################################### | 100% 
    qt-5.6.3             | 67.0 MB   | ##################################### | 100% 
    boost-1.64.0         | 307 KB    | ##################################### | 100% 
    perl-archive-tar-2.1 | 27 KB     | ##################################### | 100% 
    samtools-0.1.19      | 1.8 MB    | ##################################### | 100% 
    perl-5.22.0.1        | 14.6 MB   | ##################################### | 100% 
    sqlite-3.26.0        | 2.3 MB    | ##################################### | 100% 
    freetype-2.5.5       | 738 KB    | ##################################### | 100% 
    libiconv-1.14        | 1.3 MB    | ##################################### | 100% 
    pyqt-5.6.0           | 4.4 MB    | ##################################### | 100% 
    gnutls-3.5.19        | 1.9 MB    | ##################################### | 100% 
    boost-cpp-1.64.0     | 16.2 MB   | ##################################### | 100% 
    cairo-1.14.8         | 572 KB    | ##################################### | 100% 
    matplotlib-1.5.3     | 6.0 MB    | ##################################### | 100% 
    ucsc-bedgraphtobigwi | 392 KB    | ##################################### | 100% 
    pandas-0.21.0        | 10.8 MB   | ##################################### | 100% 
    glib-2.50.2          | 2.9 MB    | ##################################### | 100% 
    fontconfig-2.12.1    | 176 KB    | ##################################### | 100% 
    libcurl-7.69.1       | 535 KB    | ##################################### | 100% 
    nettle-3.3           | 869 KB    | ##################################### | 100% 
    krb5-1.17.1          | 1.2 MB    | ##################################### | 100% 
    Preparing transaction: done
    Verifying transaction: done
    Executing transaction: done
    #
    # To activate this environment, use
    #
    #     $ conda activate interactomeseq
    #
    # To deactivate an active environment, use
    #
    #     $ conda deactivate
    


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

    Cloning into 'InteractomeSeq'...
    remote: Enumerating objects: 28, done.[K
    remote: Counting objects: 100% (28/28), done.[K
    remote: Compressing objects: 100% (21/21), done.[K
    remote: Total 575 (delta 9), reused 16 (delta 7), pack-reused 547[K
    Receiving objects: 100% (575/575), 256.33 KiB | 0 bytes/s, done.
    Resolving deltas: 100% (419/419), done.
    Checking connectivity... done.
    You have new mail in /var/mail/simone



```bash

```
