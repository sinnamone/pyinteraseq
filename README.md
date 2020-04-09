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

Test if miniconda is installed


```bash
which conda
```

    /home/spuccio/miniconda3/condabin/conda



```bash
conda create -n interactomeseq python=2.7.14 matplotlib=1.5.3 trimmomatic=0.36 kallisto=0.46.0 cutadapt=1.12 blast r-base=3.4.1 pandas=0.21.0 samtools=0.1.19  cython biom-format biopython pybedtools bioconductor-edger ucsc-bedgraphtobigwig -c conda-forge -c bioconda -c r -c anaconda -y
```

    Collecting package metadata (current_repodata.json): done
    Solving environment: failed with repodata from current_repodata.json, will retry with next repodata source.
    Collecting package metadata (repodata.json): done
    Solving environment: done
    
    ## Package Plan ##
    
      environment location: /home/spuccio/miniconda3/envs/interactomeseq
    
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
        - python=2.7.14
        - r-base=3.4.1
        - samtools=0.1.19
        - trimmomatic=0.36
        - ucsc-bedgraphtobigwig
    
    
    The following packages will be downloaded:
    
        package                    |            build
        ---------------------------|-----------------
        python-2.7.14              |                5        13.8 MB  conda-forge
        ------------------------------------------------------------
                                               Total:        13.8 MB
    
    The following NEW packages will be INSTALLED:
    
      _libgcc_mutex      conda-forge/linux-64::_libgcc_mutex-0.1-conda_forge
      _openmp_mutex      conda-forge/linux-64::_openmp_mutex-4.5-1_llvm
      _r-mutex           conda-forge/noarch::_r-mutex-1.0.1-anacondar_1
      bedtools           bioconda/linux-64::bedtools-2.29.2-hc088bd4_0
      bioconductor-edger bioconda/linux-64::bioconductor-edger-3.22.5-r341hfc679d8_0
      bioconductor-limma bioconda/linux-64::bioconductor-limma-3.36.5-r341h470a237_0
      biom-format        conda-forge/linux-64::biom-format-2.1.7-py27h3010b51_1002
      biopython          conda-forge/linux-64::biopython-1.76-py27h516909a_0
      blast              bioconda/linux-64::blast-2.6.0-boost1.61_1
      boost              anaconda/linux-64::boost-1.61.0-py27_0
      bz2file            conda-forge/noarch::bz2file-0.98-py_0
      bzip2              conda-forge/linux-64::bzip2-1.0.8-h516909a_2
      ca-certificates    conda-forge/linux-64::ca-certificates-2020.4.5.1-hecc5488_0
      cairo              anaconda/linux-64::cairo-1.14.8-0
      certifi            conda-forge/linux-64::certifi-2019.11.28-py27h8c360ce_1
      click              conda-forge/noarch::click-7.1.1-pyh8c360ce_0
      curl               conda-forge/linux-64::curl-7.64.0-h646f8bb_0
      cutadapt           bioconda/linux-64::cutadapt-1.12-py27_1
      cycler             conda-forge/noarch::cycler-0.10.0-py_2
      cython             conda-forge/linux-64::cython-0.29.15-py27haf22ab1_1
      dbus               conda-forge/linux-64::dbus-1.13.0-h4e0c4b3_1000
      expat              conda-forge/linux-64::expat-2.2.9-he1b5a44_2
      fontconfig         anaconda/linux-64::fontconfig-2.12.1-3
      freetype           anaconda/linux-64::freetype-2.5.5-2
      future             conda-forge/linux-64::future-0.18.2-py27h8c360ce_1
      glib               anaconda/linux-64::glib-2.50.2-1
      gsl                conda-forge/linux-64::gsl-2.6-h294904e_0
      gst-plugins-base   conda-forge/linux-64::gst-plugins-base-1.8.0-0
      gstreamer          conda-forge/linux-64::gstreamer-1.8.0-1
      harfbuzz           anaconda/linux-64::harfbuzz-0.9.39-2
      hdf5               conda-forge/linux-64::hdf5-1.10.5-nompi_h3c11f04_1104
      icu                anaconda/linux-64::icu-54.1-0
      jpeg               conda-forge/linux-64::jpeg-9c-h14c3975_1001
      kallisto           bioconda/linux-64::kallisto-0.46.0-h4f7b962_1
      krb5               conda-forge/linux-64::krb5-1.16.3-hc83ff2d_1000
      libblas            conda-forge/linux-64::libblas-3.8.0-16_openblas
      libcblas           conda-forge/linux-64::libcblas-3.8.0-16_openblas
      libcurl            conda-forge/linux-64::libcurl-7.64.0-h01ee5af_0
      libdeflate         conda-forge/linux-64::libdeflate-1.5-h516909a_0
      libedit            conda-forge/linux-64::libedit-3.1.20170329-0
      libffi             conda-forge/linux-64::libffi-3.2.1-he1b5a44_1007
      libgcc             conda-forge/linux-64::libgcc-7.2.0-h69d50b8_2
      libgcc-ng          conda-forge/linux-64::libgcc-ng-9.2.0-h24d8f2e_2
      libgfortran-ng     conda-forge/linux-64::libgfortran-ng-7.3.0-hdf63c60_5
      libiconv           conda-forge/linux-64::libiconv-1.14-4
      liblapack          conda-forge/linux-64::liblapack-3.8.0-16_openblas
      libopenblas        conda-forge/linux-64::libopenblas-0.3.9-h5ec1e0e_0
      libpng             conda-forge/linux-64::libpng-1.6.37-hed695b0_1
      libssh2            conda-forge/linux-64::libssh2-1.8.0-h1ad7b7a_1003
      libstdcxx-ng       conda-forge/linux-64::libstdcxx-ng-9.2.0-hdf63c60_2
      libtiff            conda-forge/linux-64::libtiff-4.1.0-hc7e4089_6
      libuuid            conda-forge/linux-64::libuuid-2.32.1-h14c3975_1000
      libwebp-base       conda-forge/linux-64::libwebp-base-1.1.0-h516909a_3
      libxcb             conda-forge/linux-64::libxcb-1.13-h14c3975_1002
      libxml2            anaconda/linux-64::libxml2-2.9.9-hea5a465_1
      llvm-openmp        conda-forge/linux-64::llvm-openmp-9.0.1-hc9558a2_2
      lz4-c              conda-forge/linux-64::lz4-c-1.8.3-he1b5a44_1001
      matplotlib         anaconda/linux-64::matplotlib-1.5.3-np111py27_1
      mysql-connector-c  conda-forge/linux-64::mysql-connector-c-6.1.11-hab6429c_1002
      ncurses            conda-forge/linux-64::ncurses-5.9-10
      numpy              conda-forge/linux-64::numpy-1.11.3-py27he5ce36f_1207
      openjdk            conda-forge/linux-64::openjdk-11.0.1-h516909a_1016
      openssl            conda-forge/linux-64::openssl-1.0.2u-h516909a_0
      pandas             conda-forge/linux-64::pandas-0.21.0-py27_0
      pango              anaconda/linux-64::pango-1.40.3-1
      pcre               conda-forge/linux-64::pcre-8.39-0
      pip                conda-forge/noarch::pip-20.0.2-py_2
      pixman             conda-forge/linux-64::pixman-0.34.0-h14c3975_1003
      pthread-stubs      conda-forge/linux-64::pthread-stubs-0.4-h14c3975_1001
      pybedtools         bioconda/linux-64::pybedtools-0.8.1-py27he513fc3_0
      pycairo            conda-forge/linux-64::pycairo-1.16.3-py27_0
      pyparsing          conda-forge/noarch::pyparsing-2.4.7-pyh9f0ad1d_0
      pyqi               conda-forge/noarch::pyqi-0.3.2-py_0
      pyqt               conda-forge/linux-64::pyqt-5.6.0-py27h13b7fb3_1008
      pysam              bioconda/linux-64::pysam-0.15.4-py27hbcae180_0
      python             conda-forge/linux-64::python-2.7.14-5
      python-dateutil    conda-forge/noarch::python-dateutil-2.8.1-py_0
      python_abi         conda-forge/linux-64::python_abi-2.7-1_cp27mu
      pytz               conda-forge/noarch::pytz-2019.3-py_0
      qt                 anaconda/linux-64::qt-5.6.2-5
      r-base             r/linux-64::r-base-3.4.1-1
      r-lattice          conda-forge/linux-64::r-lattice-0.20_35-r341hc070d10_0
      r-locfit           conda-forge/linux-64::r-locfit-1.5_9.1-r341h470a237_2
      r-rcpp             conda-forge/linux-64::r-rcpp-1.0.0-r341h9d2a408_0
      readline           conda-forge/linux-64::readline-7.0-0
      samtools           bioconda/linux-64::samtools-0.1.19-3
      scipy              conda-forge/linux-64::scipy-1.2.1-py27h09a28d5_1
      setuptools         conda-forge/linux-64::setuptools-44.0.0-py27_0
      sip                conda-forge/linux-64::sip-4.18.1-py27hf484d3e_1000
      six                conda-forge/noarch::six-1.14.0-py_1
      sqlite             conda-forge/linux-64::sqlite-3.20.1-2
      tk                 conda-forge/linux-64::tk-8.6.10-hed695b0_0
      trimmomatic        bioconda/noarch::trimmomatic-0.36-6
      ucsc-bedgraphtobi~ bioconda/linux-64::ucsc-bedgraphtobigwig-377-h199ee4e_0
      wheel              conda-forge/noarch::wheel-0.34.2-py_1
      xopen              conda-forge/linux-64::xopen-0.8.4-py27h8c360ce_1
      xorg-libxau        conda-forge/linux-64::xorg-libxau-1.0.9-h14c3975_0
      xorg-libxdmcp      conda-forge/linux-64::xorg-libxdmcp-1.1.3-h516909a_0
      xz                 conda-forge/linux-64::xz-5.2.5-h516909a_0
      zlib               conda-forge/linux-64::zlib-1.2.11-h516909a_1006
      zstd               conda-forge/linux-64::zstd-1.4.4-h3b9ef0a_2
    
    
    
    Downloading and Extracting Packages
    python-2.7.14        | 13.8 MB   | ##################################### | 100% 
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

    (interactomeseq) 




```bash
conda list | less 
```

If you decide later that you don’t want the environment or its packages anymore, deactivate the environment and then run this command:


```bash
# Commented to avoid deletion
#conda remove --name interactomeseq --all
```

### Installing additional dependencies accessible in your interactomeseq environment

Then install qiime1 version 1.91.1

This should be launched as control 


```bash
conda install -c bioconda pysam -y
```

    Collecting package metadata (current_repodata.json): done
    Solving environment: done
    
    # All requested packages already installed.
    
    (interactomeseq) 




```bash
conda install -c conda-forge readline=6.2 -y
```

    Collecting package metadata (current_repodata.json): done
    Solving environment: failed with initial frozen solve. Retrying with flexible solve.
    Solving environment: failed with repodata from current_repodata.json, will retry with next repodata source.
    Collecting package metadata (repodata.json): done
    Solving environment: done
    
    ## Package Plan ##
    
      environment location: /home/spuccio/miniconda3/envs/interactomeseq
    
      added / updated specs:
        - readline=6.2
    
    
    The following packages will be REMOVED:
    
      krb5-1.16.3-hc83ff2d_1000
      libcurl-7.64.0-h01ee5af_0
    
    The following packages will be DOWNGRADED:
    
      curl                                    7.64.0-h646f8bb_0 --> 7.52.1-0
      python                                           2.7.14-5 --> 2.7.14-2
      readline                                            7.0-0 --> 6.2-0
      sqlite                                           3.20.1-2 --> 3.13.0-1
      tk                                      8.6.10-hed695b0_0 --> 8.5.19-2
    
    
    Preparing transaction: done
    Verifying transaction: done
    Executing transaction: done
    (interactomeseq) 




```bash
 pip install qiime 
```

### Clone the Github repository


```bash
git clone https://github.com/sinnamone/InteractomeSeq
```
