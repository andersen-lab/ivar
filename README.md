iVar
===========

![C/C++ CI](https://github.com/andersen-lab/ivar/workflows/C/C++%20CI/badge.svg) [![DOI](https://zenodo.org/badge/143471288.svg)](https://zenodo.org/badge/latestdoi/143471288) [![install with conda](https://anaconda.org/bioconda/ivar/badges/version.svg)](https://anaconda.org/bioconda/ivar) [![install with conda](https://anaconda.org/bioconda/ivar/badges/platforms.svg)](https://anaconda.org/bioconda/ivar)


iVar is a computational package that contains functions broadly useful for viral amplicon-based sequencing. Additional tools for metagenomic sequencing are actively being incorporated into iVar. While each of these functions can be accomplished using existing tools, iVar contains an intersection of functionality from multiple tools that are required to call iSNVs and consensus sequences from viral sequencing data across multiple replicates. We implemented the following functions in iVar: (1) trimming of primers and low-quality bases, (2) consensus calling, (3) variant calling - both iSNVs and insertions/deletions, and (4) identifying mismatches to primer sequences and excluding the corresponding reads from alignment files.

[An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar](https://doi.org/10.1186/s13059-018-1618-7)

*Genome Biology* 2019 **20**:8

Nathan D Grubaugh, Karthik Gangavarapu, Joshua Quick, Nathaniel L Matteson, Jaqueline Goes De Jesus, Bradley J Main, Amanda L Tan, Lauren M Paul, Doug E Brackney, Saran Grewal, Nikos Gurfield, Koen KA Van Rompay, Sharon Isern, Scott F Michael, Lark L Coffey, Nicholas J Loman, Kristian G Andersen

bioRxiv doi: [https://doi.org/10.1101/383513](https://doi.org/10.1101/383513)

## Manual

Manual for iVar is available [here](https://andersen-lab.github.io/ivar/html/).

## Installation

### Dependencies

* [HTSlib](http://www.htslib.org/download/)
* [GCC](https://gcc.gnu.org/) any version after v5.0. Support for C++11 standard required.

Note:
* It is highly recommended that [samtools](https://github.com/samtools/samtools) also be installed alongside iVar. iVar uses the output of samtools mpileup to call variants and generate consensus sequences. In addition, samtools `sort` and `index` commands are very useful to setup a pipeline using iVar.

Installing via conda
====================

iVar is available on bioconda. To install conda, please use the [miniconda](https://conda.io/miniconda.html) package. After intalling conda please add the following channels,

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

To install iVar,

```
conda install ivar
```


Installing via homebrew
=======================

iVar can be installed using [Homebrew](https://brew.sh/).

```
brew install brewsci/bio/ivar
```


Installing on Mac
=================

Installing build tools
----------------------

[Xcode](https://developer.apple.com/xcode/) from Apple is required to compile iVar (and other tools) from source. If you don't want to install the full Xcode package from the AppStore, you can install the Xcode command line tools,

```
xcode-select --install
```

[GNU Autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html#Autotools-Introduction) is required to compile iVar from source.

To install Autotools using [homebrew](https://brew.sh/) please use the command below,

```
brew install autoconf automake libtool
```

HTSlib installed using conda
-----------------------------

HTSlib can be installed with [conda](https://conda.io/docs/) using the command,

```
conda install -c bioconda htslib
```

The conda binary is by default installed at /opt/. You can check the installation location by running the following command,

```
which conda
```

The output of the command will be in this format - /opt/conda/bin/conda or /opt/anaconda2/bin/conda or /opt/anaconda3/bin/conda depending on whether you installed miniconda or anaconda.

If the output is for example, /opt/conda/bin/conda, then you can add the path to the lib folder to $LD_LIBRARY_PATH using the command below.
You can add this to your ~/.bash_profile or ~/.bashrc to avoid rerunning the command everytime a new bash session starts.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/lib
```

HTSlib installed from source
----------------------------

Installation instructions and downloads for HTSlib can be found at http://www.htslib.org/download/.

If HTSlib is installed in a non standard location, please add the following to your .bash_profile so that iVar can find HTSlib dynamic libraries during runtime.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/hts/lib/folder
```

Installing iVar
---------------

To install iVar, run the following commands.

```
./autogen.sh
./configure
make
make install
```

If HTSlib was installed using conda, please run the following commands by supplying the prefix to the bin folder of the conda binary.

The prefix to the bin folder can be found using the command `which conda`. The output of the command will be in this format - /opt/conda/bin/conda or /opt/anaconda2/bin/conda or /opt/anaconda3/bin/conda depending on whether you installed miniconda or anaconda. For example, if the output of the command is /opt/conda/bin/conda, the prefix to the htslib bin folder will be /opt/conda. This can be supplied to ./configure --with-hts=/opt/conda.

```
./autogen.sh
./configure --with-hts=/prefix/to/bin/folder/with/HTSlib
make
make install
```

If HTSlib was installed in a non standard location, please run the following commands,

```
./autogen.sh
./configure --with-hts=/prefix/to/bin/folder/with/HTSlib
make
make install
```

To test installation just run, `ivar version` and you should get the following output,

```
iVar version 1.0

Please raise issues and bug reports at https://github.com/andersen-lab/ivar/
```

Installing on Linux
===================

Installing build tools
----------------------

[GNU Autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html#Autotools-Introduction) is required to compile iVar from source.

To install Autotools using [APT](https://help.ubuntu.com/lts/serverguide/apt.html) please use the command below,

```
apt-get install autotools-dev
```

HTSlib installed using conda
-----------------------------

HTSlib can be installed with [conda](https://conda.io/docs/) using the command,

```
conda install -c bioconda htslib
```

The conda binary is by default installed at /opt/. You can check the installation location by running the following command,

```
which conda
```

The output of the command will be in this format - /opt/conda/bin/conda or /opt/anaconda2/bin/conda or /opt/anaconda3/bin/conda depending on whether you installed miniconda or anaconda.

If the output is for example, /opt/conda/bin/conda, then you can add the path to the lib folder to $LD_LIBRARY_PATH using the command below.
You can add this to your ~/.bash_profile or ~/.bashrc to avoid rerunning the command everytime a new bash session starts.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/lib
```

HTSlib installed from source
----------------------------

Installation instructions and downloads for HTSlib can be found at http://www.htslib.org/download/.

If HTSlib is installed in a non standard location, please add the following to your .bash_profile so that iVar can find HTSlib dynamic libraries during runtime.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/hts/lib/folder
```

Installing iVar
---------------

To install iVar, run the following commands.

```
./autogen.sh
./configure
make
make install
```

If HTSlib was installed using conda, please run the following commands by supplying the prefix to the bin folder of the conda binary.

The prefix to the bin folder can be found using the command `which conda`. The output of the command will be in this format - /opt/conda/bin/conda or /opt/anaconda2/bin/conda or /opt/anaconda3/bin/conda depending on whether you installed miniconda or anaconda. For example, if the output of the command is /opt/conda/bin/conda, the prefix to the htslib bin folder will be /opt/conda. This can be supplied to ./configure --with-hts=/opt/conda.

```
./autogen.sh
./configure --with-hts=/prefix/to/bin/folder/with/HTSlib
make
make install
```

If HTSlib was installed in a non standard location, please run the following commands,

```
./autogen.sh
./configure --with-hts=/prefix/to/bin/folder/with/HTSlib
make
make install
```

To test installation just run, `ivar version` and you should get the following output,

```
iVar version 1.0

Please raise issues and bug reports at https://github.com/andersen-lab/ivar/
```


Running from Docker
===================

iVar can also be run via [Docker](https://www.docker.com/). Pull the docker image from [Docker Hub](https://hub.docker.com/) using the following command,

```
docker pull andersenlabapps/ivar
```

This docker image contains all the required dependencies to run iVar and the [pipelines](@ref cookbookpage) developed using iVar.
You will have to attach a docker volume to get data into the docker container. Instructions to do so are in the [Docker docs](https://docs.docker.com/storage/volumes/).

[iVar on Docker Hub](https://hub.docker.com/r/andersenlabapps/ivar/)

Contact
=======

For bug reports please email gkarthik[at]scripps.edu or raise an issue on Github.

Acknowledgements
=======

This work was supported in part by NIH grants U19AI135995, R21AI137690, and UL1TR002550.

