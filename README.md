[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/tracy/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/tracy/badges/downloads.svg)](https://anaconda.org/bioconda/tracy)
[![Build Status](https://travis-ci.org/gear-genomics/tracy.svg?branch=master)](https://travis-ci.org/gear-genomics/tracy)
[![Docker Build](https://img.shields.io/docker/build/geargenomics/tracy.svg)](https://hub.docker.com/r/geargenomics/tracy/)
[![GitHub license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://raw.githubusercontent.com/gear-genomics/tracy/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/gear-genomics/tracy.svg)](https://github.com/gear-genomics/tracy/releases)
[![GitHub Issues](https://img.shields.io/github/issues/gear-genomics/tracy.svg)](https://github.com/gear-genomics/tracy/issues)


Installing Tracy
----------------

The easiest way to get Tracy is to download the statically linked binary or the singularity container (SIF file) from the [Tracy release page](https://github.com/gear-genomics/tracy/releases). Alternatively, you can download Tracy from [Bioconda](https://anaconda.org/bioconda/tracy) or pull the [Tracy docker container](https://hub.docker.com/r/geargenomics/tracy/).

Building from Source
--------------------

`git clone --recursive https://github.com/gear-genomics/tracy.git`

`cd tracy/`

`make all`

`make install`

Tracy requires some system libraries such as bzip2, zlib and boost. For Ubuntu Linux you install these using:

`apt-get install -y build-essential g++ cmake zlib1g-dev libbz2-dev liblzma-dev libboost-all-dev`

The Mac OSX versions of these packages are:

`brew install cmake zlib readline xz bzip2 gsl libtool pkg-config boost`

For Mac OSX you also often need to set the library path to HTSlib.

`cd tracy/`

`export DYLD_LIBRARY_PATH=`pwd`/src/htslib/`


Running Tracy
-------------

`tracy -h`


Basecalling a Trace File
------------------------

To get the primary sequence (highest peak) of a trace file in FASTA or FASTQ format.

`tracy basecall -f fasta -o out.fasta input.ab1`

`tracy basecall -f fastq -o out.fastq input.ab1`

To get full trace information, including primary and secondary basecalls for heterozygous variants.

`tracy basecall -f tsv -o out.tsv input.ab1`


Alignment to a Fasta Slice
--------------------------

Alignment of a trace file to a FASTA reference slice.

`tracy align -f align -o out.align -g ref_slice.fa input.ab1`


Alignment to an indexed reference genome
----------------------------------------

Alignment to a large reference genome requires a pre-built index on a bgzip compressed genome.

`tracy index -o hg38.fa.fm9 hg38.fa.gz`

`samtools faidx hg38.fa.gz`

Once the index has been built you can align to the indexed genome.

`tracy align -g hg38.fa.gz input.ab1`


Separating heterozygous variants
--------------------------------

Double-peaks in the Chromatogram can cause alignment issues. Tracy supports deconvolution of heterozygous variants into two separate alleles.

`tracy decompose -g hg38.fa.gz -f align -o outprefix input.ab1`

The two alleles are then separately aligned.

`cat outprefix.align1 outprefix.align2`

You can also use a wildtype chromatogram for decomposition.

`tracy decompose -g wildtype.ab1 -f align -o outprefix mutated.ab1`


SNV & InDel Variant Calling and Annotation
--------------------------------------------

Tracy can call and annotate variants with respect to a reference genome.

`tracy decompose -v -a homo_sapiens -g hg38.fa.gz -f align -o outprefix input.ab1`

This command produces a variant call file in binary BCF format. It can be converted to VCF using bcftools.

`bcftools view outprefix.bcf`


Using forward & reverse ab1 files to improve variant calling
------------------------------------------------------------

If you do have forward and reverse trace files for the same expected genomic variant you can merge variant files and check consistency of calls and genotypes.

Forward trace decomposition:

`tracy decompose -f align -o forward -a homo_sapiens -g hg38.fa.gz forward.ab1`

Reverse trace decomposition:

`tracy decompose -f align -o reverse -a homo_sapiens -g hg38.fa.gz reverse.ab1`

Left-alignment of InDels:

`bcftools norm -O b -o forward.norm.bcf -f hg38.fa.gz forward.bcf`

`bcftools norm -O b -o reverse.norm.bcf -f hg38.fa.gz reverse.bcf`

Merging of normalized variant files:

`bcftools merge --force-samples forward.norm.bcf reverse.norm.bcf`


Trace assembly
--------------

If you tiled a genomic region with multiple chromatogram files you can assemble all of these with tracy. 

`tracy assemble -r <reference.fa> file1.ab1 file2.ab1 fileN.ab1`

Tracy also supports de novo assembly if chromatogram trace files overlap sufficiently with each other.

`tracy assemble file1.ab1 file2.ab1 fileN.ab1`

Questions
---------

In case of questions feel free to send us an [email](https://www-db.embl.de/EMBLPersonGroup-PersonPicture/MailForm/?recipient=ggenomics).
