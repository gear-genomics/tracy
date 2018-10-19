[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/tracy/README.html)
[![Build Status](https://travis-ci.org/gear-genomics/tracy.svg?branch=master)](https://travis-ci.org/gear-genomics/tracy)
[![Docker Build](https://img.shields.io/docker/build/geargenomics/tracy.svg)](https://hub.docker.com/r/geargenomics/tracy/)
[![GitHub license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://raw.githubusercontent.com/gear-genomics/tracy/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/gear-genomics/tracy.svg)](https://github.com/gear-genomics/tracy/releases)
[![GitHub Issues](https://img.shields.io/github/issues/gear-genomics/tracy.svg)](https://github.com/gear-genomics/tracy/issues)


Installing Tracy
----------------

The easiest way to get Tracy is to download a statically linked binary from the [Tracy release page](https://github.com/gear-genomics/tracy/releases) or to download Tracy from [Bioconda](https://anaconda.org/bioconda/tracy). Building from source is also possible:

`apt-get install -y build-essential g++ cmake zlib1g-dev libbz2-dev liblzma-dev libboost-all-dev`

`git clone --recursive https://github.com/gear-genomics/tracy.git`

`cd tracy/`

`make all`

`make install`

Running Tracy
-------------

`./tracy -h`


Basecalling a Trace File
------------------------

To get the primary sequence (highest peak) of a trace file in FASTA or FASTQ format.

`./tracy basecall -f fasta -o out.fasta input.ab1`

`./tracy basecall -f fastq -o out.fastq input.ab1`

To get full trace information, including primary and secondary basecalls for heterozygous variants.

`./tracy basecall -f tsv -o out.tsv input.ab1`


Alignment to a Fasta Slice
--------------------------

Alignment of a trace file to a FASTA reference slice.

`./tracy align -f align -o out.align -g ref_slice.fa input.ab1`


Alignment to an indexed reference genome
----------------------------------------

Alignment to a large reference genome requires a pre-built index on a bgzip compressed genome.

`./tracy index -o hg19.fa.fm9 hg19.fa.gz`

`samtools faidx hg19.fa.gz`

Once the index has been built you can align to the indexed genome.

`./tracy align -g hg19.fa.gz input.ab1`


Separating heterozygous variants
--------------------------------

Double-peaks in the Chromatogram can cause alignment issues. Tracy supports deconvolution of heterozygous variants into two separate alleles.

`./tracy decompose -g hg19.fa.gz -f align -o outprefix input.ab1`

The two alleles are then separately aligned.

`cat outprefix.align1 outprefix.align2`

You can also use a wildtype chromatogram for decomposition.

`./tracy decompose -g wildtype.ab1 -f align -o outprefix mutated.ab1`


Questions
---------

In case of questions feel free to send us an [email](https://www-db.embl.de/EMBLPersonGroup-PersonPicture/MailForm/?recipient=ggenomics).
