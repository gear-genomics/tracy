[![Build Status](https://travis-ci.org/gear-genomics/tracy.svg?branch=master)](https://travis-ci.org/gear-genomics/tracy)
[![Docker Build](https://img.shields.io/docker/build/geargenomics/tracy.svg)](https://hub.docker.com/r/geargenomics/tracy/)
[![GitHub license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://raw.githubusercontent.com/gear-genomics/tracy/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/gear-genomics/tracy.svg)](https://github.com/gear-genomics/tracy/releases)
[![GitHub Issues](https://img.shields.io/github/issues/gear-genomics/tracy.svg)](https://github.com/gear-genomics/tracy/issues)


Installing Tracy
----------------

The easiest way to get Tracy is to download a statically linked binary from the [Tracy release page](https://github.com/gear-genomics/tracy/releases). Building from source is also possible:

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
