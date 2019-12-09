# Installation

Tracy is available as a [Bioconda package](https://anaconda.org/bioconda/tracy), as a pre-compiled statically linked binary from [Tracy's github release page](https://github.com/gear-genomics/tracy/releases), as a singularity container [SIF file](https://github.com/gear-genomics/tracy/releases) or as a minimal [Docker container](https://hub.docker.com/r/geargenomics/tracy/).


## Installation from Source

To build Tracy from source you need some build essentials and the Boost libraries, i.e. for Ubuntu:

```bash
apt install \
    build-essential g++ \
    cmake \
    git-all \
    liblzma-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev
```

Once you have installed these system libraries you can compile and link Tracy.

```bash
git clone --recursive https://github.com/gear-genomics/tracy.git
cd tracy/
make all
make install
./bin/tracy -h
```

## Installation for Mac OSX

To build Tracy from source you need some system libraries.

```bash
brew install \
     cmake \
     zlib \
     readline \
     xz \
     bzip2 \
     gsl \
     libtool \
     pkg-config \
     boost
```

For Mac OSX you also often need to set the library path to HTSlib.

```bash
git clone --recursive https://github.com/gear-genomics/tracy.git
cd tracy/
make all
export DYLD_LIBRARY_PATH=`pwd`/src/htslib/
./src/tracy -h
```
