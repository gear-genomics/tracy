# Installation

Tracy is available as a [Bioconda package](https://anaconda.org/bioconda/tracy), a [Homebrew formula](https://formulae.brew.sh/formula/tracy-genomics), a [Singularity container](https://github.com/gear-genomics/tracy/releases) or as a minimal [Docker container](https://hub.docker.com/r/geargenomics/tracy/).

If desired, you can compile Tracy from source on Ubuntu or macOS to suit your specific needs.

## Installation from Source

### Ubuntu

To build Tracy from source, you need some build essentials and the Boost libraries:

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

Once you have these libraries installed, compile and link Tracy:

```bash
git clone --recursive https://github.com/gear-genomics/tracy.git
cd tracy/
make all
make install
./bin/tracy -h
```

### macOS

To build Tracy from source, you will need the following libraries:

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

For macOS, you also often need to set the `DYLD_LIBRARY_PATH` to `htslib`.

```bash
git clone --recursive https://github.com/gear-genomics/tracy.git
cd tracy/
make all
export DYLD_LIBRARY_PATH=`pwd`/src/htslib/
./src/tracy -h
```
