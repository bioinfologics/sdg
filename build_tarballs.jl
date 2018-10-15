# Note that this script can accept some limited command-line arguments, run
# `julia build_tarballs.jl --help` to see a usage message.
using BinaryBuilder

name = "BSG"
version = v"0.0.1"

# Collection of sources required to build BSGBuilder
sources = [ "./" ]

# Bash recipe for building across all platforms
script = raw"""
cd $WORKSPACE/srcdir/
cmake -DCMAKE_INSTALL_PREFIX=$prefix -DCMAKE_TOOLCHAIN_FILE=/opt/$target/$target.toolchain
make -j${nproc}
make install

"""

# These are the platforms we will build for by default, unless further
# platforms are passed in on the command line
platforms = [
    Linux(:x86_64, libc = :glibc, compiler_abi = CompilerABI(:gcc7)),
    Linux(:aarch64, libc = :glibc, compiler_abi = CompilerABI(:gcc7)),
    Linux(:powerpc64le, libc = :glibc, compiler_abi = CompilerABI(:gcc7)),
    Linux(:x86_64, libc = :musl, compiler_abi = CompilerABI(:gcc7)),
    Linux(:aarch64, libc = :musl, compiler_abi = CompilerABI(:gcc7)),
    MacOS(:x86_64, compiler_abi = CompilerABI(:gcc7)),
    FreeBSD(:x86_64, compiler_abi = CompilerABI(:gcc7))
]

# The products that we will ensure are always built
products(prefix) = [
    ExecutableProduct(prefix, "bsg-workspace", Symbol("bsg-workspace")),
    ExecutableProduct(prefix, "bsg-lhapassembler", Symbol("bsg-lhapassembler")),
    ExecutableProduct(prefix, "seq-thread", Symbol("seq-thread")),
    ExecutableProduct(prefix, "bsg-flowmaker", Symbol("bsg-flowmaker")),
    ExecutableProduct(prefix, "bsg-kmerspectra", Symbol("bsg-kmerspectra")),
    LibraryProduct(prefix, "libbsg", :libbsg),
    ExecutableProduct(prefix, "gfa-extract", Symbol("gfa-extract")),
    ExecutableProduct(prefix, "bsg-mapper", Symbol("bsg-mapper")),
    ExecutableProduct(prefix, "bsg-untangler", Symbol("bsg-untangler")),
    ExecutableProduct(prefix, "bsg-datastore", Symbol("bsg-datastore"))
]

# Dependencies that must be installed before this package can be built
dependencies = [
    "https://github.com/BenJWard/ZLibBuilder/releases/download/v1.2.11/build_ZLibBuilder.v1.2.11.jl",
    "https://github.com/staticfloat/Bzip2Builder/releases/download/v1.0.6-1/build_Bzip2.v1.0.6.jl"
]

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies)

