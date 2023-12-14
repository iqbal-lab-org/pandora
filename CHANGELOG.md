# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.12.0-alpha.0]

### Fixed

- `pandora` mapping has been improved by doing a better detection of conflicting clusters and selection [[#344][344]];

### Added

- Parameter `--min-gene-coverage-proportion` to `pandora` `map`, `compare` and `discover` subcommands [[#351][351]];
- Parameter `--no-gene-coverage-filtering` to `pandora` `map`, `compare` and `discover` subcommands [[#352][352]];
- Parameter `--partial-matching-lower-bound` to `pandora` `map`, `compare` and `discover` subcommands [[#353][353]];

## [0.11.0-alpha.0]

This version is a major release that breaks backwards compatibility with previous versions of `pandora`.
It improves `pandora` runtime performance by 15x and RAM usage by 20x;

### Changed
- The `pandora` index changed from a set of files in a directory structure to a single, compressible and indexable `zip`
file (`pandora` indexes now have the suffix `.panidx.zip`). This is now the single file that is produced by the
`pandora index` command and is required as argument to all the other `pandora` commands. This index is self contained in 
the sense that it encodes all the information and metadata about it (e.g. which PRGs were used to create it, window and 
kmer size, etc). This new index provide the infrastructure for the next features and simplifies working with large 
reference pangenome collections, with a few million PRGs. This new index breaks backwards compatibility with previous 
`pandora` versions. The structure of this zip archive is as follows:
  * `_prg_names`: The names of the PRGs used as input to create this index;
  * `_prg_max_path_lengths`: the length of the longest path through each PRG;
  * `_prg_lengths`: the length of the string representation of each PRG;
  * `_minhash`: the minimizer hash data structure;
  * `_metadata`: metadata about the index;
  * `*.gfa`: the several GFA files describing the minimizing kmer graph for each PRG;
  * `*.fa`: the string representation of each PRG;
- Minimum C++ standard upgraded from `C++11` to `C++14`;
- We now test whether the genotype confidence of a variant is greater than or equal to the threshold provided by
`--gt-conf`. Previously we only tested if it was greater than;

### Removed
- Removed CLI parameters `-w`, `-k` and `--clean` from the following `pandora` subcommands: `compare`, `discover`, `map`,
`seq2path`;   
- Removed `merge_index` subcommand;
- Removed gene-DBG and noise-filtering modules;

### Fixed
- Fixed a major bug on finding the longest path through PRGs;
- Several refactorings to the `pandora` index implementation;
- Optimisation of the `pandora` index data structure;
 
### Added
- A memory-efficient way to load PRGs when indexing and mapping, where we don't need to load all PRGs at once to process 
them, but just load on demand (also known as lazy loading). This is particularly useful when working with very large 
PanRGs;
- Random multimapping of reads if they map equally well to several graphs, reducing mapping bias. Added parameter
`--rng-seed` to `pandora map/compare/discover` commands to make multimapping deterministic, if required;
- A new parameter to deal with auto-updating error rate and kmer model (see `--dont-auto-update-params` parameter in 
`pandora map/compare/discover` commands);
- Three new parameters to control when a gene should be filtered out due to too low or too high coverage (see
`--min-abs-gene-coverage`, `--min-rel-gene-coverage` and `--max-rel-gene-coverage` parameters in
`pandora map/compare/discover` commands);


## [0.10.0-alpha.0]

### Changed

- Denovo discovery is now done by repeatedly polishing the loci's maximum likelihood sequences using the regions of the
reads that mapped to the loci through [Racon](Racon);
- Pandora `discover` CLI heavily changed: parameters `-M,--mapped-reads`, `--clean-dbg`, `--discover-k`, `--max-ins`,
`--covg-threshold`, `-l`, `-L`, `-d,--merge`, `-N`, `--min-dbg-dp` removed;

### Added
- Pandora `map`, `compare` and `discover` commands now produce [SAM](SAM) files;
- Parameter `-K`/`--debugging-files` to pandora `map`, `compare` and `discover` commands to create extra
debugging files, which are able to describe completely the mapping process of `pandora`.


## [0.9.2]

### Changed

- The VCF INFO field `SVTYPE` has now been changed to `VC` [[#249][249]]

### Fixed

- More robust TSV file parsing. Empty line no longer required at end [[#213][213]]
- Handle ambiguous bases properly instead of skipping to next read once we reach one [[#294][294]]

## [0.9.1]

### Added

- `pandora` is now installable through `conda`;
- A script to archive the `pandora` repository with git submodules;

### Changed

- Improved the sample example so now we can assert that the output produced is the
  expected one;
- Changes to the build process that enables `pandora` to be compiled in the `conda`
  environment;

## [0.9.0]

### Changed

- Version bump from `0.9.0-rc2` to `0.9.0`.

## [0.9.0-rc2]

### Changed

- `pandora discover` now processes one sample at a time, but runs with several threads
  on the heavy tasks, i.e. when mapping reads, finding candidate regions, and finding
  denovo variants. The result is that it now takes a lot less RAM to run on multiple
  samples.

## [0.9.0-rc1]

### Changed

- `pandora discover` now receives read index files describing samples and reads, and
  discover denovo sequences in these samples. To improve performance on discovering
  denovo sequences on several samples, `pandora discover` is now multithreaded, but the
  performance is still the same as the previous version, i.e. each sample is processed
  in a single-threaded way;
- `pandora discover` output changed to a proprietary format. See [example](example) for
  the new output;
- `pandora` can now communicate with a
  [`make_prg` prototype](https://github.com/leoisl/make_prg) that is able to update PRGs
  without needing to realign and remake the PRG. This provides major performance
  upgrades to running the full `pandora` pipeline with denovo discovery enabled, and
  there is no need anymore to use a `snakemake` pipeline (see
  [this example](example/run_pandora.sh) to how to run the full pipeline);
- We now use [musl libc](https://musl.libc.org/) instead of
  [Holy Build Box](https://github.com/phusion/holy-build-box) to build a precompiled
  portable binary, removing the dependency on `OpenMP 4.0+` or `GCC 4.9+`, and `GLIBC`;

## [0.8.0]

### Added

- We now provide a script to build a portable precompiled binary as another option to
  run `pandora` easily. The portable binary is now provided with the release;
- `pandora` can now provide a meaningful stack trace in case of errors, to facilitate
  debugging (need to pass flag `-DPRINT_STACKTRACE` to `CMake`). Due to this, we now add
  debug symbols (`-g` flag) to every `pandora` build type, but this
  [does not impact performance](https://stackoverflow.com/a/39223245). The precompiled
  binary has this enabled.

### Changed

- We now use the [Hunter](https://github.com/cpp-pm/hunter) package manager, removing
  the requirement of having `ZLIB` and `Boost` system-wide installations;
- `GATB` is now a git submodule instead of an external project downloaded and compiled
  during compilation time. This means that when git cloning `pandora`, `cgranges` and
  `GATB` are also downloaded/cloned, and when preparing the build (running `cmake`),
  `Hunter` downloads and installs `Boost`, `GTest` and `ZLIB`. Thus we still need
  internet connection to prepare the build (running `cmake`) but not for compiling
  (running `make`).
- We now use a GATB fork that accepts a `ZLIB` custom installation;
- Refactored all thirdparty libraries (`cgranges`, `GATB`, `backward`, `CLI11`,
  `inthash`) into their own directory `thirdparty`.

### Fixed

- Refactored asserts into exceptions, and now `pandora` can be compiled correctly in the
  `Release` mode. The build process is thus able to create a more optimized binary,
  resulting in improved performance.
- Don't assume Nanopore reads are longer than loci [[#265][265]]

## [v0.7.0]

There is a significant amount of changes to the project between version 0.6 and this
release. Only major things are listed here. Future releases from this point will have
their changes meticulously documented here.

### Added

- `discover` subcommand for de novo variant discovery [[#234][234]]
- many more tests

### Changed

- FASTA/Q files are now parsed with `klib` [[#223][223]]
- command-line interface is now overhauled with many breaking changes [[#224][224]]
- global genotyping has been made default [[#220][220]]
- Various improvements to VCF-related functions

### Fixed

- k-mer coverage underflow bug in `LocalPRG` [[#183][183]]

[Unreleased]: https://github.com/rmcolq/pandora/compare/0.12.0-alpha.0...HEAD
[0.12.0-alpha.0]: https://github.com/rmcolq/pandora/compare/0.12.0-alpha.0...0.11.0-alpha.0
[0.11.0-alpha.0]: https://github.com/rmcolq/pandora/compare/0.11.0-alpha.0...0.10.0-alpha.0
[0.10.0-alpha.0]: https://github.com/rmcolq/pandora/compare/0.10.0-alpha.0...0.9.2
[0.9.2]: https://github.com/rmcolq/pandora/compare/0.9.2...0.9.1
[0.9.1]: https://github.com/rmcolq/pandora/releases/tag/0.9.1
[0.9.0]: https://github.com/rmcolq/pandora/releases/tag/0.9.0
[0.9.0-rc2]: https://github.com/rmcolq/pandora/releases/tag/0.9.0-rc2
[0.9.0-rc1]: https://github.com/rmcolq/pandora/releases/tag/0.9.0-rc1
[0.8.0]: https://github.com/rmcolq/pandora/releases/tag/0.8.0
[183]: https://github.com/rmcolq/pandora/issues/183
[213]: https://github.com/rmcolq/pandora/issues/213
[220]: https://github.com/rmcolq/pandora/pull/220
[223]: https://github.com/rmcolq/pandora/pull/223
[224]: https://github.com/rmcolq/pandora/pull/224
[234]: https://github.com/rmcolq/pandora/pull/234
[249]: https://github.com/rmcolq/pandora/issues/249
[265]: https://github.com/rmcolq/pandora/pull/265
[294]: https://github.com/rmcolq/pandora/issues/294
[320]: https://github.com/rmcolq/pandora/issues/320
[v0.7.0]: https://github.com/rmcolq/pandora/releases/tag/v0.7.0
[Racon]: https://github.com/lbcb-sci/racon
[SAM]: https://samtools.github.io/hts-specs/SAMv1.pdf
