# Changelog

All notable changes to this project will be documented in this file.

The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this
project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.9.1]

### Added
- `pandora` is now installable through `conda`;
- A script to archive the `pandora` repository with git submodules;

### Changed
- Improved the sample example so now we can assert that the output produced is the expected one;
- Changes to the build process that enables `pandora` to be compiled in the `conda` environment;

## [0.9.0]

### Changed
- Version bump from `0.9.0-rc2` to `0.9.0`.

## [0.9.0-rc2]

### Changed
- `pandora discover` now processes one sample at a time, but runs with several threads on the heavy tasks, i.e. when
mapping reads, finding candidate regions, and finding denovo variants. The result is that it now takes a lot less RAM to
run on multiple samples.

## [0.9.0-rc1]

### Changed
- `pandora discover` now receives read index files describing samples and reads, and discover denovo sequences in these samples.
  To improve performance on discovering denovo sequences on several samples, `pandora discover` is now multithreaded, but
  the performance is still the same as the previous version, i.e. each sample is processed in a single-threaded way;
- `pandora discover` output changed to a proprietary format. See [example](example) for the new output;
- `pandora` can now communicate with a [`make_prg` prototype](https://github.com/leoisl/make_prg) that is able to update PRGs
without needing to realign and remake the PRG. This provides major performance upgrades to running the full `pandora` pipeline
with denovo discovery enabled, and there is no need anymore to use a `snakemake` pipeline
(see [this example](example/run_pandora.sh) to how to run the full pipeline);
- We now use [musl libc](https://musl.libc.org/) instead of [Holy Build Box](https://github.com/phusion/holy-build-box)
to build a precompiled portable binary, removing the dependency on `OpenMP 4.0+` or `GCC 4.9+`, and `GLIBC`;

## [0.8.0]

### Added

- We now provide a script to build a portable precompiled binary as
  another option to run `pandora` easily. The portable binary is now
  provided with the release;
- `pandora` can now provide a meaningful stack trace in case of errors,
  to facilitate debugging (need to pass flag `-DPRINT_STACKTRACE` to
  `CMake`). Due to this, we now add debug symbols (`-g` flag) to every
  `pandora` build type, but this
  [does not impact performance](https://stackoverflow.com/a/39223245).
  The precompiled binary has this enabled.

### Changed

- We now use the [Hunter](https://github.com/cpp-pm/hunter) package
  manager, removing the requirement of having `ZLIB` and `Boost`
  system-wide installations;
- `GATB` is now a git submodule instead of an external project
  downloaded and compiled during compilation time. This means that when
  git cloning `pandora`, `cgranges` and `GATB` are also
  downloaded/cloned, and when preparing the build (running `cmake`),
  `Hunter` downloads and installs `Boost`, `GTest` and `ZLIB`. Thus we
  still need internet connection to prepare the build (running `cmake`)
  but not for compiling (running `make`).
- We now use a GATB fork that accepts a `ZLIB` custom installation;
- Refactored all thirdparty libraries (`cgranges`, `GATB`, `backward`,
  `CLI11`, `inthash`) into their own directory `thirdparty`.

### Fixed

- Refactored asserts into exceptions, and now `pandora` can be compiled
  correctly in the `Release` mode. The build process is thus able to
  create a more optimized binary, resulting in improved performance.
- Don't assume Nanopore reads are longer than loci [[#265][265]]



## [v0.7.0]

There is a significant amount of changes to the project between version
0.6 and this release. Only major things are listed here. Future releases
from this point will have their changes meticulously documented here.

### Added

- `discover` subcommand for de novo variant discovery [[#234][234]]
- many more tests

### Changed

- FASTA/Q files are now parsed with `klib` [[#223][223]]
- command-line interface is now overhauled with many breaking changes
  [[#224][224]]
- global genotyping has been made default [[#220][220]]
- Various improvements to VCF-related functions

### Fixed

- k-mer coverage underflow bug in `LocalPRG` [[#183][183]]

[Unreleased]: https://github.com/rmcolq/pandora/compare/0.9.1...HEAD
[0.9.1]: https://github.com/rmcolq/pandora/releases/tag/0.9.1
[0.9.0]: https://github.com/rmcolq/pandora/releases/tag/0.9.0
[0.9.0-rc2]: https://github.com/rmcolq/pandora/releases/tag/0.9.0-rc2
[0.9.0-rc1]: https://github.com/rmcolq/pandora/releases/tag/0.9.0-rc1
[0.8.0]: https://github.com/rmcolq/pandora/releases/tag/0.8.0
[v0.7.0]: https://github.com/rmcolq/pandora/releases/tag/v0.7.0

[183]: https://github.com/rmcolq/pandora/issues/183
[220]: https://github.com/rmcolq/pandora/pull/220
[223]: https://github.com/rmcolq/pandora/pull/223
[224]: https://github.com/rmcolq/pandora/pull/224
[234]: https://github.com/rmcolq/pandora/pull/234
[265]: https://github.com/rmcolq/pandora/pull/265


