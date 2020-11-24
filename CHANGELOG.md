# Changelog

All notable changes to this project will be documented in this file.

The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this
project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/olivierlacan/keep-a-changelog/compare/v0.7.0...HEAD
[v0.7.0]: https://github.com/rmcolq/pandora/releases/tag/v0.7.0

[183]: https://github.com/rmcolq/pandora/issues/183
[220]: https://github.com/rmcolq/pandora/pull/220
[223]: https://github.com/rmcolq/pandora/pull/223
[224]: https://github.com/rmcolq/pandora/pull/224
[234]: https://github.com/rmcolq/pandora/pull/234
