# Changelog

All notable BactiPipe application changes will be documented in this file.
Releases follow [Semantic Versioning](VERSIONING.md), and Git tags use the
`vMAJOR.MINOR.PATCH` format.

## Unreleased

## 0.1.0 - 2026-08-31

Initial formally versioned development release.

- Provides Illumina and Nanopore QC and assembly workflows.
- Provides cohort-level `relate` strain-relatedness analysis.
- Provides `detect` virulence and antimicrobial-resistance analysis.
- Records the authoritative BactiPipe application version in generated reports
  and downstream provenance.
- Exposes the application version through `bactipipe -v` and
  `bactipipe --version`.
