# Changelog

All notable BactiPipe application changes will be documented in this file.
Releases follow [Semantic Versioning](VERSIONING.md), and Git tags use the
`vMAJOR.MINOR.PATCH` format.

## Unreleased

- Isolate per-sample QC, assembly, and organism-identification failures so other batch samples can finish.
- Add repeatable `--only-sample` selectors, transactional `--update-existing`
  publication, and `sample_failures.tsv` for targeted QC reanalysis.
- Keep Nanopore sample-folder/barcode resolution aligned with the web launcher.
- Record deterministic application source and formal database identities in downstream result provenance.
- Include the formal BactiPipe database identity in QC and downstream reports.

## 0.1.0 - 2026-08-31

Initial formally versioned development release.

- Provides Illumina and Nanopore QC and assembly workflows.
- Provides cohort-level `relate` strain-relatedness analysis.
- Provides `detect` virulence and antimicrobial-resistance analysis.
- Records the authoritative BactiPipe application version in generated reports
  and downstream provenance.
- Exposes the application version through `bactipipe -v` and
  `bactipipe --version`.
