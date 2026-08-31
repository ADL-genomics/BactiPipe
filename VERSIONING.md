# BactiPipe versioning policy

BactiPipe uses Semantic Versioning in the form `MAJOR.MINOR.PATCH`. Release Git
tags add a leading `v`, for example `v0.1.0`. The authoritative application
version is defined once in `bactipipe/__version__.py`; package metadata, CLI
output, provenance, and reports must read that value.

While the major version is zero, the command and output contracts are still
being stabilized:

- `PATCH` releases contain compatible fixes, documentation, and internal
  changes with no intended analytical-contract change.
- `MINOR` releases add commands, options, output fields, or intentional
  analytical changes. Any pre-1.0 incompatibility must be highlighted in the
  changelog.
- `MAJOR` releases identify a stable-contract incompatible release after
  BactiPipe reaches 1.0.0.

Application, output-schema, tool, and database versions are separate. A
BactiPipe release does not imply a new CheckM, KmerFinder, typing, AMR, or
virulence database release.

## Release checklist

1. Confirm the web-integrated and standalone analytical packages are in sync.
2. Run the focused and contract test suites.
3. Update `bactipipe/__version__.py` and `CHANGELOG.md` in the same commit.
4. Verify `bactipipe -v`, `bactipipe --version`, and package metadata agree.
5. Create an annotated `vMAJOR.MINOR.PATCH` tag on the release commit.
6. Publish the identical commit and tag to the authorized GitHub remotes.
