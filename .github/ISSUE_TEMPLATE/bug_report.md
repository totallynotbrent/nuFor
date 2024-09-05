name: Bug report
description: Report a reproducible bug in nuFor
labels: [bug]
body:
  - type: textarea
    id: summary
    attributes:
      label: Summary
      description: What went wrong, and what did you expect instead?
    validations:
      required: true
  - type: textarea
    id: repro
    attributes:
      label: Reproduction
      description: Exact commands, case files, or input needed to reproduce.
    validations:
      required: true
  - type: textarea
    id: environment
    attributes:
      label: Environment
      description: OS, compiler (gfortran/rustc) versions, CMake, HDF5, git revision, thread count.
    validations:
      required: true
  - type: dropdown
    id: area
    attributes:
      label: Area
      options:
        - CLI
        - Web UI
        - Fortran kernel
        - Rust runtime
        - Build / CI
        - Data formats
        - Verification
        - Documentation
        - Other
  - type: input
    id: revision
    attributes:
      label: Git revision
      description: Output of `git describe --always` or the commit hash you tested.