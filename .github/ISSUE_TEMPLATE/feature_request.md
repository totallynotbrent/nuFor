name: Feature request
description: Propose work against the roadmap or a new idea
labels: [enhancement]
body:
  - type: textarea
    id: problem
    attributes:
      label: Problem or gap
      description: What capability is missing or weak? Reference a roadmap step or spec section if relevant.
    validations:
      required: true
  - type: textarea
    id: proposal
    attributes:
      label: Proposed change
      description: Sketch of the approach and what files/area it touches.
    validations:
      required: true
  - type: textarea
    id: verification
    attributes:
      label: Verification
      description: How should this be tested or verified before it counts as done?
    validations:
      required: true
  - type: dropdown
    id: scope
    attributes:
      label: Scope
      options:
        - Year 1 foundations
        - Year 2 2D / performance
        - Year 3+ architecture
        - Tooling / DX
        - Documentation
        - Not sure