# Repository guidance

## Verification scope

- Use the smallest relevant validation after ordinary code or
  documentation changes, such as a targeted test file or a single
  affected Quarto render.
- Do not routinely run the complete `devtools::test()`,
  `devtools::check()`, package build, vignette build, and pkgdown
  validation suite after each change.
- Run the full package and documentation validation suite without asking
  only when the package version is clearly being bumped for a release.
- For work that does not include a clear version bump, ask the user
  before running the full validation suite. An explicit user request to
  run it counts as approval.
- Aggregate related edits before validation and avoid repeating a
  successful expensive check unless later changes could affect what it
  covered.
