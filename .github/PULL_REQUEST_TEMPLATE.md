<!--
Thank you for contributing to theiagene!

Fill in the sections below and delete anything that doesn't apply.
-->

<!-- Indicate the issue number if applicable; otherwise, delete -->
This PR closes #

<!-- Delete the < > around "NOT" if your branch should be retained after merging -->
🗑️ This dev branch should <NOT> be deleted after merging to main.

## :brain: Summary

<!-- What does this PR do, and why? A paragraph or two is plenty. -->

## :hammer_and_wrench: Technical

### Changes
<!-- Describe your changes. Bullets are fine. Call out anything a reviewer would
     otherwise have to reverse-engineer from the diff. -->

### Behavior changes
<!-- Could this change the output of an existing command? Anything that alters
     coverage values, which variants are extracted, report line formatting, or
     output file names/columns needs an explicit callout here, since downstream
     WDL tasks parse these outputs. -->

This PR may lead to different results for existing workflows: **Yes/No**

This PR is backwards incompatible (CLI arguments, output files, or library API): **Yes/No**

### Dependencies
<!-- New or bumped dependencies in pyproject.toml? Write "None" if unchanged. -->

## :test_tube: Testing

<!-- Describe how you tested this. Include the commands you ran and the input
     data (test fixtures, a real BAM/VCF/GFF, a VEP-annotated VCF, etc.). -->

### Verification
<!-- Mark each box [X]. Where a check genuinely doesn't apply, write "N/A —
     <reason>" next to it rather than leaving it unchecked. -->

- [ ] **Pytests** — `pytest tests/ --cov=src/theiagene --cov-report=term-missing` passes locally, and the *Theiagene Tests* check is green on this PR
  <!-- That workflow only triggers on changes under src/theiagene/**, tests/**,
       conftest.py, or pyproject.toml. If it didn't run for this PR, say so here. -->
- [ ] **Test coverage** — new or updated tests added under `tests/` for the changed behavior
- [ ] **Version bump** — `version` bumped in **both** `pyproject.toml` and `src/theiagene/__init__.py`, or intentionally left alone (state why)
  <!-- The Docker workflow fails if the release tag doesn't match pyproject.toml,
       and `theiagene --version` reads __init__.py, so the two must agree. -->
- [ ] **Docker build** — `docker build --target test .` succeeds (runs pytest inside the image) and `docker build -t theiagene:pr . && docker run --rm theiagene:pr theiagene --version` prints the expected version
- [ ] **Live run** — exercised on real data or inside the consuming WDL task, where applicable (describe above)

### Suggested scenarios for the reviewer to test
<!-- Edge cases, test data, or commands the reviewer should try. -->

## :microscope: Final Developer Checklist

- [ ] The change has been run and the outputs are as anticipated
- [ ] No credentials or non-public sample data were committed (including in test fixtures)
- [ ] `README.md` has been updated for any user-facing change (new subcommands, arguments, or output formats)
- [ ] All CI checks are passing

## 🎯 Reviewer Checklist
<!-- Indicate NA when not applicable -->

- [ ] The change does what the summary claims and the approach is sound
- [ ] Test coverage is adequate for the change
- [ ] Confirmed that changes work as expected with existing workflows
- [ ] Documentation is accurate
- [ ] The PR author has addressed all comments
