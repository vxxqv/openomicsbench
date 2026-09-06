# Contributing to OpenOmicsBench

Start with a small proposal: name the source accession, the intended task and the reason a compact version would help. A useful contribution can be a corrected sample label or a clearer explanation; it does not need to be a new dataset.

## Intake before processing

Create a record under `curation/intake` before downloading large files. Include the source URL, publication, access date, assay, study design, sample count, available layers and estimated download size. Leave unverified fields null and explain what needs checking. Do not assign a numerical score to a design that has only been read at series level.

Record the redistribution decision with dated evidence. GREEN requires specific support for the material being packaged. AMBER permits a source link and project-written recipe while the data remain outside the release. RED excludes the source. RECHECK pauses previously cleared material until the changed terms have been reviewed. Public download access alone is insufficient.

## Build a reviewable object

Use the stable ID assigned during curation. Resolve every source file to a version or accession, record its hash and preserve exact sample order. Distinguish biological replicates from repeated libraries. Keep raw counts, normalized expression and quantification estimates in separate files with explicit roles.

Freeze the design, reduction algorithm, seed, feature universe and candidate grid before evaluating reductions. Retain the full candidate curve, including failed sizes. Explain why the smallest passing size is useful for its task. Do not select genes using a desired publication conclusion or manually remove awkward samples to improve a figure.

Run the baseline from a pinned environment. Record the code revision, package versions, resolved parameters and file hashes. Keep source and transformed files separate. For paired reads, check both mate count and read identity. For a count matrix, verify labels and exact integer values after selection.

## Review and tests

Run the core tests and repository checker. Add a regression test for a new failure mode; avoid tests that only repeat an assignment in the implementation. Schema changes must regenerate the JSON Schema and explain how old manifests migrate. A biological adapter also needs a small successful regeneration and a meaningful expected-output check.

Have another person follow the dataset card from a clean environment. Record the command, outcome and any correction needed, with their consent. An automated agent run does not satisfy the independent-user gate.

## Writing and comments

Write as a maintainer explaining the files to another person. Name the dataset, state the task, give a command that works and explain a limitation where it affects use. Remove promotional claims and repeated boilerplate. Do not use em dashes in documentation, help text or release notes.

Comments should explain a reason or a constraint that the code does not make obvious. Keep public docstrings concise. There is no maximum comment length; a necessary explanation is better than an unexplained trick.

Dataset cards need a source citation, included files and sizes, the transformation, expected behavior, specific limitations and the exact release citation once it exists. Release notes should explain changes by user impact. Do not publish a placeholder DOI, author or successful validation result.
