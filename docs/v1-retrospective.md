# Version 1 development retrospective

Version 1 started as a software and curation plan, then became a 12-object biological collection. The strongest decision was to keep the full selected-sample matrix beside each pocket. That made every reduction testable and exposed weak candidates early. Two studies were rejected for different reasons: one could not meet the preservation thresholds, while another did not match its stated annotation. Keeping those failures prevented quiet changes to the quality bar.

Factorial studies required more care than simple two-group comparisons. A study can support several useful objects, but each object needs its own contrast, sample subset, title, identifier and evidence. The overlap audit now distinguishes legitimate same-study reuse from duplicate content. This should remain a release gate as the collection grows.

Licensing also changed the collection. Public access did not provide enough support for redistributing every candidate source. Expression Atlas supplied a clear licence and provider-credit requirement, so v1 packages only material that passed that review. Link-only and rejected records remain visible without being presented as certified data.

The next collection release should prioritize three improvements:

1. Add read-level objects with paired-read identity checks, recorded strandedness and alignment quality evidence.
2. Add a second assay only after its object contract, task-specific metrics and reference checks are as explicit as the RNA-seq path.
3. Test installation and retrieval on more operating systems and record actual user friction before changing the command interface.

Post-release feedback may change that order. Any v2 proposal should cite observed use, preserve exact v1 records and introduce new evidence without moving the published v1 tag.
