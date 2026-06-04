# Nature Methods compliance checklist

Checked on 2026-06-04 against current official Nature Methods pages.

## Chosen article type

Recommended framing: Article, with an explicit note that the present draft is a proof-of-concept unless expanded to multi-dataset validation.

Rationale: Nature Methods states that an Article describes a novel method or tool and should include full technical description, strong validation, reproducibility, general applicability and potential for biological discovery. Spatial omics and computational methods are within scope.

Alternative framing: Brief Communication if the submission is narrowed to a technical critique and preliminary diagnostic. That format is much shorter and allows at most 2 display items.

## Article-format constraints used

- Abstract: up to 150 words, unreferenced.
- Main text: 3,000 words, up to 5,000 words at editorial discretion, excluding abstract, Methods, references and figure legends.
- Display items: up to 6 figures and/or tables.
- Structure: introduction without heading, Results, Discussion and Online Methods.
- Results and Online Methods should use topical subheadings; Discussion should not.
- References: Nature Methods typically recommends up to 50.
- Supplementary Information is allowed.

Official source: https://www.nature.com/nmeth/content

## Initial-submission requirements used

- Initial submission does not need special formatting if it is suitable for editorial assessment and peer review.
- PDF, Word or TeX/LaTeX initial submissions are accepted.
- Submission should include a manuscript file, cover letter and optional Supplementary Information.
- The manuscript file should include author names and affiliations unless double-blind review is chosen, enough methods detail for replication, and a reference list.
- Methods should include all elements necessary to interpret and replicate the results.
- LLM use should be documented in Methods when used.

Official sources:

- https://www.nature.com/nmeth/submission-guidelines/initial-formatting
- https://www.nature.com/nmeth/submission-guidelines/preparing-your-submission

## Writing guidance used

Nature Methods asks that submissions be clear to non-specialists, minimize unexplained jargon, define abbreviations, and make the background, rationale and main conclusions clear.

Official source: https://www.nature.com/nmeth/submission-guidelines/writing-and-language

## Current draft status

This package is formatted as a Nature Methods Article draft and contains 4 proposed display figures plus a graphical abstract. It is not yet a complete submission-ready Article because it currently relies on one public Atera WTA breast dataset.

## Evidence gaps to close before real submission

- Validate the uRNA-only COSTE diagnostic across at least 3-5 tissues or platforms.
- Compare against segmentation-aware methods such as Baysor-like transcript-aware segmentation or newer graph-based reassignment methods.
- Add synthetic segmentation perturbations to quantify sensitivity to over-segmentation, under-segmentation and transcript diffusion.
- Add runtime and memory benchmarks for uRNA-only COSTE on whole-transcriptome data.
- Deposit exact reproducibility code and a small test fixture, plus a DOI-backed release.
- Decide whether the final message is a new method extension of COSTE or a focused technical diagnostic for spatial transcriptomics QC.
