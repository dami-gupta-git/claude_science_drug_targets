---
name: docs-style
description: Tone and writing conventions for any document produced in this project — target dossiers, summaries, reports, skill docs, README files, or any other written output. Load this skill whenever producing or editing a document, summary, or report.
---
## Guidelines
- Write sentences, not phrases
- Keep documents terse, professional, and concise
- No double negatives
- Avoid the use of em-dashes

## Strict Rules for TONE

Professional and descriptive. State data; leave the reader to judge significance.

Avoid these:

- Editorialising and value judgements — "the interesting part", "worth knowing
  about", "the row that settles it", "a tempting substitute".
- Self-narration about the development process — "an earlier pass had this
  wrong", "found in review", "I rejected this approach".
- Marketing register — "powerful", "seamless", "robust" as unquantified
  adjectives.
- Emphasis by adjective where a number would do. "13.5 A RMSD and 0.18
  recovery" carries more than "poor agreement".
- Technical-sounding labels standing in for the property that makes them true —
  "the helpers here are pure", "a robust parser", "an idempotent writer". Name
  the actual property instead: "none of these call the network", "re-running
  it on the same input overwrites rather than appends". A label can be true
  and still read as a claim about the code's quality rather than a fact about
  its behaviour; the reader should not have to take the label on faith.


Prefer the passive or impersonal voice for design rationale: "a LIGSITE-style
detector was evaluated and not adopted" over "I tried it and it did not work".

## Numbers and claims
- If any numbers are cited, a context needs to be provided for what they are.
- In summary documents like README.md, reduce how many numbers are shown.
- Avoid statements like '5 helper functions', or '2 files'. These change frequently, do not
hardcode.

## Prose vs tables
- The intro/summary paragraph is prose only: classification, dominant finding, in plain English. No statistics. It should be about 3 sentences long. If any numbers are included, explain what they are, don't just put them in brackets.
- Each section opens with one plain-English sentence stating the finding. Numbers follow after that, in a table or a supporting sentence — not in the opening sentence.