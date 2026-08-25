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

## Verify the prose against the table before writing it

The kernels compute correctly; the summarising step is where documents go wrong.
Every error found in review of these runs has been a hand-written sentence that
the accompanying CSV did not support, never a wrong number in the CSV. Treat the
prose as the untrusted part.

Four rules, each addressing a failure that has actually occurred:

- **Derive every count, never type it.** A sentence saying "49 of 54 tissues are
  typical" or "five lineages reach significance" must come from
  `df.level.value_counts()` or `int((df.p < 0.05).sum())` evaluated in the same
  cell that writes the text. A tally read off a printed table by eye will
  sometimes be the count of the rows that happened to be visible.
- **Derive every ranking, and state what was excluded.** "The highest ranking
  tissues are X, Y, Z" must be a slice of the sorted frame. If the ranking
  excludes some rows — cell lines among tissues, untested lineages, off-scope
  papers — the exclusion is part of the claim and belongs in the sentence.
- **Never write a path or identifier you have not printed.** A run directory, a
  timestamp, an accession or a PMID goes into prose only by reading the variable
  that holds it. A plausible-looking timestamp typed from memory will be wrong
  and will point a reader at a directory that does not exist.
- **Count the sources that bear on the claim, not the sources in the run.** "All
  five sources agree" is false when two of them are silent on the point. Name
  the sources that speak to it and say the others do not.

Two habits behind these failures are worth naming, because they recur:

- A claim that makes two independent results line up neatly is the one most
  likely to have been trimmed to fit. When a ranking, count or overlap appears to
  corroborate another run, re-derive it before writing the sentence.
- A qualitative assertion is as checkable as a numeric one. "Already used
  clinically", "well characterised", "routinely tested" either trace to a run in
  the report or come out. Absence of a curated record is not evidence, and
  neither is background knowledge that no run supports.

For a report collating several runs, write the self-audit as a table: one row per
claim, with the value in the prose, the value re-derived from source, and the
verdict. Ship it beside the report. It is cheap, it catches this class of error
before a reader does, and a claim that cannot be expressed as a row in it is
usually a claim that should not be in the report.