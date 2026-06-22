# How to publish this — the steps only you can do

I (the AI assistant) prepared everything below but **cannot create accounts, verify your email, or
click "submit"** — posting a preprint is a permanent public act under your name, so each step here is
yours. Estimated total time: Zenodo ~15 min, bioRxiv ~30–45 min. All free. **Do these in order.**

## 0. One-time prerequisites (10 min, do first)
- [x] **ORCID registered** — `0009-0000-5520-915X` (Seth Cope). Already written into `CITATION.cff`,
      `.zenodo.json`, `ABSTRACT.md`, and `paper/MANUSCRIPT.md`.
- [x] **Real name in all metadata** (`LICENSE`, `CITATION.cff`, `.zenodo.json`, `paper/MANUSCRIPT.md`,
      `ABSTRACT.md`) — done. Affiliation: "Independent researcher". **Do not list an institutional
      affiliation you aren't part of.**
- [x] **bioRxiv account created.** (Registered without the ORCID link — fine; add the ORCID later in
      bioRxiv account settings so it attaches to the preprint.)

## 1. Zenodo — get a citable DOI for the code (DONE)
- [x] Repo is **public**; Zenodo↔GitHub toggle ON; release **v1.0.0** published.
- [x] **DOI minted: `10.5281/zenodo.20799761`** (https://doi.org/10.5281/zenodo.20799761) — written
      into `paper/MANUSCRIPT.md` (Data availability), `CITATION.cff`, and `ABSTRACT.md`.

## 2. Manuscript PDF (DONE — no pandoc needed)
- [x] `paper/manuscript.html` is a self-contained, figure-embedded render (Figs 1 & 2 + the ATI supp
      are base64-embedded). **To make the PDF:** open `paper/manuscript.html` in Firefox →
      **File → Print → Save as PDF** (A4, default margins) → save as `paper/manuscript.pdf`.
      (Rebuild the HTML anytime with `/tmp/mdvenv/bin/python /tmp/build_manuscript_html.py`.)
- [x] Both figures (`analysis/p14_coupling_phase.png` Fig 1, `analysis/p16_analytic.png` Fig 2) are
      embedded in the HTML/PDF; you can also upload them to bioRxiv as separate files if it asks.
- [ ] **Re-read the whole thing in your own voice.** Fix anything that doesn't sound like you. Keep
      every caveat — the honesty is the credibility.

## 3. bioRxiv — post the preprint (SUBMITTED 2026-06-22)
- [x] Submitted: **MS ID `BIORXIV/2026/733776`**, in screening (~24–72 h). Filed as Biological
      research / "Research article with data", subject **Systems Biology**, license **CC-BY**, sole
      author Seth Cope (corresponding) + ORCID, no funding, no competing interests, AI disclosed.
      (Abstract web field required ASCII — em-dash `—` and middle-dot `·` rejected; Greek χ/κ ok.)
- [ ] **Wait for the "posted" email**, which gives the public preprint link + bioRxiv DOI.
- [ ] To edit before it posts: author area → "Manuscripts Undergoing Screening" → "Request Return of
      Manuscript". Do NOT resubmit (creates a duplicate).

## 4. After it's live
- [ ] Put the bioRxiv link in the GitHub README.
- [ ] Send the outreach emails (de Boer/Dodd first; then Meyerhans, Conway) — see chat drafts — each
      linking the bioRxiv preprint and the repo.
- [ ] Optionally submit a poster abstract to **SMB 2026** or apply for the **CROI New Researcher
      Scholarship**.

## Hard rules (protect yourself)
- **Never pay a journal/"publisher" that emails you offering to publish an HIV cure.** That is the
  predatory-journal trap; "HIV cure" is the bait. bioRxiv and Zenodo are free.
- **Frame it as a modeling hypothesis, never as a cure or a finding.** It generates a testable idea;
  it does not demonstrate anything in a person. Overclaiming is the fastest way to be dismissed.
- **Fix the factual flags before any email:** it is *Leor* Weinberger (HIV TIP, Gladstone), not Ariel;
  no human HIV-TIP trial exists (say "preclinical"); don't repeat the unconfirmed "VxBiosciences"
  partnership.
