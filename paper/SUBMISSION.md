# How to publish this — the steps only you can do

I (the AI assistant) prepared everything below but **cannot create accounts, verify your email, or
click "submit"** — posting a preprint is a permanent public act under your name, so each step here is
yours. Estimated total time: Zenodo ~15 min, bioRxiv ~30–45 min. All free. **Do these in order.**

## 0. One-time prerequisites (10 min, do first)
- [ ] **Register a free ORCID** at https://orcid.org — a researcher ID. Add it to `CITATION.cff` and
      `paper/MANUSCRIPT.md` (replace the TODOs). Most preprint servers ask for it.
- [ ] **Put your real name** in `LICENSE`, `CITATION.cff`, `.zenodo.json`, and `paper/MANUSCRIPT.md`
      (they currently say "Seth C."). Decide the affiliation you'll use — "Independent researcher" is
      fine and respected; if you want a citable institutional-style affiliation, the Ronin Institute
      (~$100/yr) provides one. **Do not list a Meharry lab affiliation you aren't part of.**

## 1. Zenodo — get a citable DOI for the code (do this BEFORE bioRxiv, ~15 min)
- [ ] Go to https://zenodo.org → **Log in with GitHub** (top right) → authorize.
- [ ] Menu → **GitHub** → find `sethc555/hiv-aids-research` → flip the toggle **ON**.
- [ ] Back on GitHub, make a release: repo → **Releases** → **Draft a new release** → tag `v1.0.0`,
      title it, publish. (`.zenodo.json` in the repo auto-fills the metadata.)
- [ ] Zenodo auto-archives the release and issues a **DOI** within a minute or two. Copy it.
- [ ] Paste that DOI into `paper/MANUSCRIPT.md` (the "DOI: TODO" in Data availability) and re-release,
      or just cite it in the bioRxiv submission.

## 2. Manuscript PDF (5 min)
- [ ] Convert the draft to PDF: `pandoc paper/MANUSCRIPT.md -o paper/manuscript.pdf`
      (install pandoc + a LaTeX engine, or paste the markdown into Google Docs and export PDF).
- [ ] Make sure the two figures are embedded or supplied separately:
      `analysis/p14_coupling_phase.png` (Fig 1) and a plot of `analysis/p16_analytic.npz` (Fig 2).
- [ ] **Re-read the whole thing in your own voice.** Fix anything that doesn't sound like you. Keep
      every caveat — the honesty is the credibility.

## 3. bioRxiv — post the preprint (~30 min)
- [ ] Go to https://www.biorxiv.org → **Submit** → create an account (personal email is OK; an
      institutional email just speeds the automated check).
- [ ] **Subject area:** *Systems Biology* (or *Microbiology*). NOT medRxiv — this is a modeling study.
- [ ] **Authors/affiliation:** your name + ORCID; affiliation "Independent researcher" (or Ronin).
- [ ] **License:** CC-BY-4.0 (standard, lets people cite/reuse) — or CC-BY-NC if you prefer.
- [ ] Upload `manuscript.pdf` + figure files.
- [ ] **Competing interests:** none. **Funding:** none.
- [ ] **Disclose the AI assistance** in the comments/methods (the manuscript already does).
- [ ] Submit. Screening (not peer review) takes ~24–48 h; then it's public with a DOI.

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
