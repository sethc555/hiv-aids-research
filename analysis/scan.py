#!/usr/bin/env python3
"""Semantic Scholar literature scan for the HIV-cure frontier.

Mirrors campus_twin/research/s2_scan.py conventions (urllib, x-api-key header,
429 backoff). Reads S2_API_KEY from the environment. Queries map onto the
non-standard-pathway list + the web-research gaps (CAR-T, block-and-lock,
consensus/timeline) + the "predict post-intervention control" modeling target.

  python3 scan.py            # runs all topics, writes scan_results.md
"""
import json, os, sys, time, urllib.parse, urllib.request

KEY = os.environ.get("S2_API_KEY", "")
BASE = "https://api.semanticscholar.org/graph/v1/paper/search"
FIELDS = ("title,year,authors,citationCount,influentialCitationCount,abstract,"
          "externalIds,venue,openAccessPdf,publicationTypes")
OUT = "/home/seth/dev/hiv-aids-research/analysis/scan_results.md"

# (stem, heading, query, year_filter, sort, limit)
#   sort: "recent" -> (year desc, citations desc);  "cite" -> citations desc
TOPICS = [
    ("post_treatment_control", "Post-treatment control & immune correlates (the prediction TARGET)",
     "HIV post-treatment control immune correlates remission ART interruption", "2017-2026", "recent", 15),
    ("vaccinal_effect_bnab", "bNAb vaccinal effect (CD8 T-cell control, 3BNC117/10-1074)",
     "HIV broadly neutralizing antibody vaccinal effect CD8 T cell 3BNC117 10-1074", "2018-2026", "recent", 15),
    ("ls_bnab_remission", "Long-acting LS-bNAbs & ATI remission (RIO-class)",
     "long-acting LS broadly neutralizing antibodies HIV remission treatment interruption", "2021-2026", "recent", 12),
    ("transplant_mechanism", "Allo-transplant cures: graft-versus-reservoir mechanism",
     "allogeneic stem cell transplant HIV remission graft versus reservoir CCR5 wild-type", "2019-2026", "recent", 12),
    ("adcc_reservoir", "ADCC / Fc-effector reservoir clearance",
     "antibody dependent cellular cytotoxicity HIV reservoir clearance Fc effector", "2018-2026", "recent", 10),
    ("reservoir_measurement", "Reservoir measurement: IPDA / intact proviral DNA",
     "HIV intact proviral DNA assay IPDA reservoir quantification", "2018-2026", "cite", 12),
    ("nfl_sequencing", "Near-full-length proviral sequencing & defectivation",
     "HIV near full length proviral sequencing intact defective reservoir", "2017-2026", "cite", 10),
    ("clonal_expansion", "Clonal expansion & integration sites (BACH2/STAT5B) — control FEATURES",
     "HIV reservoir clonal expansion integration site BACH2 STAT5B proliferation", "2014-2026", "cite", 12),
    ("latency_reversal", "Shock-and-kill / LRAs (is it dead?)",
     "HIV latency reversing agent shock and kill HDAC inhibitor reservoir clinical trial", "2014-2026", "cite", 12),
    ("block_and_lock", "Block-and-lock / dCA-Tat deep latency (WEB-PASS GAP)",
     "HIV block and lock didehydro-cortistatin Tat inhibitor deep latency", "2015-2026", "cite", 10),
    ("car_t_hiv", "Anti-HIV CAR-T cell therapy (WEB-PASS GAP)",
     "chimeric antigen receptor CAR T cell HIV reservoir cure", "2019-2026", "recent", 12),
    ("crispr_excision", "CRISPR provirus excision (EBT-101-class)",
     "CRISPR Cas9 excision HIV provirus gene editing functional cure", "2019-2026", "recent", 10),
    ("tips_dip", "Therapeutic interfering particles / DIPs (heterodox, single-dose)",
     "HIV therapeutic interfering particles defective interfering single administration", "2011-2026", "cite", 8),
    ("tat_latency_circuit", "Tat positive-feedback latency circuit / stochastic bistability (dynamical control)",
     "HIV Tat positive feedback latency stochastic bistable noise gene expression circuit", "2005-2026", "cite", 10),
    ("reservoir_dynamics_model", "Mathematical models of reservoir dynamics/decay",
     "mathematical model HIV latent reservoir decay dynamics clonal", "2013-2026", "cite", 10),
    ("elite_controllers", "Elite controllers: chromatin/integration of intact provirus",
     "HIV elite controller deep latency chromatin integration intact provirus", "2018-2026", "cite", 10),
    ("predict_rebound", "Predictors of time-to-rebound at ATI (the MODELING problem)",
     "predictors time to viral rebound analytical treatment interruption biomarker HIV", "2018-2026", "recent", 12),
    ("cure_roadmap", "Cure strategy reviews / roadmap & consensus (WEB-PASS GAP)",
     "HIV cure strategies functional cure remission roadmap review", "2023-2026", "recent", 12),
]


def search(q, year, limit, retries=4):
    params = {"query": q, "limit": min(limit * 2, 60), "fields": FIELDS}
    if year:
        params["year"] = year
    url = f"{BASE}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"x-api-key": KEY} if KEY else {})
    for a in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.load(r).get("data", []) or []
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(3 * (a + 1)); continue
            return [{"_err": f"HTTP {e.code}"}]
        except Exception as e:
            return [{"_err": str(e)}]
    return [{"_err": "rate-limited"}]


def srt(papers, mode):
    def y(p): return p.get("year") or 0
    def c(p): return p.get("citationCount") or 0
    if mode == "recent":
        return sorted(papers, key=lambda p: (y(p), c(p)), reverse=True)
    return sorted(papers, key=c, reverse=True)


def authstr(p):
    a = p.get("authors") or []
    names = [x.get("name", "") for x in a[:3]]
    if len(a) > 3:
        names.append("et al.")
    return ", ".join(names)


def link(p):
    ext = p.get("externalIds") or {}
    if ext.get("DOI"):
        return f"https://doi.org/{ext['DOI']}"
    if ext.get("PubMed"):
        return f"https://pubmed.ncbi.nlm.nih.gov/{ext['PubMed']}"
    return ""


def fmt(p):
    if p.get("_err"):
        return [f"- _query error: {p['_err']}_"]
    title = p.get("title", "Untitled")
    yr = p.get("year", "n.d.")
    ven = p.get("venue") or ""
    cc = p.get("citationCount") or 0
    icc = p.get("influentialCitationCount") or 0
    oa = "  ·  OA-PDF" if (p.get("openAccessPdf") or {}).get("url") else ""
    ab = (p.get("abstract") or "").strip().replace("\n", " ")
    ab = (ab[:420].rstrip() + "…") if len(ab) > 420 else (ab or "_no abstract_")
    L = link(p)
    head = f"- **{title}** ({yr})"
    meta = f"  — {authstr(p)} · _{ven}_ · **{cc:,} cites** (infl {icc}){oa}"
    out = [head, meta]
    if L:
        out.append(f"  {L}")
    out.append(f"  > {ab}")
    return out


def main():
    if not KEY:
        print("warning: S2_API_KEY not set — unauthenticated (will likely 429)", file=sys.stderr)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    seen, all_papers, sections = {}, [], []
    for stem, heading, q, year, sort, limit in TOPICS:
        print(f"querying: {heading[:60]}…", file=sys.stderr)
        papers = srt(search(q, year, limit), sort)[:limit]
        sections.append((heading, q, year, papers))
        for p in papers:
            pid = (p.get("externalIds") or {}).get("DOI") or p.get("title")
            if pid and pid not in seen and not p.get("_err"):
                seen[pid] = p
                all_papers.append(p)
        time.sleep(1.1)

    top_cite = sorted(all_papers, key=lambda p: p.get("citationCount") or 0, reverse=True)[:20]
    newest = sorted([p for p in all_papers if (p.get("year") or 0) >= 2025],
                    key=lambda p: (p.get("year") or 0, p.get("citationCount") or 0), reverse=True)[:25]

    lines = [f"# HIV-cure frontier — Semantic Scholar scan",
             f"_Generated {time.strftime('%Y-%m-%d')} via S2 Graph API · {len(TOPICS)} topics · "
             f"{len(all_papers)} unique papers._", "",
             "## ★ Most-cited across all topics (seminal anchors)", ""]
    for p in top_cite:
        lines += fmt(p)
    lines += ["", "## ☆ Newest (2025–2026) across all topics (the live frontier)", ""]
    for p in newest:
        lines += fmt(p)
    lines += ["", "---", ""]
    for heading, q, year, papers in sections:
        lines += [f"## {heading}", f"`{q}`  · years={year}", ""]
        for p in papers:
            lines += fmt(p)
        lines += [""]
    with open(OUT, "w") as f:
        f.write("\n".join(lines))
    print(f"\nwrote {OUT}  ({len(all_papers)} unique papers, {len(newest)} from 2025-26)", file=sys.stderr)


if __name__ == "__main__":
    main()
