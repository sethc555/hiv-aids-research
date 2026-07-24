#!/usr/bin/env python3
"""Second-pass S2 analysis for the HIV-cure corpus.

Three jobs (the user's "do all nexts"):
  1. citation_graph.md — references (parents) + citations (children) of the
     anchor papers, with an aggregated "active labs" ranking and the 2025-26
     live follow-on work.
  2. bibliography.md   — the topic scan, venue/citation-gated (predatory review
     mills removed), deduped, themed: a clean foundation doc.
Raw S2 JSON is cached to disk so re-runs are free (self-curating corpus).

Reads S2_API_KEY from env. Mirrors scan.py conventions (urllib, x-api-key, 429 backoff).
"""
import os
import json, os, sys, time, urllib.parse, urllib.request

KEY = os.environ.get("S2_API_KEY", "")
G = "https://api.semanticscholar.org/graph/v1"
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "raw_cache.json")

# ---- topic scan (same 18 as scan.py) -------------------------------------
TOPICS = [
    ("post-treatment control & immune correlates", "HIV post-treatment control immune correlates remission ART interruption", "2017-2026"),
    ("bNAb vaccinal effect", "HIV broadly neutralizing antibody vaccinal effect CD8 T cell 3BNC117 10-1074", "2018-2026"),
    ("long-acting LS-bNAb remission", "long-acting LS broadly neutralizing antibodies HIV remission treatment interruption", "2021-2026"),
    ("allo-transplant mechanism", "allogeneic stem cell transplant HIV remission graft versus reservoir CCR5 wild-type", "2019-2026"),
    ("ADCC / Fc-effector clearance", "antibody dependent cellular cytotoxicity HIV reservoir clearance Fc effector", "2018-2026"),
    ("reservoir measurement (IPDA)", "HIV intact proviral DNA assay IPDA reservoir quantification", "2018-2026"),
    ("near-full-length sequencing", "HIV near full length proviral sequencing intact defective reservoir", "2017-2026"),
    ("clonal expansion / integration sites", "HIV reservoir clonal expansion integration site BACH2 STAT5B proliferation", "2014-2026"),
    ("shock-and-kill / LRAs", "HIV latency reversing agent shock and kill HDAC inhibitor reservoir clinical trial", "2014-2026"),
    ("block-and-lock / dCA-Tat", "HIV block and lock didehydro-cortistatin Tat inhibitor deep latency", "2015-2026"),
    ("anti-HIV CAR-T", "chimeric antigen receptor CAR T cell HIV reservoir cure", "2019-2026"),
    ("CRISPR provirus excision", "CRISPR Cas9 excision HIV provirus gene editing functional cure", "2019-2026"),
    ("TIPs / DIPs", "HIV therapeutic interfering particles defective interfering single administration", "2011-2026"),
    ("Tat latency circuit (dynamical)", "HIV Tat positive feedback latency stochastic bistable noise gene expression circuit", "2005-2026"),
    ("reservoir dynamics models", "mathematical model HIV latent reservoir decay dynamics clonal", "2013-2026"),
    ("elite controllers", "HIV elite controller deep latency chromatin integration intact provirus", "2018-2026"),
    ("predictors of time-to-rebound", "predictors time to viral rebound analytical treatment interruption biomarker HIV", "2018-2026"),
    ("cure roadmap / consensus", "HIV cure strategies functional cure remission roadmap review", "2023-2026"),
]

# ---- citation-graph anchors ----------------------------------------------
ANCHORS = [
    ("TIPs — NHP proof-of-concept (Pitchai 2024 Science)", "DOI:10.1126/science.adn5866"),
    ("TIPs — platform (Chaturvedi 2021 Cell)",             "DOI:10.1016/j.cell.2021.11.004"),
    ("Reservoir dynamics model (Barbehenn 2024 Nat Commun)","DOI:10.1038/s41467-024-54116-1"),
    ("ATI meta-analysis (Gunst 2025 Nat Commun)",          "DOI:10.1038/s41467-025-56116-1"),
    ("Tat-circuit fate decision (Weinberger 2007 Nat Genet)","DOI:10.1038/ng.116"),
]

PREDATORY = ["iaa journal", "iaajas", "iaajsr", "idosr", "inosr", "nijpp", "rojphm",
             "research output journal", "newport international", "integralize",
             "biocaster", "jscrte", "bioscientia", "science insights",
             "international integralize", "iaa j", "inosr "]
REPUTABLE = ["nature", "science translational", "science advances", "science immunology",
             "cell", "immunity", "lancet", "pnas", "proceedings of the national academy",
             "journal of clinical investigation", "jci insight", "plos", "mbio",
             "journal of virology", "journal of infectious", "clinical infectious",
             "ebiomedicine", "frontiers in", "viruses", "clinical chemistry",
             "communications medicine", "scientific reports", "methods", "virology journal",
             "current opinion in hiv", "aids", "retrovirology", "elife", "blood",
             "chemical reviews", "nature genetics", "nature biomedical", "value in health",
             "pathogens and immunity", "mucosal immunology", "cellular & molecular",
             "journal of medical virology", "ieee", "journal of biological dynamics",
             "mathematical biosciences", "virus evolution", "star protocols",
             "journal of virus eradication", "cells", "journal of immunology", "vaccines",
             "archives of virology", "pharmaceutics", "topics in antiviral",
             "biorxiv", "medrxiv", "research square", "tropical medicine"]
PREPRINTS = ["biorxiv", "medrxiv", "research square"]


def _get(url):
    req = urllib.request.Request(url, headers={"x-api-key": KEY} if KEY else {})
    for a in range(5):
        try:
            with urllib.request.urlopen(req, timeout=45) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(3 * (a + 1)); continue
            return {"_err": f"HTTP {e.code}", "data": []}
        except Exception as e:
            return {"_err": str(e), "data": []}
    return {"_err": "rate-limited", "data": []}


def search(query, year, limit=24):
    f = "title,year,authors,citationCount,influentialCitationCount,abstract,externalIds,venue,openAccessPdf"
    p = urllib.parse.urlencode({"query": query, "year": year, "limit": limit, "fields": f})
    return _get(f"{G}/paper/search?{p}").get("data", []) or []


def graph(anchor_id, kind):  # kind: 'citations' | 'references'
    f = "title,year,authors,venue,citationCount,externalIds"
    p = urllib.parse.urlencode({"fields": f, "limit": 200})
    d = _get(f"{G}/paper/{anchor_id}/{kind}?{p}").get("data", []) or []
    key = "citingPaper" if kind == "citations" else "citedPaper"
    return [x[key] for x in d if x.get(key)]


def load_cache():
    if os.path.exists(CACHE):
        with open(CACHE) as f:
            return json.load(f)
    return {"topics": {}, "graph": {}}


def venue_ok(v):
    v = (v or "").lower()
    if any(p in v for p in PREDATORY):
        return False, "predatory"
    if any(r in v for r in REPUTABLE):
        return True, ("preprint" if any(p in v for p in PREPRINTS) else "ok")
    return None, "other"


def keep(p):
    ok, tag = venue_ok(p.get("venue"))
    cc = p.get("citationCount") or 0
    if ok is False:
        return False, tag
    if ok is True:
        return True, tag
    return (cc >= 15), ("highcite" if cc >= 15 else "other-lowcite")


def auth(p, n=3):
    a = [x.get("name", "") for x in (p.get("authors") or [])]
    return ", ".join(a[:n]) + (" et al." if len(a) > n else "")


def doi(p):
    e = p.get("externalIds") or {}
    return f"https://doi.org/{e['DOI']}" if e.get("DOI") else (
        f"https://pubmed.ncbi.nlm.nih.gov/{e['PubMed']}" if e.get("PubMed") else "")


def line(p, tag=""):
    cc, icc = p.get("citationCount") or 0, p.get("influentialCitationCount") or 0
    t = f" _[{tag}]_" if tag == "preprint" else ""
    L = f"\n  {doi(p)}" if doi(p) else ""
    return f"- **{p.get('title','?')}** ({p.get('year','n.d.')}){t} — {auth(p)} · _{p.get('venue','')}_ · {cc:,} cites (infl {icc}){L}"


def main():
    if not KEY:
        print("warning: no S2_API_KEY", file=sys.stderr)
    cache = load_cache()
    fetched = 0

    # ---- 1. topic scan (cached) ----
    all_papers, seen = [], {}
    for label, q, year in TOPICS:
        if q not in cache["topics"]:
            print(f"fetch topic: {label}", file=sys.stderr)
            cache["topics"][q] = search(q, year); fetched += 1; time.sleep(1.1)
        for p in cache["topics"][q]:
            k = (p.get("externalIds") or {}).get("DOI") or p.get("title")
            if k and k not in seen:
                seen[k] = p; p["_topic"] = label; all_papers.append(p)

    # ---- 2. citation graph (cached) ----
    for label, aid in ANCHORS:
        for kind in ("references", "citations"):
            ck = f"{aid}::{kind}"
            if ck not in cache["graph"]:
                print(f"fetch graph: {label} {kind}", file=sys.stderr)
                cache["graph"][ck] = graph(aid, kind); fetched += 1; time.sleep(1.1)

    with open(CACHE, "w") as f:
        json.dump(cache, f)
    print(f"fetched {fetched} new calls; cache has "
          f"{len(cache['topics'])} topics, {len(cache['graph'])} graph sets", file=sys.stderr)

    # ===== bibliography.md =====
    kept, dropped = [], {"predatory": 0, "other-lowcite": 0}
    for p in all_papers:
        ok, tag = keep(p)
        (kept.append((p, tag)) if ok else dropped.__setitem__(tag, dropped.get(tag, 0) + 1))
    by_topic = {}
    for p, tag in kept:
        by_topic.setdefault(p["_topic"], []).append((p, tag))

    bib = [f"# HIV-cure — filtered bibliography (foundation doc)",
           f"_S2 Graph API · {time.strftime('%Y-%m-%d')} · {len(all_papers)} scanned → "
           f"**{len(kept)} kept**, {dropped['predatory']} predatory-venue + "
           f"{dropped.get('other-lowcite',0)} low-cite-obscure dropped._",
           "", "_Venue-gated to reputable journals/preprints OR ≥15 citations. "
           "Preprints tagged. Grouped by theme (first match)._", ""]
    for label, _, _ in TOPICS:
        ps = by_topic.get(label, [])
        if not ps:
            continue
        ps.sort(key=lambda x: (x[0].get("year") or 0, x[0].get("citationCount") or 0), reverse=True)
        bib.append(f"## {label}  ({len(ps)})")
        for p, tag in ps:
            bib.append(line(p, tag))
        bib.append("")
    with open(os.path.join(HERE, "bibliography.md"), "w") as f:
        f.write("\n".join(bib))

    # ===== citation_graph.md =====
    cg = [f"# HIV-cure — citation-graph expansion of anchors",
          f"_S2 Graph API · {time.strftime('%Y-%m-%d')}_", ""]
    lab_counts = {}
    for label, aid in ANCHORS:
        refs = cache["graph"].get(f"{aid}::references", [])
        cites = cache["graph"].get(f"{aid}::citations", [])
        recent = sorted([c for c in cites if (c.get("year") or 0) >= 2025],
                        key=lambda c: (c.get("year") or 0, c.get("citationCount") or 0), reverse=True)
        for c in cites:  # count active labs from CHILDREN (who builds on it)
            for a in (c.get("authors") or []):
                nm = a.get("name", "")
                if nm:
                    lab_counts[nm] = lab_counts.get(nm, 0) + 1
        cg += [f"## {label}",
               f"`{aid}` — {len(refs)} references (parents), {len(cites)} citations (children), "
               f"{len(recent)} children from 2025-26.", "",
               f"### Live follow-on (2025-26 children)"]
        cg += [line(c) for c in recent[:14]] or ["- _none indexed yet_"]
        top_refs = sorted(refs, key=lambda r: r.get("citationCount") or 0, reverse=True)[:6]
        cg += ["", "### Most-cited intellectual parents (references)"]
        cg += [line(r) for r in top_refs] or ["- _none_"]
        cg += [""]
    top_labs = sorted(lab_counts.items(), key=lambda kv: kv[1], reverse=True)[:30]
    cg += ["## ★ Most active authors across all anchor-citing papers",
           "_(frequency = how often they appear among papers building on the anchors → active labs)_", ""]
    cg += [f"- **{nm}** — {n}" for nm, n in top_labs if n >= 2]
    with open(os.path.join(HERE, "citation_graph.md"), "w") as f:
        f.write("\n".join(cg))

    print(f"\nwrote bibliography.md ({len(kept)} papers) and citation_graph.md", file=sys.stderr)
    print(f"dropped: {dropped}", file=sys.stderr)


if __name__ == "__main__":
    main()
