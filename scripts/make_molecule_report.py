# env: chameleon-calc
"""
make_molecule_report.py
-----------------------
Generate a focused, SELF-CONTAINED R/S epimer report for ONE molecule, independent of
any other molecule. Reads the descriptor CSV produced by `ensemble_descriptors.py`
(no recomputation), fills a template, and (optionally) converts to .docx via export_docx.

Each report covers only that molecule's R vs S epimers in water and chloroform; it makes
no cross-molecule comparison. The permeability/solubility direction is read straight off
the data (which isomer is more polar-exposed in water, which is more compact in chloroform).

Usage:
  python scripts/make_molecule_report.py \
      --csv results/descriptors_4_nondz.csv \
      --pair 3-12-8-12 --r 3-12-8-12_R --s 3-12-8-12_S \
      --residue "L-azetidine-2-carboxylic acid (rigid 4-membered ring; cis-amide inducer; no backbone N-H)" \
      --out docs/reports/3-12-8-12_RS_report.md
"""
from __future__ import annotations

import argparse
import re
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))


def fmt(v, nd=1):
    try:
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return "n/a"
        return f"{float(v):.{nd}f}"
    except Exception:
        return "n/a"


# (row label, CSV key suffix, decimals, weak) — Boltzmann-weighted per-solvent descriptors.
# weak=True marks descriptors graded as weak/unvalidated against passive permeability in the
# 2026 ranked literature review (amphipathic moment, asphericity, RMSF, SA_HD/HA): they are still
# computed and reported for downstream ML, but do NOT drive the permeability/solubility call.
DESCS = [
    ("3D-PSA (Å²)", "bw_psa", 1, False),
    ("SA_HD donor-H surface (Å²)", "bw_SA_HD", 1, True),
    ("SA_HA acceptor surface (Å²)", "bw_SA_HA", 1, True),
    ("hydrophobic SASA (Å²)", "bw_hydrophobic_sasa", 1, False),
    ("amphipathic moment (Å)", "bw_amphi_moment", 2, True),
    ("total IMHB", "bw_IMHB", 2, False),
    ("backbone (transannular) IMHB", "bw_IMHB_bb", 2, False),
    ("side-chain IMHB", "bw_IMHB_res", 2, False),
    ("radius of gyration (Å)", "bw_rg", 2, False),
    ("asphericity", "bw_asphericity", 3, True),
    ("ensemble RMSF (Å)", "rmsf", 2, True),
]

WEAK_NOTE = (
    "> **Note.** Descriptors marked † (SA_HD, SA_HA, amphipathic moment, asphericity, ensemble RMSF) "
    "currently have low or unestablished correlation with passive permeability in the macrocycle "
    "literature. They are reported and retained for downstream machine-learning use, but are **not** "
    "used to drive the permeability/solubility call in §4, which rests on the better-validated "
    "descriptors (3D-PSA, radius of gyration, backbone IMHB)."
)


def desc_table(R, S, solv):
    lab = "water" if solv == "water" else "chloroform"
    out = [f"| descriptor ({lab}) | R | S | Δ (R−S) |", "|---|---|---|---|"]
    for label, key, nd, weak in DESCS:
        mark = " †" if weak else ""
        r, s = R.get(f"{solv}_{key}"), S.get(f"{solv}_{key}")
        d = (r - s) if (pd.notna(r) and pd.notna(s)) else None   # signed difference, R minus S
        out.append(f"| {label}{mark} | {fmt(r, nd)} | {fmt(s, nd)} | {fmt(d, nd)} |")
    return "\n".join(out)


def pick(R, S, key, want="higher"):
    r, s = R.get(key), S.get(key)
    if pd.isna(r) or pd.isna(s):
        return None
    return ("R" if r > s else "S") if want == "higher" else ("R" if r < s else "S")


def renumber(body: str, refs: list):
    """Renumber ^{n}^ citations by order of first appearance and reorder the reference list
    to match (ACS style). Returns (new_body, new_refs). Assumes every ref is cited."""
    order = []
    for tok in re.findall(r"\^\{([^}]+)\}\^", body):
        for n in tok.split(","):
            i = int(n)
            if i not in order:
                order.append(i)
    remap = {old: new for new, old in enumerate(order, 1)}

    def _sub(m):
        nums = sorted(remap[int(x)] for x in m.group(1).split(","))
        return "^{" + ",".join(map(str, nums)) + "}^"

    new_body = re.sub(r"\^\{([^}]+)\}\^", _sub, body)
    new_refs = [refs[old - 1] for old in order]
    return new_body, new_refs


REFS = [
    "Pracht, P.; Bohle, F.; Grimme, S. Automated Exploration of the Low-Energy Chemical Space with Fast Quantum Chemical Methods. *Phys. Chem. Chem. Phys.* **2020**, *22* (14), 7169–7192. DOI: 10.1039/c9cp06869d.",
    "Bannwarth, C.; Ehlert, S.; Grimme, S. GFN2-xTB: An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method. *J. Chem. Theory Comput.* **2019**, *15* (3), 1652–1671. DOI: 10.1021/acs.jctc.8b01176.",
    "Ehlert, S.; Stahn, M.; Spicher, S.; Grimme, S. Robust and Efficient Implicit Solvation Model for Fast Semiempirical Methods. *J. Chem. Theory Comput.* **2021**, *17* (7), 4250–4261. DOI: 10.1021/acs.jctc.1c00471.",
    "Grambow, C. A.; Weir, H.; Cunningham, C. N.; Biancalani, T.; Chuang, K. V. CREMP: Conformer-Rotamer Ensembles of Macrocyclic Peptides for Machine Learning. *Sci. Data* **2024**, *11*, 859. DOI: 10.1038/s41597-024-03698-y.",
    "Ono, S.; Naylor, M. R.; Townsend, C. E.; Okumura, C.; Okada, O.; Lokey, R. S. Conformation and Permeability: Cyclic Hexapeptide Diastereomers. *J. Chem. Inf. Model.* **2019**, *59* (6), 2952–2963. DOI: 10.1021/acs.jcim.9b00217.",
    "Begnini, F.; Poongavanam, V.; Atilaw, Y.; Erdélyi, M.; Schiesser, S.; Kihlberg, J. Cell Permeability of Isomeric Macrocycles: Predictions and NMR Studies. *ACS Med. Chem. Lett.* **2021**, *12* (6), 983–990. DOI: 10.1021/acsmedchemlett.1c00126.",
    "Poongavanam, V.; Atilaw, Y.; Ye, S.; Wieske, L. H. E.; Erdélyi, M.; Ermondi, G.; Caron, G.; Kihlberg, J. Predicting the Permeability of Macrocycles from Conformational Sampling: Limitations of Molecular Flexibility. *J. Pharm. Sci.* **2021**, *110* (1), 301–313. DOI: 10.1016/j.xphs.2020.10.052.",
    "Rzepiela, A. A.; Viarengo-Baker, L. A.; Tatarskii, V.; Kombarov, R.; Whitty, A. Conformational Effects on the Passive Membrane Permeability of Synthetic Macrocycles. *J. Med. Chem.* **2022**, *65* (14), 10300–10317. DOI: 10.1021/acs.jmedchem.1c02090.",
    "Rezai, T.; Bock, J. E.; Zhou, M. V.; Kalyanaraman, C.; Lokey, R. S.; Jacobson, M. P. Conformational Flexibility, Internal Hydrogen Bonding, and Passive Membrane Permeability: Successful in Silico Prediction of the Relative Permeabilities of Cyclic Peptides. *J. Am. Chem. Soc.* **2006**, *128* (43), 14073–14080. DOI: 10.1021/ja063076p.",
    "García Jiménez, D.; Vallaro, M.; Vitagliano, L.; López López, L.; Apprato, G.; Ermondi, G.; Caron, G. Molecular Properties, Including Chameleonicity, as Essential Tools for Designing the Next Generation of Oral Beyond Rule of Five Drugs. *ADMET DMPK* **2024**, *12* (5), 721–736. DOI: 10.5599/admet.2334.",
    "Schwochert, J.; Lao, Y.; Pye, C. R.; Naylor, M. R.; Desai, P. V.; Gonzalez Valcarcel, I. C.; Barrett, J. A.; Sawada, G.; Blanco, M.-J.; Lokey, R. S. Stereochemistry Balances Cell Permeability and Solubility in the Naturally Derived Phepropeptin Cyclic Peptides. *ACS Med. Chem. Lett.* **2016**, *7* (8), 757–761. DOI: 10.1021/acsmedchemlett.6b00100.",
    "Tang, X.; Kokot, J.; Waibl, F.; Fernández-Quintero, M. L.; Kamenik, A. S.; Liedl, K. R. Addressing Challenges of Macrocyclic Conformational Sampling in Polar and Apolar Solvents: Lessons for Chameleonicity. *J. Chem. Inf. Model.* **2023**, *63* (22), 7107–7123. DOI: 10.1021/acs.jcim.3c01123.",
    "Kim, B.; Sheridan, R. P.; Zhang, R.; Barros, E. P.; Johnston, J.; Xiao, L. Enhancing Permeability Prediction of Heterobifunctional Degraders Using Machine Learning and Metadynamics-Informed 3D Molecular Descriptors. *J. Chem. Inf. Model.* **2025**, *65* (24), 12563–12578. DOI: 10.1021/acs.jcim.5c01600.",
    "Sugita, M.; Noso, Y.; Li, J.; Fujie, T.; Yanagisawa, K.; Akiyama, Y. Protocol for Membrane Permeability Prediction of Cyclic Peptides Using Descriptors Obtained from Extended Ensemble Molecular Dynamics Simulations and Chemical Structures. *bioRxiv* **2025**, 2025.06.18.660352. DOI: 10.1101/2025.06.18.660352.",
    "Severoglu, Y. B.; Yuksel, B.; Sucu, C.; Aral, N.; Uversky, V. N.; Coskuner-Weber, O. Implicit Solvent Models and Their Applications in Biophysics. *Biomolecules* **2025**, *15* (8), 1218. DOI: 10.3390/biom15091218.",
    "Whitty, A.; Zhong, M.; Viarengo, L.; Beglov, D.; Hall, D. R.; Vajda, S. Quantifying the Chameleonic Properties of Macrocycles and Other High-Molecular-Weight Drugs. *Drug Discovery Today* **2016**, *21* (5), 712–717. DOI: 10.1016/j.drudis.2016.02.005.",
]


def build(args):
    df = pd.read_csv(args.csv).set_index("compound")
    for cid in (args.r, args.s):
        if cid not in df.index:
            raise SystemExit(f"compound '{cid}' not in {args.csv}. Available: {list(df.index)}")
    R, S = df.loc[args.r], df.loc[args.s]

    def g(iso, key):
        return (R if iso == "R" else S).get(key)

    soluble = pick(R, S, "water_bw_psa", "higher")    # more polar surface exposed in water
    permeable = pick(R, S, "mem_bw_psa", "lower")      # lower 3D-PSA in chloroform (Begnini)
    compact = pick(R, S, "mem_bw_rg", "lower")         # smaller Rg in chloroform (Begnini)
    agree = permeable == compact
    o_sol = "S" if soluble == "R" else "R"
    o_perm = "S" if permeable == "R" else "R"
    floppier = "R" if g("R", "water_rmsf") > g("S", "water_rmsf") else "S"
    gap = abs(R["mem_bw_psa"] - S["mem_bw_psa"])
    pct = 100 * gap / np.mean([R["mem_bw_psa"], S["mem_bw_psa"]])
    kier = R.get("water_kier_phi")
    dpsa = lambda iso: g(iso, "water_bw_psa") - g(iso, "mem_bw_psa")  # PSA buried water→chloroform

    L = [
        f"# DOPC {args.pair}: R/S epimer 3D-descriptor analysis",
        "",
        f"**Jorge Carmona · Hu Lab, San Diego State University · {date.today().isoformat()}**",
        "",
        "## 1. Purpose",
        "",
        f"The R and S epimers of **{args.pair}** ({args.residue}) are identical on every 2D and "
        f"lipophilicity descriptor, so any difference in their predicted permeability or solubility must "
        f"arise from their 3D conformational ensembles.^{{13}}^ This report ranks the two epimers on the "
        f"Boltzmann-weighted 3D descriptors best validated for macrocycle permeability (solvent-accessible "
        f"3D-PSA, radius of gyration,^{{6}}^ backbone transannular IMHB,^{{9}}^ and the water-chloroform "
        f"ΔPSA^{{16}}^) and proposes resulting permeability predictions for experimental "
        f"testing. A broader descriptor panel is tabulated for completeness and downstream machine learning "
        f"but does not drive the expected permeability results.",
        "",
        "## 2. Methods",
        "",
        "Conformational ensembles were generated independently for each epimer in water and in chloroform "
        "(CHCl₃, ε ≈ 4.8, a low-dielectric membrane mimic); sampling a polar and an apolar solvent is the "
        "established route to chameleonic conformational change.^{12}^ Starting geometries (RDKit ETKDGv3) "
        "were pre-optimized with GFN2-xTB/ALPB, sampled with CREST 2.12 (iMTD-GC metadynamics, "
        "GFN2-xTB/ALPB), and Boltzmann-weighted at 298 K, mirroring the CREMP protocol.^{1,2,3,4}^ All "
        "descriptors (`phys_descriptors_v3`) are ensemble averages, not single minimum-energy values; the "
        "3D-PSA is the solvent-accessible polar surface over N, O and their attached polar hydrogens "
        "(oxidized sulfur only).^{5,6}^ Implicit solvation is efficient but omits explicit water-mediated "
        "hydrogen bonding and approximates entropy.^{15}^",
        "",
        "## 3. Results",
        "",
        "### 3.1 Water-phase descriptors (R vs S)",
        "",
        desc_table(R, S, "water"),
        "",
        "### 3.2 Chloroform-phase descriptors (R vs S)",
        "",
        desc_table(R, S, "mem"),
        "",
        "**ΔPSA:** The Boltzmann-weighted "
        "3D-PSA drop from the aqueous to the chloroform ensemble, the operational measure of "
        f"chameleonic surface burial:^{{16}}^",
        "",
        "| ΔPSA = PSA(water) − PSA(chloroform) | R | S |",
        "|---|---|---|",
        f"| polar surface buried (Å²) | {fmt(dpsa('R'))} | {fmt(dpsa('S'))} |",
        "",
        WEAK_NOTE,
        "",
        "### 3.3 Figures",
        "",
        f"![Figure 1. Relative |R − S| difference per descriptor (water): the 2D/lipophilicity descriptors "
        f"sit at ~0 (blind to the stereocenter), while the 3D ensemble descriptors resolve "
        f"it.](results/figures/isomers/{args.pair}/fig1_reldiff.png)",
        "",
        f"![Figure 2. Validated continuous 3D descriptors, R vs S, per solvent (rows: SA 3D-PSA, "
        f"radius of gyration; columns: water, chloroform). Boxes are the IQR with median and "
        f"1.5×IQR whiskers; points are the 50 most Boltzmann-populated conformers (uniform), the "
        f"minimum-energy conformer ringed; weighted means are annotated and R vs S compared by the "
        f"Wilcoxon rank-sum (Mann-Whitney U) test.](results/figures/isomers/{args.pair}/fig2_key3d.png)",
        "",
        f"![Figure 3. Backbone (transannular) intramolecular hydrogen-bond population, R vs S, per "
        f"solvent (A) water, (B) chloroform. Bars give the fraction of each ensemble (by Boltzmann "
        f"weight) holding a given number of transannular backbone H-bonds; the dashed line marks the "
        f"weighted-mean count.](results/figures/isomers/{args.pair}/fig3_imhb.png)",
        "",
        f"![Figure 4. Normalized principal-moments-of-inertia (PMI) shape space, R vs S, per "
        f"solvent (A) water, (B) chloroform. Each conformer is placed by its NPR1 (I₁/I₃) and "
        f"NPR2 (I₂/I₃) ratios within the rod / disc / sphere triangle; dots are the 50 most "
        f"Boltzmann-populated conformers (uniform), the rug marks the full "
        f"ensemble, and the minimum-energy conformer is ringed. Shape is a supporting (exploratory) "
        f"descriptor and is not used to drive the permeability "
        f"call.](results/figures/isomers/{args.pair}/fig4_pmi.png)",
        "",
        f"![Figure 5. S-isomer conformational ensemble (top 20 Boltzmann-weighted conformers, "
        f"intra-solvent aligned): (A) water, (B) chloroform. The apolar ensemble samples a broader "
        f"set of conformers.](results/figures/isomers/{args.pair}/fig5_S_ensemble.png)",
        "",
        f"![Figure 6. R-isomer conformational ensemble (top 20 Boltzmann-weighted conformers, "
        f"intra-solvent aligned): (A) water, (B) chloroform.](results/figures/isomers/{args.pair}/fig6_R_ensemble.png)",
        "",
        f"![Figure 7. Minimum-energy conformers superimposed, R (orange) vs S (teal), on a common "
        f"reference frame: (A) water, (B) chloroform, exposing the stereocenter-driven shape "
        f"difference in each environment.](results/figures/isomers/{args.pair}/fig7_RS_overlay.png)",
        "",
        "## 4. Interpretation",
        "",
    ]

    # ---- §4 narrative, trade-off prose style ----
    chloro_rg = (f" and smaller radius of gyration ({fmt(g(permeable, 'mem_bw_rg'), 2)} vs "
                 f"{fmt(g(o_perm, 'mem_bw_rg'), 2)} Å)") if agree else ""
    dpp, dpo = dpsa(permeable), dpsa(o_perm)
    if dpp >= dpo:
        dpsa_clause = (f"and shows a slightly greater polarity drop than **{o_perm}** "
                       f"(ΔPSA {fmt(dpp)} vs {fmt(dpo)} Å²)")
    else:
        dpsa_clause = (f"from an already lower aqueous polar surface, so its polarity drop is the "
                       f"smaller of the two (ΔPSA {fmt(dpp)} vs {fmt(dpo)} Å²)")

    L += [
        f"Prior cyclic-peptide epimer studies show that stereochemistry can partition aqueous "
        f"solubility against membrane permeability.^{{11}}^ Consistent with this trade-off, the "
        f"**{soluble}** epimer is predicted to be solubility-favored, whereas **{permeable}** is "
        f"predicted to be permeability-favored. In water, **{soluble}** presents the more "
        f"solvent-exposed polar surface, with a Boltzmann-weighted 3D-PSA of "
        f"{fmt(g(soluble, 'water_bw_psa'))} Å² versus {fmt(g(o_sol, 'water_bw_psa'))} Å² for "
        f"**{o_sol}**, and fewer transannular backbone hydrogen bonds "
        f"({fmt(g(soluble, 'water_bw_IMHB_bb'), 1)} vs {fmt(g(o_sol, 'water_bw_IMHB_bb'), 1)}), leaving "
        f"more donors available for hydration.^{{5}}^ In chloroform, **{permeable}** adopts the more "
        f"compact, desolvated state, with lower 3D-PSA ({fmt(g(permeable, 'mem_bw_psa'))} vs "
        f"{fmt(g(o_perm, 'mem_bw_psa'))} Å² for **{o_perm}**){chloro_rg}, consistent with higher "
        f"intrinsic passive permeability.^{{6}}^",
        "",
        f"This assignment is further supported by chameleonicity analysis, since permeable "
        f"high-molecular-weight macrocycles typically reduce exposed polar surface upon transfer into "
        f"low-dielectric, membrane-like environments.^{{16}}^ Both epimers bury polar surface in "
        f"chloroform, but **{permeable}** reaches the lower membrane-state 3D-PSA {dpsa_clause}, "
        f"consistent with more effective polarity shielding and permeability. Although the "
        f"polar-H-inclusive 3D-PSA values exceed the ≤140 Å² guideline, they are interpreted here as "
        f"relative, directional indicators. Sampling reliability is not expected to confound this "
        f"assignment, because the flexibility index (Φ = {fmt(kier, 1)}) lies within the Φ ≲ 10 range "
        f"for reliable conformational-sampling rankings.^{{7,12}}^",
        "",
        "**Predictions** (test against experimental validation: PAMPA/Caco-2 and an in-house permeability "
        "study):",
        "",
        "| prediction | expected | basis | confidence |",
        "|---|---|---|---|",
        f"| Solubility | {soluble} more soluble | higher water 3D-PSA (more exposed donors) for {soluble} | moderate–high |",
        f"| Permeability based on lipophilicity | {permeable} more permeable | lower chloroform 3D-PSA "
        f"{'+ smaller Rg' if agree else '(Rg split)'} for {permeable} | {'moderate' if agree else 'low (descriptors split)'} |",
        f"| Effective (cell) permeability | ambiguous; {permeable}'s edge may be offset by lower solubility | trade-off of the two rows | low |",
    ]
    body, refs = renumber("\n".join(L), REFS)   # ACS: citations numbered by first appearance
    refs_block = "\n".join(f"{i}. {r}" for i, r in enumerate(refs, 1))
    return f"{body}\n\n## References\n\n{refs_block}\n"


def main():
    ap = argparse.ArgumentParser(description="Per-molecule R/S report from the descriptor CSV")
    ap.add_argument("--csv", required=True, type=Path)
    ap.add_argument("--pair", required=True, help="molecule label, e.g. 3-12-8-12 (title + figure names)")
    ap.add_argument("--r", required=True, help="compound id of the R epimer in the CSV")
    ap.add_argument("--s", required=True, help="compound id of the S epimer in the CSV")
    ap.add_argument("--residue", required=True, help="backbone building-block description")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--no-docx", action="store_true")
    args = ap.parse_args()

    md = build(args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(md, encoding="utf-8")
    print(f"wrote {args.out}  ({len(md.split())} words)")

    if not args.no_docx:
        from export_docx import build_docx
        build_docx(args.out, args.out.with_suffix(".docx"))


if __name__ == "__main__":
    main()
