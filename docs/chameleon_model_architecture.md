# CyclicPermeabilityModel — Architecture Specification

**Status: Design phase — not yet implemented**
*Last updated: 2026-03-23*

---

## Scientific Motivation

Membrane permeability of cyclic peptides is governed by three entangled phenomena:

1. **Chirality and residue-level sequence patterning** — stereochemistry (d/l residues, N-methylation, β-residues) and positional motifs drive which conformational states are accessible. This is the PI's core hypothesis.

2. **Permeability-relevant conformational transitions** — specifically whether chameleonic switching (aqueous-extended ↔ membrane-collapsed) is energetically accessible or forbidden. This is the ΔPSA / chameleonicity direction from January 2026.

3. **Membrane environment** — the membrane being crossed matters. Assay type, lipid composition, pH, and cell context all modulate apparent permeability independently of the molecule's intrinsic properties.

The model is designed to disentangle these three contributions and quantify each one's role.

---

## Module Layout

```
CyclicPermeabilityModel
├── SequenceEncoder
├── StaticDescriptorEncoder
├── DynamicEnsembleEncoder
├── AssayContextEncoder
├── ModalityFusion
├── SharedTrunk
├── RegressionHead
├── ClassificationHead
└── MechanisticHeads
    ├── DeltaPSAHead
    ├── IMHBHead
    └── StatePopulationHead
```

---

## Encoder Descriptions

### SequenceEncoder
**Captures:** residue identity, chirality, positional effects, motif logic

- Input: HELM notation or residue sequence (e.g. `[meA].[dL].[bHph].[P].[dL].[F]`)
- Requires custom tokenizer — non-standard residues (meA, bHph, dP, Bn_Gly, etc.) are not covered by standard protein language models
- Encodes: d/l chirality per position, N-methylation, β-residue type, ring closure position
- This is where stereochemistry and residue-level sequence combinations are modeled

**Featurization pipeline (pre-model):**
- HELM parser → residue vocabulary
- Per-residue features: amino acid identity, chirality flag, backbone modification (NMe, β, γ), position index, ring size

---

### StaticDescriptorEncoder
**Captures:** molecular-level 2D properties — size, polarity, lipophilicity, H-bond capacity

- Input: RDKit 2D descriptors (MolWt, MolLogP, TPSA, NumHDonors, NumHAcceptors, RingCount, FractionCSP3, etc.)
- These are conformer-independent and always available
- Serves as the baseline modality — every compound has this

---

### DynamicEnsembleEncoder
**Captures:** solvent-dependent conformational switching, compact vs exposed states, IMHB behavior, transition accessibility

- Input: conformational ensemble descriptors computed from CREMP (CHCl₃) or OpenMM GBSA-OBC (dual-solvent)
- Key features: `delta_psa3d`, `psa3d_std`, `delta_hb`, `delta_Rg`, `bw_psa3d`, `norm_delta_psa`
- Represents the January 2026 dynamic/chameleonicity direction
- **Note:** if input is raw 3D conformer geometries → requires GNN or SE(3)-equivariant transformer (future phase)
- **Current phase:** scalar descriptors only (same interface as StaticDescriptorEncoder but physically distinct features)

**Featurization pipeline (pre-model):**
- CREMP pickle → `cremp_deltapsa.py` → scalar features (current)
- OpenMM GBSA-OBC dual-solvent → Tier-2 pipeline (future)

---

### AssayContextEncoder
**Captures:** membrane/cell context, assay conditions, pH, environment dependence

- Input: structured metadata per experimental measurement
- Key fields: assay type (PAMPA / Caco-2 / MDCK / RRCK), membrane composition, pH, temperature, detection limit, source lab/protocol
- Operationalizes the insight that the membrane being crossed matters independently of the molecule
- **Risk:** on heterogeneous datasets this encoder may learn protocol identity rather than biology — only valid on source-stratified or well-characterized data

**Featurization pipeline (pre-model):**
- Structured metadata schema (to be designed — see data schema doc)
- One-hot or embedding for assay type
- Continuous values for pH, temperature, lipid composition

---

## Fusion and Trunk

### ModalityFusion
- Combines outputs of all four encoders
- Strategy TBD: concatenation, cross-attention, or gated fusion
- Must handle missing modalities gracefully (masking) — not all compounds will have all four modalities populated

### SharedTrunk
- 2–3 layer MLP operating on the fused representation
- Shared across all heads — forces the model to learn a single unified latent representation of permeability

---

## Output Heads

### RegressionHead
- Predicts continuous PAMPA permeability (log scale)
- Loss: MSE or Huber

### ClassificationHead
- Predicts binary permeable/impermeable
- Threshold: PAMPA ≥ −6.0 (standard CycPeptMPDB convention)
- Loss: BCE with class weighting

### MechanisticHeads (auxiliary — physics-informed regularization)

#### DeltaPSAHead
- Predicts ΔPSA from sequence/structure alone
- Supervision: CREMP-computed `delta_psa3d` labels
- Scientific purpose: forces the model to learn conformational switching as an intermediate representation
- **Gate:** only train once CREMP benchmark validates ΔPSA label quality

#### IMHBHead
- Predicts intramolecular H-bond count difference (membrane vs aqueous conformer)
- Supervision: `delta_hb` from conformer pipeline

#### StatePopulationHead
- Predicts Boltzmann population of lowest-energy conformer (`pop_lowest_pct`)
- Supervision: CREMP `poplowestpct`
- Scientific purpose: compounds with high population in one dominant state are likely rigid; flat distributions indicate conformational flexibility
- **Gate:** only available for CREMP overlap compounds

---

## Phased Implementation Plan

| Phase | Encoders active | Data required | Gate condition |
|-------|----------------|---------------|----------------|
| 1 | StaticDescriptor + DynamicEnsemble (scalars) | Full 7k feature_matrix.csv | Baseline — start here |
| 2 | + SequenceEncoder | Full 7k + custom HELM tokenizer | Custom residue vocabulary built |
| 3 | + MechanisticHeads | CREMP overlap subset | CREMP benchmark validates ΔPSA labels |
| 4 | + AssayContextEncoder | Source-stratified or single-protocol data | Furukawa-only or controlled experimental set |

---

## Ablation Plan

The following ablations are scientifically required (not just ML hygiene):

| Run | Encoders | Scientific question |
|-----|----------|-------------------|
| A | Static only | How much does 2D chemistry explain? |
| B | Dynamic only | How much does chameleonicity explain alone? |
| C | Sequence only | How much does sequence pattern explain? |
| D | Static + Dynamic | Does conformational info add to 2D? |
| E | Static + Sequence | Does sequence add to 2D? |
| F | Dynamic + Sequence | Can you predict without static descriptors? |
| G | Full (no context) | Best physics-grounded model |
| H | Full (with context) | Does assay context help or overfit? |

---

## Visualization Stack

| Tool | Use case |
|------|----------|
| Torchview | Development — catch architectural mistakes, verify tensor shapes |
| VisualTorch | PI slides / papers — one block per modality, clean publication figure |
| TensorBoard | Training monitoring — total loss, RMSE/MAE, AUROC/F1, per-modality ablation runs |

---

## Current Status (2026-03-23)

- [x] DynamicEnsembleEncoder featurization pipeline — `scripts/cremp_deltapsa.py` complete, `results/cremp_deltapsa.csv` generated (2,457 compounds)
- [x] StaticDescriptorEncoder features — available in `results/feature_matrix.csv`
- [ ] CREMP benchmark notebook — validates ΔPSA label quality (gate for Phase 3)
- [ ] HELM tokenizer / custom residue vocabulary
- [ ] AssayContext metadata schema
- [ ] Phase 1 PyTorch implementation
- [ ] New repo: `Chameleon_Model` (production implementations)

---

## Open Questions

1. DynamicEnsembleEncoder: scalar descriptors now, or invest in GNN on raw 3D conformer geometries?
2. ModalityFusion strategy: concatenation vs cross-attention (cross-attention is more expressive but needs more data)
3. What is the minimum overlap size between CREMP and CycPeptMPDB to make Phase 3 heads meaningful?
4. Can StatePopulationHead generalize beyond CREMP compounds if trained on the overlap subset?
