# CyclicPermeabilityModel — Data Schema

**Status: Active design**
*Last updated: 2026-03-23*

---

## Design Decisions (locked)

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Database | SQLite | Single researcher, integrates cleanly with pandas |
| Transitions table | Placeholder — do not populate | State definition not yet pinned. Populating prematurely creates scientifically incoherent rows |
| residue_contributions | Phase 3 only | Expensive to compute; only meaningful once you know which peptides matter |
| Versioning | method_version on every computed table, never silently overwrite | Reproducibility: must always be able to trace what was computed, how, with which assumptions |

---

## Phased Schema Scope

| Phase | Tables active | What you can do |
|-------|--------------|-----------------|
| 1 | compounds, residues, static_descriptors, conformational_ensembles, environments, dynamic_features, assay_types, assay_contexts, permeability_results, orthogonal_measurements | Train Phase 1 model, full ablation |
| 2 | + conformational_states | State-level features, better mechanistic labels |
| 3 | + residue_contributions, transitions (if state definition pinned) | Residue-resolved dynamic contributions, MechanisticHeads training |

---

## Full Normalized Schema

### compounds
Central registry. One row per unique molecule (stereochemistry-preserving).

```sql
CREATE TABLE compounds (
    compound_id        INTEGER PRIMARY KEY,
    smiles_canonical   TEXT NOT NULL UNIQUE,   -- RDKit canonical, stereo-preserving
    smiles_input       TEXT,                   -- original as deposited
    inchi              TEXT UNIQUE,
    inchikey           TEXT UNIQUE,
    helm               TEXT,                   -- full HELM string
    ring_size          INTEGER,                -- number of residues in macrocycle
    mw                 REAL,
    formula            TEXT,
    source_db          TEXT,                   -- 'CycPeptMPDB','CREMP','experimental','literature'
    source_id          TEXT,                   -- original ID in source DB
    date_added         DATE,
    notes              TEXT
);
```

---

### residues
Per-position residue data. One row per residue per compound.

```sql
CREATE TABLE residues (
    residue_id         INTEGER PRIMARY KEY,
    compound_id        INTEGER NOT NULL REFERENCES compounds(compound_id),
    position           INTEGER NOT NULL,       -- 1-indexed, N→C direction
    residue_code       TEXT NOT NULL,          -- HELM monomer code e.g. 'meA','dL','bHph'
    amino_acid_parent  TEXT,                   -- canonical AA parent e.g. 'Ala','Leu'
    chirality          TEXT CHECK (chirality IN ('L','D','achiral','unknown')),
    n_methylated       INTEGER DEFAULT 0,      -- boolean (0/1)
    backbone_type      TEXT CHECK (backbone_type IN ('alpha','beta','gamma','peptoid','other')),
    sidechain_smiles   TEXT,                   -- SMILES of sidechain only
    is_proline_type    INTEGER DEFAULT 0,      -- Pro and dP constrain phi; boolean (0/1)
    UNIQUE (compound_id, position)
);
```

---

### static_descriptors
2D molecular descriptors. Conformer-independent. Always available. One row per compound.

```sql
CREATE TABLE static_descriptors (
    descriptor_id      INTEGER PRIMARY KEY,
    compound_id        INTEGER NOT NULL REFERENCES compounds(compound_id),
    method_version     TEXT NOT NULL,          -- e.g. 'rdkit_2024.03' — NEVER silently overwrite
    mw                 REAL,
    mol_logp           REAL,
    tpsa_2d            REAL,
    num_hbd            INTEGER,
    num_hba            INTEGER,
    num_rotatable      INTEGER,
    fraction_csp3      REAL,
    ring_count         INTEGER,
    num_aromatic_rings INTEGER,
    molar_refractivity REAL,
    qed                REAL,
    computed_at        TEXT,                   -- ISO timestamp
    UNIQUE (compound_id, method_version)       -- allows re-versioning without deleting old rows
);
```

---

### environments
Definitions of solvent/membrane environments used for conformer generation.
Insert once, reference forever.

```sql
CREATE TABLE environments (
    environment_id     INTEGER PRIMARY KEY,
    name               TEXT NOT NULL UNIQUE,   -- 'vacuum','CHCl3_implicit','water_implicit',
                                               -- 'DMPC_explicit','POPC_explicit'
    solvent            TEXT,
    dielectric         REAL,                   -- relative permittivity
    temperature_K      REAL,
    method             TEXT,                   -- 'ETKDGv3_MMFF94s','CREST_GFN2','OpenMM_GBSA'
    notes              TEXT
);

-- Seed data
INSERT INTO environments VALUES
    (1,'vacuum',          'none',  1.0,   298.15,'ETKDGv3_MMFF94s','RDKit vacuum conformer generation'),
    (2,'CHCl3_implicit',  'CHCl3', 4.8,   298.15,'CREST_GFN2',     'CREMP dataset — CHCl3 implicit solvent'),
    (3,'water_implicit',  'H2O',   80.0,  298.15,'OpenMM_GBSA_OBC','GBSA-OBC water phase'),
    (4,'membrane_implicit','lipid', 4.8,  490.0, 'OpenMM_GBSA_OBC','GBSA-OBC membrane phase — 490K accelerated sampling');
```

---

### conformational_ensembles
One ensemble per compound per environment per source.

```sql
CREATE TABLE conformational_ensembles (
    ensemble_id        INTEGER PRIMARY KEY,
    compound_id        INTEGER NOT NULL REFERENCES compounds(compound_id),
    environment_id     INTEGER NOT NULL REFERENCES environments(environment_id),
    source             TEXT NOT NULL,          -- 'CREMP','local_ETKDGv3','OpenMM_tier2'
    n_confs_total      INTEGER,
    n_confs_unique     INTEGER,
    ensemble_energy    REAL,                   -- kcal/mol
    pop_lowest_pct     REAL,                   -- Boltzmann pop of lowest energy state (%)
    method_version     TEXT NOT NULL,          -- e.g. 'cremp_deltapsa_v1.0'
    parameter_json     TEXT,                   -- JSON of all method parameters
    computed_at        TEXT,                   -- ISO timestamp
    UNIQUE (compound_id, environment_id, source, method_version)
);
```

---

### dynamic_features
Conformational descriptors derived from an ensemble.
One row per ensemble. Never silently overwrite — bump method_version.

```sql
CREATE TABLE dynamic_features (
    feature_id         INTEGER PRIMARY KEY,
    ensemble_id        INTEGER NOT NULL REFERENCES conformational_ensembles(ensemble_id),
    method_version     TEXT NOT NULL,          -- e.g. 'dynfeat_v1.0'
    parameter_json     TEXT,                   -- JSON: psa_probe_radius, n_confs_sampled, etc.
    -- PSA descriptors
    aq_psa3d           REAL,                   -- max-PSA conformer polar SASA (Å²)
    mem_psa3d          REAL,                   -- min-PSA conformer polar SASA (Å²)
    delta_psa3d        REAL,                   -- aq_psa3d - mem_psa3d
    psa3d_std          REAL,                   -- std across sampled ensemble
    psa3d_spread       REAL,                   -- max - min across sampled ensemble
    bw_psa3d           REAL,                   -- Boltzmann-weighted mean PSA
    norm_delta_psa     REAL,                   -- delta_psa3d / SASA_total (Yu 2026)
    -- H-bond descriptors
    aq_hb_count        INTEGER,
    mem_hb_count       INTEGER,
    delta_hb           INTEGER,                -- mem - aq (positive = more HB buried in membrane)
    -- Shape descriptors
    aq_rg              REAL,
    mem_rg             REAL,
    delta_rg           REAL,
    aq_npr1            REAL,
    mem_npr1           REAL,
    aq_npr2            REAL,
    mem_npr2           REAL,
    -- Metadata
    n_confs_sampled    INTEGER,                -- how many conformers were actually evaluated
    computed_at        TEXT,
    UNIQUE (ensemble_id, method_version)
);
```

---

### conformational_states
*(Phase 2 — do not populate until state definition is pinned)*

A state is a metastable region of conformational space.
Definition must be locked before any rows are inserted.

```sql
CREATE TABLE conformational_states (
    state_id           INTEGER PRIMARY KEY,
    ensemble_id        INTEGER NOT NULL REFERENCES conformational_ensembles(ensemble_id),
    state_label        TEXT NOT NULL,          -- e.g. 'open','closed','intermediate'
    state_definition   TEXT NOT NULL,          -- human-readable: what rule assigns a frame here
    method_version     TEXT NOT NULL,
    parameter_json     TEXT,                   -- thresholds used for state assignment
    -- Descriptor summary for this state
    mean_psa3d         REAL,
    mean_rg            REAL,
    mean_hb_count      INTEGER,
    population_pct     REAL,                   -- Boltzmann population %
    n_frames           INTEGER,
    computed_at        TEXT,
    UNIQUE (ensemble_id, state_label, method_version)
);
```

---

### transitions
*(Phase 3 — PLACEHOLDER. Do not populate until ALL of the following are defined:)*
- *What is a state? (locked definition in conformational_states)*
- *What counts as a transition? (minimum residence time, probability threshold)*
- *At what level is it reported? (per trajectory, per replicate, aggregated)*

```sql
CREATE TABLE transitions (
    transition_id      INTEGER PRIMARY KEY,
    compound_id        INTEGER NOT NULL REFERENCES compounds(compound_id),
    state_from_id      INTEGER NOT NULL REFERENCES conformational_states(state_id),
    state_to_id        INTEGER NOT NULL REFERENCES conformational_states(state_id),
    environment_id     INTEGER NOT NULL REFERENCES environments(environment_id),
    -- Energetics
    delta_g_kcal       REAL,                   -- free energy difference
    barrier_kcal       REAL,                   -- activation energy (if computed)
    switching_accessible INTEGER,              -- boolean: accessible at 310K?
    -- Statistical
    transition_prob    REAL,                   -- observed probability from trajectory
    mean_first_passage_ns REAL,                -- MFPT if computed
    -- Versioning
    method_version     TEXT NOT NULL,
    parameter_json     TEXT,
    computed_at        TEXT,
    UNIQUE (compound_id, state_from_id, state_to_id, method_version)
);
```

---

### Assay and Experimental Tables

```sql
CREATE TABLE assay_types (
    assay_type_id      INTEGER PRIMARY KEY,
    name               TEXT NOT NULL UNIQUE,   -- 'PAMPA','Caco2','MDCK','RRCK','EPSA'
    membrane_type      TEXT,                   -- 'artificial_lipid','cell_monolayer','none'
    lipid_composition  TEXT,
    direction          TEXT CHECK (direction IN ('A_to_B','B_to_A','bidirectional','unknown')),
    measures_efflux    INTEGER DEFAULT 0,
    notes              TEXT
);

CREATE TABLE assay_contexts (
    context_id         INTEGER PRIMARY KEY,
    assay_type_id      INTEGER NOT NULL REFERENCES assay_types(assay_type_id),
    lab_source         TEXT NOT NULL,          -- 'Furukawa_2016','internal_2026', etc.
    ph                 REAL,
    temperature_c      REAL,
    incubation_h       REAL,
    detection_method   TEXT,                   -- 'UV','LC-MS','fluorescence'
    detection_limit    REAL,                   -- log Papp lower bound
    buffer_composition TEXT,
    notes              TEXT
);

CREATE TABLE permeability_results (
    result_id              INTEGER PRIMARY KEY,
    compound_id            INTEGER NOT NULL REFERENCES compounds(compound_id),
    context_id             INTEGER NOT NULL REFERENCES assay_contexts(context_id),
    log_papp               REAL,              -- primary value (log scale)
    papp_raw               REAL,              -- raw Papp cm/s if available
    efflux_ratio           REAL,              -- Caco-2/MDCK efflux ratio
    permeable              INTEGER,           -- binary label (0/1)
    threshold_used         REAL,             -- log Papp threshold applied
    replicate_n            INTEGER DEFAULT 1,
    std_dev                REAL,
    is_at_detection_limit  INTEGER DEFAULT 0, -- boolean
    data_quality           TEXT CHECK (data_quality IN ('high','medium','low','flagged')),
    reference              TEXT,
    UNIQUE (compound_id, context_id)
);

CREATE TABLE orthogonal_measurements (
    measurement_id     INTEGER PRIMARY KEY,
    compound_id        INTEGER NOT NULL REFERENCES compounds(compound_id),
    measurement_type   TEXT NOT NULL CHECK (measurement_type IN (
                           'kinetic_solubility','thermodynamic_solubility',
                           'plasma_stability','microsomal_stability',
                           'aggregation_flag','pka','logd')),
    value              REAL,
    unit               TEXT,
    conditions         TEXT,
    lab_source         TEXT,
    reference          TEXT,
    measured_at        TEXT
);
```

---

### residue_contributions
*(Phase 3 only — expensive to compute)*

```sql
CREATE TABLE residue_contributions (
    contribution_id    INTEGER PRIMARY KEY,
    ensemble_id        INTEGER NOT NULL REFERENCES conformational_ensembles(ensemble_id),
    position           INTEGER NOT NULL,
    residue_code       TEXT NOT NULL,
    method_version     TEXT NOT NULL,
    -- SASA per residue in each conformational state
    aq_sasa            REAL,                   -- Å², max-PSA conformer
    mem_sasa           REAL,                   -- Å², min-PSA conformer
    delta_sasa         REAL,
    aq_burial          REAL,                   -- 0=exposed, 1=buried
    mem_burial         REAL,
    -- HB participation
    aq_hb_donor        INTEGER,               -- boolean
    aq_hb_acceptor     INTEGER,
    mem_hb_donor       INTEGER,
    mem_hb_acceptor    INTEGER,
    computed_at        TEXT,
    UNIQUE (ensemble_id, position, method_version)
);
```

---

## ML Feature Assembly View

One row per compound per assay context. Missingness flags included — do not drop NULLs, mask them during training.

```sql
CREATE VIEW ml_feature_matrix AS
SELECT
    c.compound_id,
    c.smiles_canonical,
    c.ring_size,
    c.source_db,
    -- static
    sd.mw, sd.mol_logp, sd.tpsa_2d,
    sd.num_hbd, sd.num_hba, sd.fraction_csp3, sd.qed,
    -- dynamic — CHCl3 ensemble (membrane environment, CREMP)
    df_mem.delta_psa3d          AS delta_psa3d_chcl3,
    df_mem.norm_delta_psa       AS norm_delta_psa_chcl3,
    df_mem.psa3d_std            AS psa3d_std_chcl3,
    df_mem.delta_hb             AS delta_hb_chcl3,
    df_mem.delta_rg             AS delta_rg_chcl3,
    df_mem.bw_psa3d             AS bw_psa3d_chcl3,
    -- dynamic — water ensemble (if available)
    df_aq.delta_psa3d           AS delta_psa3d_water,
    df_aq.norm_delta_psa        AS norm_delta_psa_water,
    -- missingness flags — mask these during training, never impute silently
    CASE WHEN sd.descriptor_id   IS NOT NULL THEN 1 ELSE 0 END AS has_static,
    CASE WHEN df_mem.feature_id  IS NOT NULL THEN 1 ELSE 0 END AS has_dynamic_chcl3,
    CASE WHEN df_aq.feature_id   IS NOT NULL THEN 1 ELSE 0 END AS has_dynamic_water,
    -- assay context
    at.name                     AS assay_type,
    ac.lab_source,
    ac.ph,
    ac.temperature_c,
    -- labels
    pr.log_papp,
    pr.permeable,
    pr.data_quality,
    pr.is_at_detection_limit,
    -- provenance
    sd.method_version           AS static_method_version,
    df_mem.method_version       AS dynamic_method_version
FROM compounds c
LEFT JOIN static_descriptors sd
    ON c.compound_id = sd.compound_id
LEFT JOIN conformational_ensembles ce_mem
    ON c.compound_id = ce_mem.compound_id
    AND ce_mem.environment_id = 2             -- CHCl3_implicit
LEFT JOIN dynamic_features df_mem
    ON ce_mem.ensemble_id = df_mem.ensemble_id
LEFT JOIN conformational_ensembles ce_aq
    ON c.compound_id = ce_aq.compound_id
    AND ce_aq.environment_id = 3             -- water_implicit
LEFT JOIN dynamic_features df_aq
    ON ce_aq.ensemble_id = df_aq.ensemble_id
LEFT JOIN permeability_results pr
    ON c.compound_id = pr.compound_id
LEFT JOIN assay_contexts ac
    ON pr.context_id = ac.context_id
LEFT JOIN assay_types at
    ON ac.assay_type_id = at.assay_type_id;
```

---

## Versioning Rules

**Never silently replace computed outputs.** When you recompute anything:

1. Insert new rows with a new `method_version` string
2. Keep old rows queryable
3. The ML view always uses the latest version — update the view if needed
4. `parameter_json` stores all parameters used, so any run is fully reproducible

**Version naming convention:**
```
{table_prefix}_v{major}.{minor}

Examples:
  static_descriptors:    rdkit_2024.03
  dynamic_features:      dynfeat_v1.0  →  dynfeat_v1.1 (after IMHB cutoff change)
  conformational_states: states_v0.1   (placeholder — not yet defined)
  transitions:           transitions_v0.1 (placeholder — not yet defined)
```

---

## Open Items Before Phase 2

- [ ] Pin state definition: what is a conformational state in this project?
- [ ] Pin transition definition: what counts as a transition, at what level?
- [ ] Build HELM parser and residue vocabulary for SequenceEncoder
- [ ] Decide on `parameter_json` schema (what parameters are required vs optional)
