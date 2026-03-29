# hra-model-access

Map [BioModels](https://www.ebi.ac.uk/biomodels/) entries to [Human Reference Atlas](https://humanatlas.io/) ontologies (UBERON, CL, HGNC).

## Install

```bash
pip install -e .
```

## Usage

```bash
# Search BioModels and generate mapping
generate-hra-mapping --query glucose --max-results 10 --output glucose.json

# Specific model IDs
generate-hra-mapping --ids BIOMD0000000356 BIOMD0000000620

# Curated ODE models
generate-hra-mapping \
  --query 'curationstatus:"Manually curated" AND modelformat:"SBML"' \
  --max-results 100 --output curated.json
```

## Data sources

For each model the pipeline queries, in order:

1. **SBML species annotations** — UniProt, BTO, CL, GO, ChEBI, DOID identifiers from `<rdf:li>` elements
2. **SBML compartments** — compartment names mapped to UBERON via OLS
3. **BioModels model-level annotations** — BTO tissues, GO terms, disease and taxonomy
4. **PubMed abstract** — keyword extraction for anatomy and cell types
5. **UniProt tissue specificity** — tissue expression text from protein entries
6. **Description keywords** — model name and description scanned against anatomy/cell term tables
7. **Disease inference** — disease annotations mapped to likely organs

Ontology resolution uses the EBI [OLS](https://www.ebi.ac.uk/ols4/), [HGNC](https://rest.genenames.org/), and [UniProt](https://rest.uniprot.org/) REST APIs. All API responses are cached to disk at `~/.cache/hra_model_access/`.

## Output fields

| Field | Source |
|---|---|
| Model_ID, Model_Name, Paper_Title, Paper_Publication_Date, Paper_Link, DOI, PubMed_ID | BioModels API |
| BioModels_Status, Model_Type, Primary_Source_URL | BioModels API |
| Anatomical_Structures, Uberon_IDs | SBML + OLS UBERON |
| Cell_Types, CL_IDs | SBML + OLS CL |
| Key_Biomolecules, HGNC_IDs | SBML UniProt + HGNC |
| Temporal_Scale, Atlas_Scale | Keyword inference |
| System_Modeled | GO terms / description |
| Mapping_Confidence | High/Medium/Low based on field coverage |
| Mapping_Provenance | Per-model explanation of what sources contributed |

## Optional LLM enrichment

Install with `pip install -e ".[llm]"` and set `ANTHROPIC_API_KEY` to use Claude for filling remaining gaps. Omit `--no-llm` to enable.

## Tests

```bash
pip install -e ".[dev]"

# Unit tests (no network)
pytest tests/ -k "Unit or SBML"

# Integration tests (hits live APIs)
pytest tests/ -m integration
```
