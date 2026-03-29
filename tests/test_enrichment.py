"""Tests for the HRA enrichment pipeline.

Runs against live APIs (BioModels, OLS, HGNC, UniProt, PubMed).
Use `pytest -v` to run, `pytest -k unit` for offline-only tests.
"""

import json
import pytest

from hra_model_access import api
from hra_model_access.enrichment import (
    enrich, _clean, _keyword_hits, _infer_temporal, _infer_atlas,
    ANATOMY_TERMS, CELL_TERMS,
)
from hra_model_access.cli import extract_api_fields, FIELDS


# ---------------------------------------------------------------------------
# Unit tests (no network)
# ---------------------------------------------------------------------------

class TestUnit:
    def test_clean_html(self):
        assert _clean("<p>Hello <b>world</b></p>") == "Hello world"
        assert _clean("  multiple   spaces  ") == "multiple spaces"

    def test_keyword_hits_anatomy(self):
        hits = _keyword_hits("the liver and pancreas model", ANATOMY_TERMS)
        assert "liver" in hits
        assert "pancreas" in hits
        assert "brain" not in hits

    def test_keyword_hits_cells(self):
        hits = _keyword_hits("hepatocyte insulin signaling in beta cell", CELL_TERMS)
        assert "hepatocyte" in hits
        assert "pancreatic beta cell" in hits

    def test_infer_temporal(self):
        assert _infer_temporal("phosphorylation cascade signaling") == "1s<1min"
        assert _infer_temporal("glucose metabolic pathway") == "1min<1hr"
        assert _infer_temporal("chronic disease progression over months") == "1week<1year"
        assert _infer_temporal("abstract mathematical model") == ""

    def test_infer_atlas(self):
        assert _infer_atlas(["liver"], [], "hepatic model", False) == "organ"
        assert _infer_atlas([], ["neuron"], "neuronal signaling", False) == "cell"
        assert _infer_atlas(["liver"], ["hepatocyte"], "liver model", False) == "organ"
        assert _infer_atlas([], [], "whole-body pharmacokinetic model", False) == "whole body"
        assert _infer_atlas([], [], "pathway model", True) == "cell"

    def test_fields_no_duplicates(self):
        assert len(FIELDS) == len(set(FIELDS))

    def test_extract_api_fields_pubmed(self):
        details = {
            "name": "TestModel",
            "publication": {"title": "A paper", "link": "http://identifiers.org/pubmed/12345",
                            "year": 2020, "month": 6},
            "curationStatus": "CURATED",
            "modellingApproach": {"name": "ordinary differential equation model"},
        }
        fields = extract_api_fields("BIOMD0000000001", details)
        assert fields["PubMed_ID"] == "12345"
        assert fields["Paper_Link"] == "https://pubmed.ncbi.nlm.nih.gov/12345"
        assert fields["BioModels_Status"] == "Yes - curated"
        assert "ordinary differential equation" in fields["Model_Type"]

    def test_extract_api_fields_doi(self):
        details = {
            "name": "Test",
            "publication": {"title": "P", "link": "http://doi.org/10.1234/test"},
            "curationStatus": "NON_CURATED",
            "modellingApproach": {},
        }
        fields = extract_api_fields("MODEL001", details)
        assert fields["DOI"] == "10.1234/test"
        assert fields["Paper_Link"] == "https://doi.org/10.1234/test"


# ---------------------------------------------------------------------------
# SBML parsing tests (no network — uses synthetic XML)
# ---------------------------------------------------------------------------

class TestSBMLParsing:
    SAMPLE_SBML = b"""<?xml version="1.0" encoding="UTF-8"?>
    <sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
      <model id="test_model">
        <listOfCompartments>
          <compartment id="cytosol" name="Cytosol"/>
          <compartment id="liver" name="Liver"/>
          <compartment id="c1" name="c1"/>
        </listOfCompartments>
        <listOfSpecies>
          <species id="s1" name="Insulin" compartment="liver">
            <annotation>
              <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
                       xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
                <rdf:Description>
                  <bqbiol:is>
                    <rdf:Bag>
                      <rdf:li rdf:resource="http://identifiers.org/uniprot/P01308"/>
                      <rdf:li rdf:resource="http://identifiers.org/go/GO:0005615"/>
                      <rdf:li rdf:resource="http://identifiers.org/chebi/CHEBI:5931"/>
                    </rdf:Bag>
                  </bqbiol:is>
                </rdf:Description>
              </rdf:RDF>
            </annotation>
          </species>
        </listOfSpecies>
      </model>
    </sbml>"""

    def test_parse_uris(self):
        data = api._parse_sbml(self.SAMPLE_SBML)
        assert "P01308" in data.uniprot
        assert "GO:0005615" in data.go
        assert "CHEBI:5931" in data.chebi

    def test_parse_compartments(self):
        data = api._parse_sbml(self.SAMPLE_SBML)
        assert "Liver" in data.compartments
        assert "Cytosol" in data.subcellular
        # "c1" should be skipped (len <= 2)
        assert "c1" not in data.compartments
        assert "c1" not in data.subcellular

    def test_parse_invalid_xml(self):
        data = api._parse_sbml(b"not xml at all")
        assert data.uniprot == set()
        assert data.compartments == []

    def test_parse_bto_cl_doid(self):
        xml = b"""<?xml version="1.0"?>
        <sbml xmlns="http://www.sbml.org/sbml/level2/version4">
          <model>
            <listOfSpecies>
              <species id="s1" name="X">
                <annotation>
                  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
                    <rdf:Description>
                      <bqbiol:is xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
                        <rdf:Bag>
                          <rdf:li rdf:resource="http://identifiers.org/bto/BTO:0000783"/>
                          <rdf:li rdf:resource="http://identifiers.org/cl/CL:0000169"/>
                          <rdf:li rdf:resource="http://identifiers.org/doid/DOID:9352"/>
                        </rdf:Bag>
                      </bqbiol:is>
                    </rdf:Description>
                  </rdf:RDF>
                </annotation>
              </species>
            </listOfSpecies>
          </model>
        </sbml>"""
        data = api._parse_sbml(xml)
        assert "BTO:0000783" in data.bto
        assert "CL:0000169" in data.cl
        assert "DOID:9352" in data.doid


# ---------------------------------------------------------------------------
# Integration tests (live APIs)
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestIntegration:
    """These tests hit live APIs. Run with `pytest -m integration`."""

    def test_fetch_model(self):
        details = api.fetch_model("BIOMD0000000356")
        assert details["name"]
        assert details.get("publication", {}).get("title")
        assert len(details.get("modelLevelAnnotations", [])) > 0

    def test_fetch_sbml(self):
        details = api.fetch_model("BIOMD0000000620")
        sbml = api.fetch_sbml(details)
        assert len(sbml.uniprot) > 0, "Should find UniProt IDs in SBML species"

    def test_fetch_abstract(self):
        text = api.fetch_abstract("21572040")
        assert "insulin" in text.lower()

    def test_lookup_uberon(self):
        label, uid = api.lookup_uberon("liver")
        assert uid.startswith("UBERON:")
        assert "liver" in label.lower()

    def test_lookup_cl(self):
        label, cid = api.lookup_cl("hepatocyte")
        assert cid.startswith("CL:")

    def test_lookup_hgnc(self):
        sym, hid = api.lookup_hgnc("INS")
        assert hid == "HGNC:6081"
        assert sym == "INS"

    def test_lookup_uniprot(self):
        info = api.lookup_uniprot("P01308")
        assert info.get("gene_symbol") == "INS"
        assert "HGNC" in info.get("hgnc_id", "")

    def test_enrich_produces_all_fields(self):
        details = api.fetch_model("BIOMD0000000620")
        result = enrich(details)
        expected_keys = {
            "Temporal_Scale", "Atlas_Scale", "System_Modeled",
            "Anatomical_Structures", "Uberon_IDs",
            "Cell_Types", "CL_IDs",
            "Key_Biomolecules", "HGNC_IDs",
            "Mapping_Confidence", "Mapping_Provenance",
        }
        assert set(result.keys()) == expected_keys

    def test_enrich_insulin_model(self):
        """BIOMD0000000620 (Palmer T2DM) should find genes, anatomy, cells."""
        details = api.fetch_model("BIOMD0000000620")
        result = enrich(details)
        assert result["Key_Biomolecules"], "Should find gene symbols"
        assert "INS" in result["Key_Biomolecules"]
        assert result["Mapping_Confidence"] in ("High", "Medium")
        assert result["Mapping_Provenance"]

    def test_search_biomodels(self):
        ids = api.search_biomodels("glucose", max_results=3)
        assert len(ids) >= 1
        assert all(id_.startswith(("BIOMD", "MODEL")) for id_ in ids)

    def test_full_pipeline(self):
        """End-to-end: search → fetch → enrich → valid JSON structure."""
        from hra_model_access.cli import process_model, FIELDS
        row = process_model("BIOMD0000000356", use_llm=False)
        assert all(f in row for f in FIELDS), "All FIELDS should be present"
        # Should be valid JSON
        json.dumps(row)
