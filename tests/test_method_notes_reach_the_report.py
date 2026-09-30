"""Four facts the engine computed and no page showed.

Each is a statement about what a reading is *not*: two gene counts that are
not a flux, three protein readings that are not each other, an end product
the microbes do not make on their own, and the look-alike a difficult panel
has to reject. Losing any of them turns a careful reading into a loose one.
"""

from __future__ import annotations

from typing import Any

from openbiota import pdfsynbiotic as PS


def _render(*blocks: str, views: dict[str, Any]) -> str:
    import pypdf
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate

    from openbiota.pdfreport import _styles

    story: list[Any] = []
    st = _styles()
    for name in blocks:
        getattr(PS, name)(story, st, views)
    if not story:
        return ""
    path = "/tmp/_method_notes_render.pdf"
    SimpleDocTemplate(path, pagesize=A4).build(story)
    text = "".join(page.extract_text() for page in pypdf.PdfReader(path).pages)
    return " ".join(text.split())


# --- GABA ------------------------------------------------------------------


def test_the_two_gaba_counts_are_shown_without_being_divided() -> None:
    """A ratio here would read as a net balance. It is not one."""
    views = {"biotransformation": {"gaba_balance": {
        "synthesis": 33.5, "degradation": 12.0,
        "unit": "copies per 100 bacterial genomes",
        "ratio": None,
        "note": "Their difference is not a net flux and neither is their ratio.",
    }}}
    text = _render("_paired_readings_block", views=views)
    assert "33.5" in text and "12.0" in text
    assert "not a net flux" in text


def test_an_unmeasured_half_is_not_written_as_zero() -> None:
    """No supported reading is not the same as none present.

    Zero here would be a measurement nobody made.
    """
    views = {"biotransformation": {"gaba_balance": {
        "synthesis": 33.5, "degradation": None,
        "unit": "copies per 100 bacterial genomes", "note": "",
    }}}
    text = _render("_paired_readings_block", views=views)
    assert "no supported reading" in text
    assert "0.0" not in text.replace("33.5", "")


def test_no_gaba_block_without_either_count() -> None:
    views = {"biotransformation": {"gaba_balance": {
        "synthesis": None, "degradation": None, "unit": "", "note": "",
    }}}
    assert _render("_paired_readings_block", views=views) == ""


# --- the three protein concepts -------------------------------------------


def test_the_three_protein_readings_are_named_separately() -> None:
    views = {"nitrogen": {"three_different_concepts": {
        "protein_breakdown": "cutting proteins up",
        "bcaa_synthesis": "making branched-chain amino acids",
        "bcaa_fermentation": "fermenting them into branched-chain fatty acids",
        "note": "A community can be high in all three.",
    }}}
    text = _render("_three_concepts_block", views=views)
    assert "BCAA synthesis" in text, "an acronym must not be sentence-cased"
    assert "BCAA fermentation" in text
    assert "Protein breakdown" in text
    assert "high in all three" in text


def test_acronyms_survive_the_label_rule() -> None:
    assert PS._concept_label("bcaa_synthesis") == "BCAA synthesis"
    assert PS._concept_label("protein_breakdown") == "Protein breakdown"


# --- the host's own step ---------------------------------------------------


def test_the_host_step_is_credited_to_the_host() -> None:
    """PAGln is made by your liver. A stool gene count does not measure it."""
    views = {"nitrogen": {"host_steps": [{
        "microbial_product": "phenylacetate",
        "host_enzyme": "hepatic glutamine N-acyltransferase",
        "host_product": "phenylacetylglutamine (PAGln)",
        "note": "This measures only the microbial precursor.",
    }]}}
    text = _render("_host_steps_block", views=views)
    assert "phenylacetate" in text
    assert "PAGln" in text
    assert "only the microbial precursor" in text


# --- the decoys ------------------------------------------------------------


def test_each_difficult_panel_names_what_it_rejects() -> None:
    """A panel for a rare reaction is only as good as the near-miss it excludes."""
    views = {"biotransformation": {"negative_controls": {
        "polyphenol.urolithin_9_dehydroxylation": {
            "urdA": "UrdA acts on urocanate and makes imidazole propionate.",
        },
        "polyphenol.equol_daidzein_conversion": {
            "generic_reductase": "An unassigned oxidoreductase domain matches many things.",
        },
    }}}
    text = _render("_decoys_block", views=views)
    assert "urdA" in text
    assert "generic reductase" in text or "generic_reductase" in text
    assert "urolithin 9 dehydroxylation" in text
    assert "failed rather than found" in text


def test_no_decoy_block_when_nothing_needs_one() -> None:
    assert _render("_decoys_block", views={"biotransformation": {"negative_controls": {}}}) == ""


# --- the fermentation network ---------------------------------------------


def _edge(frm: str, to: str, route: str, support: str) -> dict[str, Any]:
    return {"from": frm, "from_label": frm, "to": to, "to_label": to,
            "route_id": route, "support": support,
            "style": "solid" if support == "supported" else "dashed"}


def test_a_supported_step_and_a_literature_step_are_told_apart() -> None:
    """An arrow that the genes support and one that only the literature
    permits must not read the same, or the gap gets filled in silently."""
    views = {"fermentation": {"network": {
        "edges": [_edge("substrate", "acetate", "acetate.pta_ack", "supported")],
        "unsupported_edges": [_edge("D-lactate", "butyrate",
                                    "lactate.utilisation_butyrate", "unsupported")],
        "limitations": ["Supported biochemical opportunities, not flows."],
    }}}
    text = _render("_network_block", views=views)
    assert "genes present" in text
    assert "literature only" in text
    assert "acetate.pta_ack" in text
    assert "not flows" in text


def test_the_network_never_claims_a_rate() -> None:
    views = {"fermentation": {"network": {
        "edges": [_edge("substrate", "acetate", "acetate.pta_ack", "supported")],
        "limitations": ["Supported biochemical opportunities, not flows. No edge carries a rate."],
    }}}
    assert "No edge carries a rate" in _render("_network_block", views=views)


# --- shared genes ----------------------------------------------------------


def test_the_shared_genes_are_named_so_the_cards_are_not_summed() -> None:
    """One gene serving three fibres is counted three times by a reader who adds."""
    views = {"substrates": {"shared_genes": {
        "n_genes": 24, "n_shared_genes": 9,
        "shared": {"amyA": ["carb.imo", "carb.resistant_starch", "carb.starch_general"]},
        "note": "Adding the cards together would count every shared gene once per card.",
    }}}
    text = _render("_shared_genes_block", views=views)
    assert "amyA" in text
    assert "resistant starch" in text
    assert "9" in text and "24" in text
    assert "once per card" in text


def test_no_shared_gene_block_when_nothing_is_shared() -> None:
    views = {"substrates": {"shared_genes": {"n_genes": 5, "n_shared_genes": 0, "shared": {}}}}
    assert _render("_shared_genes_block", views=views) == ""


def test_the_blocks_are_silent_on_an_empty_view() -> None:
    for name in ("_paired_readings_block", "_three_concepts_block",
                 "_host_steps_block", "_decoys_block",
                 "_network_block", "_shared_genes_block"):
        assert _render(name, views={}) == "", name
