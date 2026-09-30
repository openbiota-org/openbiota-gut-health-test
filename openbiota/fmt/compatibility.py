"""Distance, overlap and diversity descriptors (spec §7.4).

None of these numbers receives a universal positive or negative weight. The
evidence does not support maximising or minimising donor–recipient distance
(E06), so they are displayed and made available to a trained model — never
folded into the coverage objective.
"""

from __future__ import annotations

from typing import Any

from .inputs import MaterialData


def _species_vector(data: MaterialData) -> dict[str, float]:
    return {
        o.source_id: float(o.value or 0.0)
        for o in data.observations
        if o.feature_kind == "species" and o.value
    }


def _strain_resolution_note() -> dict[str, Any]:
    """Strain comparison needs a strain tool; none is installed (spec §6.2)."""
    return {
        "status": "not_resolved",
        "reason": (
            "no pinned strain tool is installed in this environment (SameStr, inStrain or StrainGE), "
            "and marker alignments were not retained by the run"
        ),
        "consequence": (
            "shared species between donor and recipient cannot be called same-strain or "
            "different-strain, and no donor-specific marker support is claimed"
        ),
        "needed": [
            "SameStr v1.2025.111 or inStrain 1.10.0 with per-sample marker alignments",
            "adequate callable breadth and depth in both materials for each species compared",
        ],
    }


def descriptors(recipient: MaterialData, donor: MaterialData) -> dict[str, Any]:
    """Same-lane descriptors for one donor against the recipient."""
    r = _species_vector(recipient)
    d = _species_vector(donor)
    keys = sorted(set(r) | set(d))
    num = sum(abs(r.get(k, 0.0) - d.get(k, 0.0)) for k in keys)
    den = sum(r.get(k, 0.0) + d.get(k, 0.0) for k in keys)
    bray = (num / den) if den > 0 else None
    shared = sum(1 for k in keys if r.get(k, 0.0) > 0 and d.get(k, 0.0) > 0)
    union = sum(1 for k in keys if r.get(k, 0.0) > 0 or d.get(k, 0.0) > 0)
    jaccard = (shared / union) if union else None

    def eco(data: MaterialData, key: str) -> float | None:
        obs = data.observation("ecology", key)
        return obs.value if obs else None

    r_sh, d_sh = eco(recipient, "shannon_index"), eco(donor, "shannon_index")
    return {
        "feature_universe": "MetaPhlAn 3 species relative abundance, shared denominator",
        "bray_curtis_dissimilarity": None if bray is None else round(bray, 6),
        "jaccard_overlap": None if jaccard is None else round(jaccard, 6),
        "shared_species": shared,
        "donor_unique_species": sum(1 for k in keys if d.get(k, 0.0) > 0 and r.get(k, 0.0) == 0),
        "recipient_unique_species": sum(1 for k in keys if r.get(k, 0.0) > 0 and d.get(k, 0.0) == 0),
        "species_richness": {"recipient": eco(recipient, "species_richness"), "donor": eco(donor, "species_richness")},
        "shannon": {
            "recipient": r_sh,
            "donor": d_sh,
            "ratio": (round(d_sh / r_sh, 4) if (r_sh and d_sh and r_sh > 0) else None),
        },
        "aitchison_distance": {
            "status": "not_computed",
            "reason": "needs a frozen feature universe and a declared zero-handling rule; "
            "fitting it to the current candidate list would make one pair's number depend on "
            "which other donors were submitted",
        },
        "strain_resolution": _strain_resolution_note(),
        "interpretation": (
            "descriptors only. The evidence does not support maximising or minimising donor-recipient "
            "similarity in general: the direction differs by cohort and feature type (Behling 2025), "
            "so these values carry no weight in the coverage objective."
        ),
    }


def set_descriptors(recipient: MaterialData, donors: list[MaterialData]) -> dict[str, Any]:
    """Descriptors for a candidate set: overlap among members, not a pool forecast."""
    vectors = {d.material.material_id: _species_vector(d) for d in donors}
    ids = sorted(vectors)
    pairwise = []
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            va, vb = vectors[a], vectors[b]
            keys = set(va) | set(vb)
            num = sum(abs(va.get(k, 0.0) - vb.get(k, 0.0)) for k in keys)
            den = sum(va.get(k, 0.0) + vb.get(k, 0.0) for k in keys)
            shared = sum(1 for k in keys if va.get(k, 0.0) > 0 and vb.get(k, 0.0) > 0)
            pairwise.append(
                {
                    "material_a": a,
                    "material_b": b,
                    "bray_curtis": round(num / den, 6) if den > 0 else None,
                    "shared_species": shared,
                }
            )
    rv = _species_vector(recipient)
    union_species = set()
    for v in vectors.values():
        union_species |= {k for k, x in v.items() if x > 0}
    return {
        "member_pairwise": pairwise,
        "member_species_union": len(union_species),
        "species_absent_from_recipient_in_union": len(
            {k for k in union_species if rv.get(k, 0.0) == 0}
        ),
        "attribution_ambiguity": (
            "members sharing a species cannot be told apart without strain resolution, so the source "
            "of any such feature would be ambiguous after exposure"
            if any(p["shared_species"] for p in pairwise)
            else "no shared species among members in this lane"
        ),
        "caveat": (
            "a union of donor species is not a forecast of the community a physical combination would "
            "produce: viability, competition, resource sharing and priority effects are not modelled here"
        ),
    }


def virtual_mixture(donors: list[MaterialData], weights: dict[str, float] | None) -> dict[str, Any]:
    """An explicitly abstract analytical mixture (spec §9.2), never a recipe."""
    if not weights:
        return {
            "status": "not_requested",
            "note": "no abstract contribution weights were supplied",
        }
    total = sum(max(0.0, w) for w in weights.values())
    if total <= 0:
        return {"status": "invalid_weights", "note": "weights must be non-negative and sum above zero"}
    norm = {k: max(0.0, v) / total for k, v in weights.items()}
    mixed: dict[str, float] = {}
    for d in donors:
        w = norm.get(d.material.material_id, 0.0)
        for k, v in _species_vector(d).items():
            mixed[k] = mixed.get(k, 0.0) + w * v
    return {
        "status": "scenario_only",
        "weights": norm,
        "unit": "relative_percent_same_lane",
        "top_species": sorted(mixed.items(), key=lambda kv: -kv[1])[:20],
        "label": (
            "analytical contribution scenario. These weights are NOT stool mass, viable-cell "
            "quantities, capsule counts or mixing proportions, and no preparation instruction follows "
            "from them."
        ),
    }
