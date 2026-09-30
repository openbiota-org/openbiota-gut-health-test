"""Report generation: terminal output, ``summary.txt``, ``results.json``.

Every qualitative label printed here is an arbitrary heuristic with no clinical
meaning, and says so wherever it appears. Every number in the text report also
appears in the JSON.
"""

from __future__ import annotations

import json
import textwrap
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from openbiota import __version__
from openbiota.logging_util import Style, human_duration
from openbiota.qc import SampleQC
from openbiota.tally import EntryResult, PanelResult, TallyResult
from openbiota.taxonomy import CommunityProfile

WIDTH: Final = 100

#: ARBITRARY HEURISTIC BANDS. Not derived from literature, not validated
#: against any cohort, and carrying no clinical meaning whatsoever. They exist
#: only to make a column of numbers scannable.
#:
#: The scale is interpretable for a single-copy gene: 100 copies per 100
#: bacterial genomes means roughly one copy per genome, i.e. essentially every
#: genome in the community carries the gene. 10 means roughly one genome in ten.
BANDS: Final = (
    (1.0, "trace"),
    (10.0, "low"),
    (50.0, "moderate"),
    (100.0, "high"),
    (float("inf"), "very high"),
)

#: Above this, a single-copy gene is being called more often than there are
#: genomes to carry it, which points at a normalisation or specificity problem
#: rather than at biology.
IMPLAUSIBLE_COPIES_PER_100: Final = 150.0

#: Mean translated identity below which accepted fragments are more likely to
#: be distant homologs than genuine carriage.
#:
#: Empirically, genes whose reference set actually covers the organisms
#: carrying them land at 87-98% mean identity with a median at or near 100%
#: (a read from a sequenced carrier matches its own reference exactly). Genes
#: with a handful of references land in the 67-73% range, which is where
#: unrelated members of the same protein fold sit. The threshold splits those
#: two regimes; it is a diagnostic trigger, not a filter, and no number is
#: changed by it.
LOW_IDENTITY_CAPTURE: Final = 75.0

#: Fragments below which the identity diagnostic is not worth raising.
LOW_IDENTITY_MIN_FRAGMENTS: Final = 20

#: Plausible rpoB yield, in fragments per million read pairs. rpoB is ~1,340 aa
#: (~4.0 kb) in a ~3.5 Mb gut genome, so it is ~0.11% of bacterial DNA. Two
#: reads per fragment, imperfect detection and a non-bacterial fraction put the
#: expectation near 1,000-3,000 per million pairs; the band below is
#: deliberately much wider than that.
RPOB_PER_MILLION_MIN: Final = 150.0
RPOB_PER_MILLION_MAX: Final = 12_000.0

STANDING_CAVEATS: Final = (
    "PART A measures gene capacity, not metabolite concentration. Percentiles and changes "
    "over time are more reliable than absolute values.",
    "PART B profile similarity is resemblance to a group-level published pattern, not a "
    "diagnosis and not a probability of disease. The underlying associations largely do not "
    "distinguish one condition from general gut disturbance; read every profile score beside "
    "the dysbiosis anchor.",
    "Measured performance, reference cohorts and limitations: docs/VALIDATION.md",
)


def band_for(copies_per_100: float | None) -> str:
    if copies_per_100 is None:
        return "indeterminate"
    if copies_per_100 <= 0:
        return "not detected"
    for threshold, label in BANDS:
        if copies_per_100 < threshold or threshold == float("inf"):
            return label
    return "very high"


def low_identity_capture(gene: EntryResult) -> bool:
    """Is this gene's signal more likely distant homologs than genuine carriage?"""
    return (
        gene.fragments >= LOW_IDENTITY_MIN_FRAGMENTS
        and gene.mean_identity is not None
        and gene.mean_identity < LOW_IDENTITY_CAPTURE
    )


@dataclass(frozen=True, slots=True)
class Diagnostic:
    level: str  # "info" | "warn" | "error"
    code: str
    message: str

    def to_json(self) -> dict[str, str]:
        return {"level": self.level, "code": self.code, "message": self.message}


def run_diagnostics(
    *,
    tally: TallyResult,
    qc: SampleQC | None,
    selected: Sequence[PanelResult],
    positive_control: str = "butyrate",
) -> list[Diagnostic]:
    """Cross-panel sanity checks, including the normalisation-artefact guard."""
    out: list[Diagnostic] = []
    normalizer = tally.normalizer

    if normalizer.indeterminate:
        out.append(
            Diagnostic(
                "error",
                "rpob-zero",
                "Zero rpoB fragments. Every normalised figure is INDETERMINATE — not zero. "
                "Nothing in this report can be compared to another sample. Check that the "
                "reference database built correctly and that the input really is a bacterial "
                "shotgun metagenome.",
            )
        )
    else:
        pairs = qc.read_pairs if qc is not None else None
        if pairs:
            per_million = normalizer.fragments / (pairs / 1e6)
            if per_million < RPOB_PER_MILLION_MIN:
                out.append(
                    Diagnostic(
                        "warn",
                        "rpob-low-yield",
                        f"rpoB yield is {per_million:,.0f} fragments per million read pairs, "
                        f"below the plausible floor of {RPOB_PER_MILLION_MIN:,.0f}. An "
                        "under-detected denominator inflates every panel at once. Possible "
                        "causes: a largely non-bacterial sample, or a truncated rpoB reference "
                        "set.",
                    )
                )
            elif per_million > RPOB_PER_MILLION_MAX:
                out.append(
                    Diagnostic(
                        "warn",
                        "rpob-high-yield",
                        f"rpoB yield is {per_million:,.0f} fragments per million read pairs, "
                        f"above the plausible ceiling of {RPOB_PER_MILLION_MAX:,.0f}. An "
                        "over-counted denominator deflates every panel at once.",
                    )
                )
            else:
                out.append(
                    Diagnostic(
                        "info",
                        "rpob-yield-ok",
                        f"rpoB yield is {per_million:,.0f} fragments per million read pairs, "
                        "within the plausible range for a bacterial shotgun metagenome.",
                    )
                )

    control = next((p for p in selected if p.panel.name == positive_control), None)
    if control is not None:
        if control.accepted_fragments == 0:
            out.append(
                Diagnostic(
                    "error",
                    "positive-control-zero",
                    f"The {positive_control!r} positive-control panel returned zero fragments. "
                    "Butyrate producers are abundant in any healthy stool metagenome, so this "
                    "means the PIPELINE IS BROKEN, not that the donor lacks them. Do not read "
                    "anything into the other panels.",
                )
            )
        elif not control.stable:
            out.append(
                Diagnostic(
                    "warn",
                    "positive-control-weak",
                    f"The {positive_control!r} positive control returned only "
                    f"{control.accepted_fragments} fragments, below its stability threshold. "
                    "Treat the whole run as suspect.",
                )
            )
        else:
            out.append(
                Diagnostic(
                    "info",
                    "positive-control-ok",
                    f"Positive control {positive_control!r}: {control.accepted_fragments:,} "
                    "fragments — the pipeline is detecting an abundant, well-behaved gene family.",
                )
            )

    # Interpretation guardrail: one panel high against normal others is mildly
    # interesting; everything high at once is an artefact.
    evaluable = [
        p for p in selected if p.copies_per_100_genomes is not None and p.stable
    ]
    elevated = [
        p
        for p in evaluable
        if p.copies_per_100_genomes is not None and p.copies_per_100_genomes >= 50.0
    ]
    if len(evaluable) >= 3 and len(elevated) >= max(3, int(0.8 * len(evaluable))):
        out.append(
            Diagnostic(
                "warn",
                "all-panels-elevated",
                f"{len(elevated)} of {len(evaluable)} statistically stable panels are in the "
                "'moderate' band or above at the same time. A high value across all panels "
                "indicates a systematic normalisation artefact — most likely an under-detected "
                "rpoB denominator — and NOT a finding. Do not interpret the individual panels.",
            )
        )
    elif len(elevated) == 1 and len(evaluable) >= 3:
        out.append(
            Diagnostic(
                "info",
                "one-panel-elevated",
                f"One stable panel ({elevated[0].panel.name}) is elevated against normal values "
                "for the others. That pattern is mildly interesting rather than an artefact — "
                "but see the standing caveats before reading anything into it.",
            )
        )

    for panel in selected:
        value = panel.copies_per_100_genomes
        ceiling = panel.panel.plausible_ceiling(IMPLAUSIBLE_COPIES_PER_100)
        if value is not None and value > ceiling:
            out.append(
                Diagnostic(
                    "warn",
                    f"implausible-{panel.panel.name}",
                    f"Panel {panel.panel.name!r} reports {value:,.1f} copies per 100 genomes, "
                    f"above the {ceiling:,.0f} this panel's genes could plausibly account for. "
                    "That is more copies than there are genomes to carry them, which points at "
                    "cross-family capture by the reference set or at a denominator problem.",
                )
            )
        if panel.detected and not panel.stable:
            out.append(
                Diagnostic(
                    "info",
                    f"unstable-{panel.panel.name}",
                    f"Panel {panel.panel.name!r}: {panel.accepted_fragments} fragments is below "
                    f"the configured minimum of {panel.panel.min_fragments_for_stability}. "
                    "Poisson noise dominates; the normalised figure is not meaningfully precise.",
                )
            )
        for gene in panel.targets:
            if gene.reference_truncated:
                out.append(
                    Diagnostic(
                        "info",
                        f"refs-truncated-{gene.key}",
                        f"{gene.key}: reference set was capped at {gene.reference_count} "
                        "sequences. Sensitivity is bounded by that cap; raise 'max_sequences' "
                        "in the panel YAML to widen it.",
                    )
                )
            if gene.reference_count < 25:
                out.append(
                    Diagnostic(
                        "info",
                        f"refs-small-{gene.key}",
                        f"{gene.key}: only {gene.reference_count} reference sequences. A low "
                        "value for this gene is weak evidence of absence.",
                    )
                )
            if low_identity_capture(gene):
                assert gene.mean_identity is not None
                in_aggregate = gene.entry_id in panel.panel.aggregate_target_ids()
                out.append(
                    Diagnostic(
                        "warn",
                        f"low-identity-{gene.key}",
                        f"{gene.key}: {gene.fragments:,} fragments accepted at only "
                        f"{gene.mean_identity:.0f}% mean translated identity, against "
                        f"{gene.reference_count} reference sequences. Genes whose reference set "
                        "covers the organisms carrying them land at 87-98%; this range is where "
                        "unrelated members of the same protein fold sit. Read this value as "
                        "distant-homolog capture rather than carriage — an upper bound whose "
                        "true value may be far lower."
                        + (
                            f" It feeds the {panel.panel.name!r} panel aggregate, so that "
                            "figure inherits the same doubt."
                            if in_aggregate
                            else ""
                        ),
                    )
                )

    out.extend(_residue_concordance_diagnostics(selected))

    if tally.unresolved_sseqids:
        out.append(
            Diagnostic(
                "warn",
                "unresolved-subjects",
                f"{tally.unresolved_sseqids:,} hit lines referenced a subject id absent from the "
                "reference index. The cached hits were produced against a different database — "
                "re-run with --force.",
            )
        )

    return out


#: Relative difference within which the two residue estimates are treated as
#: concordant. They are independent measurements of the same quantity, so
#: agreement is evidence and divergence is a warning.
RESIDUE_CONCORDANCE_TOLERANCE: Final = 0.25


def _residue_concordance_diagnostics(selected: Sequence[PanelResult]) -> list[Diagnostic]:
    """Cross-check the two halves of the residue filter against each other.

    The read-side check applies the published criterion to the small subset of
    fragments that reach the active site. The reference-side check applies a
    different criterion — is the matched reference itself residue-consistent —
    to every fragment. They share no evidence, so if both point at the same
    fraction of accepted fragments, that is genuine corroboration.
    """
    out: list[Diagnostic] = []
    for panel in selected:
        for gene in panel.targets:
            residue = gene.residue
            if residue is None or not residue.available:
                continue
            read_rate = residue.pass_rate
            reference_rate = residue.reference_consistent_fraction
            if read_rate is None or reference_rate is None:
                continue
            if residue.fragments_spanning < 20:
                out.append(
                    Diagnostic(
                        "info",
                        f"residue-read-side-thin-{gene.key}",
                        f"{gene.key}: only {residue.fragments_spanning} fragments spanned "
                        f"position {residue.canonical_position}, too few to cross-check the "
                        "reference-side estimate against. The reference-side figure stands "
                        "on its own.",
                    )
                )
                continue

            read_estimate = gene.fragments * read_rate
            reference_estimate = residue.fragments_on_residue_ok_reference
            denominator = max(read_estimate, reference_estimate, 1.0)
            divergence = abs(read_estimate - reference_estimate) / denominator

            if divergence <= RESIDUE_CONCORDANCE_TOLERANCE:
                out.append(
                    Diagnostic(
                        "info",
                        f"residue-concordant-{gene.key}",
                        f"{gene.key}: the two independent residue checks agree. Extrapolating "
                        f"the read-side pass rate ({read_rate:.0%} of "
                        f"{residue.fragments_spanning:,} fragments that reached the site) gives "
                        f"~{read_estimate:,.0f} fragments; the reference-side check gives "
                        f"{reference_estimate:,} — a {divergence:.0%} difference. The two share "
                        "no evidence, so agreement is real corroboration that roughly "
                        f"{reference_rate:.0%} of the accepted fragments are genuine.",
                    )
                )
            else:
                out.append(
                    Diagnostic(
                        "warn",
                        f"residue-discordant-{gene.key}",
                        f"{gene.key}: the two residue checks disagree by {divergence:.0%}. The "
                        f"read-side pass rate extrapolates to ~{read_estimate:,.0f} fragments "
                        f"but the reference-side check gives {reference_estimate:,}. One of the "
                        "two is unrepresentative — most likely the read-side subset, which is "
                        "small and biased toward references the reads align to well. Trust "
                        "neither figure more than the headline upper bound.",
                    )
                )
    return out


# --------------------------------------------------------------------------- #
# text rendering
# --------------------------------------------------------------------------- #


class _Writer:
    def __init__(self, style: Style) -> None:
        self.style = style
        self.lines: list[str] = []

    def add(self, text: str = "") -> None:
        self.lines.append(text)

    def rule(self, char: str = "-") -> None:
        self.add(char * WIDTH)

    def title(self, text: str) -> None:
        self.add(self.style.bold("=" * WIDTH))
        self.add(self.style.bold(f" {text}"))
        self.add(self.style.bold("=" * WIDTH))

    def section(self, text: str) -> None:
        self.add()
        self.add(self.style.bold(text.upper()))
        self.rule()

    def kv(self, key: str, value: object, indent: int = 2, width: int = 40) -> None:
        label = f"{key}:"
        pad = max(width - indent, len(label) + 1)
        self.add(f"{' ' * indent}{label:<{pad}}{value}")

    def wrap(self, text: str, indent: int = 2, bullet: str = "") -> None:
        prefix = " " * indent + bullet
        body = " ".join(str(text).split())
        self.lines.extend(
            textwrap.wrap(
                body,
                width=WIDTH,
                initial_indent=prefix,
                subsequent_indent=" " * len(prefix),
                # Keep hyphenated terms and long identifiers intact: splitting
                # "distant-homolog" or an accession across lines makes the
                # report harder to read and harder to grep.
                break_on_hyphens=False,
                break_long_words=False,
            )
            or [prefix]
        )

    def text(self) -> str:
        return "\n".join(self.lines) + "\n"


def _ordinal(value: float | None) -> str:
    """1st, 2nd, 3rd, 11th, 92nd — the suffix depends on the last two digits."""
    if value is None:
        return "—"
    n = int(round(value))
    suffix = "th" if 11 <= n % 100 <= 13 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suffix}"


def _fmt_copies(value: float | None) -> str:
    if value is None:
        return "indeterminate"
    if value == 0:
        return "0"
    if value < 0.01:
        return f"{value:.2e}"
    if value < 10:
        return f"{value:.3f}"
    return f"{value:,.1f}"


def _flags(panel: PanelResult) -> str:
    flags: list[str] = []
    if panel.indeterminate:
        flags.append("INDETERMINATE")
    elif not panel.detected:
        flags.append("none")
    elif not panel.stable:
        flags.append("UNSTABLE")
    # Propagate the low-identity warning from the genes that actually set this
    # panel's headline figure. Without this, the compact view — which is what
    # most people read — would show "ok" for a value built on distant homologs.
    aggregate_ids = set(panel.panel.aggregate_target_ids())
    if any(low_identity_capture(g) for g in panel.targets if g.entry_id in aggregate_ids):
        flags.append("LOW-ID")
    if panel.panel.upper_bound_reason:
        flags.append("upper-bound")
    if panel.panel.extension:
        flags.append("extension")
    return ",".join(flags) or "ok"


def _gene_flags(gene: EntryResult, panel: PanelResult) -> str:
    flags: list[str] = []
    if gene.fragments == 0:
        flags.append("none")
    elif gene.fragments < panel.panel.min_fragments_for_stability:
        flags.append("unstable")
    if low_identity_capture(gene):
        flags.append("LOW-ID")
    if gene.reference_count < 25:
        flags.append("few-refs")
    if gene.reference_truncated:
        flags.append("refs-capped")
    return ",".join(flags) or "ok"


def render_summary(
    *,
    sample: str,
    tally: TallyResult,
    selected: Sequence[PanelResult],
    qc: SampleQC | None,
    profile: CommunityProfile | None,
    diagnostics: Sequence[Diagnostic],
    run_meta: dict[str, Any],
    style: Style | None = None,
    similarity: Any | None = None,
) -> str:
    w = _Writer(style or Style(False))
    s = w.style

    w.title(f"GUT HEALTH METAGENOMIC SCREEN — sample {sample}   (openbiota v{__version__})")
    w.add()
    w.wrap(
        "Two kinds of measurement. PART A quantifies how abundant specific bacterial pathway "
        "genes are in this sample's DNA, normalised to bacterial genome equivalents: genetic "
        "capacity, not metabolite concentrations. PART B profiles which species are present "
        "and how closely the community resembles published disease-associated patterns: "
        "resemblance to a group-level pattern, not a diagnosis and not a probability of "
        "disease. Research use only.",
        indent=1,
    )

    # ---- run ------------------------------------------------------------- #
    w.section("run")
    for key, value in run_meta.items():
        w.kv(key, value)

    # ---- input ----------------------------------------------------------- #
    if qc is not None:
        w.section("input validation")
        w.kv("read pairs", f"{qc.read_pairs:,}" if qc.read_pairs else "not counted")
        total = qc.total_bases
        w.kv("total bases (estimated)", f"{total:,} ({total / 1e9:.2f} Gbp)" if total else "n/a")
        w.kv("paired", "yes" if qc.paired else "no (single mate)")
        w.kv(
            "read ids identical across mates",
            "yes" if qc.ids_identical_between_mates else "no",
        )
        w.kv("mate suffixes (/1, /2)", "present" if qc.mate_suffixes_present else "absent")
        if qc.record_counts_match is not None:
            w.kv("mate record counts match", "yes" if qc.record_counts_match else "NO")
        w.add()
        for mate in qc.mates:
            w.add(
                f"  {mate.label}: {mate.path.split('/')[-1]}  "
                f"len {mate.min_length}-{mate.max_length} bp (mean {mate.mean_length:.1f}), "
                f"GC {mate.gc_percent:.1f}%, mean Phred {mate.mean_quality:.1f}, "
                f"{len(mate.quality_alphabet)} distinct quality chars, "
                f"dup {mate.duplicate_fraction:.1%}"
            )
        if qc.notes:
            w.add()
            for note in qc.notes:
                w.wrap(note, indent=2, bullet="- ")

    # ---- normalisation --------------------------------------------------- #
    normalizer = tally.normalizer
    w.section("normalisation — rpoB (bacterial genome equivalents)")
    w.kv("rpoB fragments accepted", f"{normalizer.fragments:,}")
    w.kv("rpoB rejected, low identity", f"{normalizer.rejected_low_identity:,}")
    w.kv("rpoB rejected, short alignment", f"{normalizer.rejected_short_alignment:,}")
    w.kv("rpoB reference sequences", f"{normalizer.reference_count:,}")
    w.kv("rpoB mean reference length", f"{normalizer.reference_mean_length_aa:,.0f} aa")
    w.kv(
        "rpoB mean identity",
        f"{normalizer.mean_identity:.1f}%" if normalizer.mean_identity else "n/a",
    )
    w.kv("fragments with any hit", f"{tally.fragments_with_hit:,}")
    w.kv("raw DIAMOND hit lines", f"{tally.total_hit_lines:,}")
    if tally.hits_per_mate:
        w.kv("hit lines per mate", ", ".join(f"{k}={v:,}" for k, v in tally.hits_per_mate.items()))
    collapsed = tally.total_hit_lines - tally.fragments_with_hit
    w.kv("collapsed by mate pairing", f"{collapsed:,} duplicate read hits removed")
    w.add()
    if normalizer.indeterminate:
        w.wrap(
            s.red(
                "rpoB fragments are ZERO, so every normalised figure below is INDETERMINATE. "
                "That is not the same as zero and no per-million-reads substitute is offered."
            ),
            indent=2,
        )
    else:
        w.wrap(
            "copies_per_genome = (target_fragments / mean_target_ref_len_aa) / "
            "(rpoB_fragments / mean_rpoB_ref_len_aa); reported x100. Depth-independent, "
            "gene-length-corrected, and immune to human read contamination because human "
            "reads carry no bacterial rpoB.",
            indent=2,
        )

    # ---- headline -------------------------------------------------------- #
    w.section("PART A — metabolite pathway gene capacity")
    w.add(
        f"  {'metabolite':<34}{'gene':<10}{'copies/100':>12}{'band':>14}"
        f"{'frags':>9}{'decoy':>8}  flags"
    )
    w.rule()
    for panel in selected:
        metabolite = panel.panel.metabolite[:33]
        w.add(
            f"  {s.bold(f'{metabolite:<34}')}{'ALL':<10}"
            f"{_fmt_copies(panel.copies_per_100_genomes):>12}"
            f"{band_for(panel.copies_per_100_genomes):>14}"
            f"{panel.accepted_fragments:>9,}{panel.decoy_fragments:>8,}  {_flags(panel)}"
        )
        for gene in panel.targets:
            label = (gene.gene or gene.entry_id)[:9]
            w.add(
                f"  {'':<34}{label:<10}{_fmt_copies(gene.copies_per_100_genomes):>12}"
                f"{band_for(gene.copies_per_100_genomes):>14}"
                f"{gene.fragments:>9,}{'':>8}  {_gene_flags(gene, panel)}"
            )
            if gene.copies_per_100_genomes_residue_consistent is not None:
                consistent = gene.copies_per_100_genomes_residue_consistent
                frags = gene.residue.fragments_on_residue_ok_reference if gene.residue else 0
                w.add(
                    f"  {'':<34}{'  ↳ res-ok':<10}{_fmt_copies(consistent):>12}"
                    f"{band_for(consistent):>14}{frags:>9,}{'':>8}  residue-consistent subset"
                )
    w.rule()
    w.wrap(
        "'copies/100' is copies per 100 bacterial genomes. For a single-copy gene, 100 means "
        "roughly one copy per genome (essentially every genome carries it) and 10 means roughly "
        "one genome in ten. The 'ALL' row uses each panel's configured aggregate.",
        indent=2,
    )
    w.wrap(
        "Bands are ARBITRARY HEURISTICS with no clinical meaning and no literature basis.",
        indent=2,
    )
    if any(low_identity_capture(g) for p in selected for g in p.targets):
        w.wrap(
            s.yellow("LOW-ID")
            + f" marks a gene whose accepted fragments average under {LOW_IDENTITY_CAPTURE:.0f}% "
            "translated identity — distant-homolog capture rather than carriage. Its value is an "
            "upper bound whose true figure may be far lower. See DIAGNOSTICS.",
            indent=2,
        )

    # ---- community ------------------------------------------------------- #
    if profile is not None and profile.available:
        w.section("community composition — rpoB-derived (marker-limited; PART B has species)")
        w.kv("rpoB fragments classified", f"{profile.total_fragments:,}")
        w.kv("distinct reference organisms hit", f"{profile.distinct_references_hit:,}")
        w.kv("Shannon index (genus level)", f"{profile.shannon:.3f}" if profile.shannon else "n/a")
        w.kv("Simpson index", f"{profile.simpson:.4f}" if profile.simpson else "n/a")
        w.kv("Pielou evenness", f"{profile.evenness:.3f}" if profile.evenness else "n/a")
        if profile.fb_ratio is not None:
            w.kv("Firmicutes : Bacteroidetes", f"{profile.fb_ratio:.2f}")
        w.add()
        w.add("  phylum-level proportions (the reliable level):")
        for item in profile.phyla[:12]:
            bar = "#" * max(1, int(item.fraction * 50)) if item.fraction > 0.005 else ""
            w.add(f"    {item.label:<38}{item.fraction * 100:>6.2f}%  {item.fragments:>7,}  {bar}")
        w.add()
        w.add("  nearest reference genera (indicative only, NOT taxonomic assignment):")
        for item in profile.genera[:15]:
            w.add(f"    {item.label:<38}{item.fraction * 100:>6.2f}%  {item.fragments:>7,}")
        w.add()
        for caveat in profile.to_json()["caveats"]:
            w.wrap(caveat, indent=2, bullet="- ")

    # ---- panel detail ---------------------------------------------------- #
    for panel in selected:
        p = panel.panel
        w.section(f"panel: {p.name}{'  [extension]' if p.extension else ''}")
        w.wrap(p.description, indent=2)
        w.add()
        w.kv("metabolite", p.metabolite)
        if p.pathway:
            w.wrap(f"pathway: {' '.join(p.pathway.split())}", indent=2)
        w.kv("fragments accepted", f"{panel.accepted_fragments:,}")
        w.kv("fragments rejected as decoys", f"{panel.decoy_fragments:,}")
        w.kv("rejected below identity floor", f"{panel.rejected_low_identity:,}")
        w.kv("rejected on short alignment", f"{panel.rejected_short_alignment:,}")
        w.kv("rpoB fragments (denominator)", f"{normalizer.fragments:,}")
        w.kv(
            f"copies per 100 genomes ({panel.aggregate_method})",
            f"{_fmt_copies(panel.copies_per_100_genomes)}  [{band_for(panel.copies_per_100_genomes)}]",
        )
        if panel.copies_per_100_genomes_residue_consistent is not None:
            w.kv(
                "  residue-consistent subset",
                f"{_fmt_copies(panel.copies_per_100_genomes_residue_consistent)}  "
                f"[{band_for(panel.copies_per_100_genomes_residue_consistent)}]"
                + (
                    f"  ({panel.residue_consistent_fragments:,} fragments)"
                    if panel.residue_consistent_fragments is not None
                    else ""
                ),
            )
        if p.aggregate_reason:
            w.wrap(f"aggregate rationale: {' '.join(p.aggregate_reason.split())}", indent=2)
        if not panel.stable and panel.detected:
            w.add()
            w.wrap(
                s.yellow(
                    f"STATISTICALLY UNSTABLE: {panel.accepted_fragments} fragments is below the "
                    f"configured minimum of {p.min_fragments_for_stability}. Poisson noise "
                    "dominates; do not read the normalised figure as precise."
                ),
                indent=2,
            )
        if not panel.detected:
            w.add()
            w.wrap(
                "Not detected. With a small reference set or a low-abundance organism this is "
                "weak evidence of absence rather than proof of it.",
                indent=2,
            )

        w.add()
        w.add(
            f"  {'gene':<12}{'frags':>8}{'copies/100':>12}{'band':>13}"
            f"{'refs':>7}{'ref aa':>8}{'id%':>7}{'medid':>7}{'lowID':>7}{'short':>7}  flags"
        )
        for gene in panel.targets:
            w.add(
                f"  {(gene.gene or gene.entry_id):<12}{gene.fragments:>8,}"
                f"{_fmt_copies(gene.copies_per_100_genomes):>12}"
                f"{band_for(gene.copies_per_100_genomes):>13}"
                f"{gene.reference_count:>7,}{gene.reference_mean_length_aa:>8.0f}"
                f"{(f'{gene.mean_identity:.0f}' if gene.mean_identity else '-'):>7}"
                f"{(f'{gene.median_identity:.0f}' if gene.median_identity else '-'):>7}"
                f"{gene.rejected_low_identity:>7,}{gene.rejected_short_alignment:>7,}"
                f"  {_gene_flags(gene, panel)}"
            )
        if any(low_identity_capture(g) for g in panel.targets):
            w.wrap(
                "A gene flagged LOW-ID is averaging under "
                f"{LOW_IDENTITY_CAPTURE:.0f}% translated identity. A read from an organism that "
                "is actually in the reference set matches it near-exactly, which is why "
                "well-covered genes show a median identity at or near 100%. A low median instead "
                "means the reads are landing on distant relatives within the same protein fold.",
                indent=2,
            )

        if panel.decoys:
            w.add()
            w.add("  decoy competition (fragments whose best hit was a decoy of this panel):")
            for decoy in panel.decoys:
                w.add(
                    f"    {decoy.label[:60]:<62}{decoy.fragments:>9,} fragments"
                    + (f"  ({decoy.mean_identity:.0f}% mean id)" if decoy.mean_identity else "")
                )
        else:
            w.add()
            w.wrap(
                "This panel declares no decoys of its own. Its specificity comes from "
                "cross-panel competition instead: every panel is searched against one combined "
                "database, so this panel's targets compete against every other panel's targets "
                "and decoys. A zero here means no decoy is attributed to this panel, not that "
                "no competition took place.",
                indent=2,
            )

        # residue check
        for gene in panel.targets:
            if gene.residue is None:
                continue
            r = gene.residue
            w.add()
            w.add(f"  residue check — {gene.entry_id}, position {r.canonical_position}:")
            if not r.available:
                w.wrap(
                    "UNAVAILABLE: no reference could be aligned to the anchor, so no fragment "
                    "could be checked.",
                    indent=4,
                )
            else:
                w.kv("anchor", f"{r.anchor_accession} (residue = {r.canonical_residue})", indent=4)
                w.kv("accepted residues", "/".join(r.accepted_residues), indent=4)

                w.add()
                w.add("    reference side — applies to every accepted fragment:")
                w.kv(
                    "references with mapped position",
                    f"{r.references_mapped:,} of {r.references_total:,}",
                    indent=6,
                )
                w.kv(
                    "of those, carrying an accepted residue",
                    f"{r.references_residue_ok:,}"
                    + (
                        f" ({r.references_residue_ok / r.references_mapped:.0%})"
                        if r.references_mapped
                        else ""
                    ),
                    indent=6,
                )
                w.kv(
                    "fragments on a residue-consistent reference",
                    f"{r.fragments_on_residue_ok_reference:,}",
                    indent=6,
                )
                w.kv(
                    "fragments on a residue-INconsistent reference",
                    f"{r.fragments_on_residue_bad_reference:,}",
                    indent=6,
                )
                if r.reference_consistent_fraction is not None:
                    w.kv(
                        "residue-consistent fraction",
                        f"{r.reference_consistent_fraction:.1%}",
                        indent=6,
                    )
                if gene.copies_per_100_genomes_residue_consistent is not None:
                    w.kv(
                        "copies/100 genomes, consistent only",
                        f"{_fmt_copies(gene.copies_per_100_genomes_residue_consistent)}  "
                        f"[{band_for(gene.copies_per_100_genomes_residue_consistent)}]",
                        indent=6,
                    )
                w.wrap(
                    "A reference whose own active-site residue fails the criterion is very "
                    "likely a misannotated member of the homologous decoy family, so fragments "
                    "matching it are separated out here. The headline figure is the upper "
                    "bound; the residue-consistent figure is the tighter estimate. The truth "
                    "is between the two, because a genuine gene whose nearest reference happens "
                    "to be inconsistent is excluded from the tighter figure.",
                    indent=6,
                )

                w.add()
                w.add("    read side — the published filter, where the read reaches the site:")
                w.add(
                    f"      {r.fragments_checked:,} accepted, of which {r.fragments_spanning:,} "
                    f"spanned position {r.canonical_position}; {r.passed:,} passed, "
                    f"{r.failed:,} failed, {r.gapped:,} aligned a gap"
                )
                w.kv("not spanning the position", f"{r.not_spanned:,}", indent=6)
                w.kv("on a reference without a mapping", f"{r.no_mapping:,}", indent=6)
                if r.observed:
                    observed = ", ".join(f"{k}={v}" for k, v in sorted(r.observed.items()))
                    w.kv("observed residues", observed, indent=6)
                if r.pass_rate is not None:
                    w.kv("pass rate on the checked subset", f"{r.pass_rate:.1%}", indent=6)
                    adjusted = gene.fragments * r.pass_rate
                    w.wrap(
                        f"Extrapolating that pass rate to all {gene.fragments:,} accepted "
                        f"fragments would give ~{adjusted:,.0f}. That extrapolation assumes the "
                        "checked subset is representative, which is unverified — it is offered "
                        "as a sensitivity check, not as a corrected result.",
                        indent=6,
                    )
                else:
                    w.wrap(
                        "No fragment spanned the mapped position. At 100-151 bp that is the "
                        "expected outcome and is exactly why the reference-side check above "
                        "carries the weight here.",
                        indent=6,
                    )

        # organisms
        contributors = [g for g in panel.targets if g.organisms]
        if contributors:
            w.add()
            w.add("  closest reference organisms (NOT taxonomic identification):")
            for gene in contributors:
                w.add(f"    {gene.gene or gene.entry_id}:")
                attributed = sum(o.fragments for o in gene.organisms)
                for org in gene.organisms:
                    w.add(
                        f"      {org.fragments:>6,}  {org.mean_identity:>5.1f}%  "
                        f"{org.organism[:66]}"
                    )
                if attributed < gene.fragments:
                    w.add(
                        f"      {gene.fragments - attributed:>6,}         "
                        f"(remaining fragments across other references)"
                    )
            w.add()
            w.wrap(
                "These are the reference sequences the reads matched best. A read matching a "
                "Clostridium reference does not establish that Clostridium is present — only "
                "that the closest sequence in this reference set came from one.",
                indent=2,
            )

        if p.organisms:
            w.add()
            w.wrap(
                f"Literature-reported producers for this pathway: {', '.join(p.organisms)}.",
                indent=2,
            )

        if p.upper_bound_reason:
            w.add()
            w.wrap(s.yellow("UPPER BOUND: ") + " ".join(p.upper_bound_reason.split()), indent=2)
        for caveat in p.caveats:
            w.add()
            w.wrap(" ".join(caveat.split()), indent=2, bullet="- ")
        w.add()
        w.wrap(f"citation: {' '.join(p.citation.split())}", indent=2)

    # ---- PART B: species + profile similarity ----------------------------- #
    if similarity is not None:
        _render_similarity(w, similarity)

    # ---- diagnostics ----------------------------------------------------- #
    w.section("diagnostics")
    order = {"error": 0, "warn": 1, "info": 2}
    colour = {"error": s.red, "warn": s.yellow, "info": s.dim}
    for diagnostic in sorted(diagnostics, key=lambda d: order.get(d.level, 3)):
        marker = {"error": "[ERROR]", "warn": "[WARN ]", "info": "[info ]"}[diagnostic.level]
        w.wrap(colour[diagnostic.level](marker) + " " + diagnostic.message, indent=2)
    if not diagnostics:
        w.add("  none")

    # ---- caveats --------------------------------------------------------- #
    w.section("standing caveats")
    for caveat in STANDING_CAVEATS:
        w.wrap(caveat, indent=2, bullet="- ")
    w.add()
    return w.text()


def _render_similarity(w: _Writer, stage: Any) -> None:
    """PART B of the text report: species composition and profile similarity."""
    s = w.style
    taxonomy = stage.taxonomy
    anchor = stage.anchor

    # ---- species composition ------------------------------------------ #
    w.section("PART B — community composition (species-level, marker-based)")
    if taxonomy is None:
        w.wrap(
            "The taxonomic engine did not run, so species-level composition and profile "
            "similarity are unavailable. See the notes at the end of this part.",
            indent=2,
        )
    else:
        w.kv("profiler", f"{taxonomy.family} {taxonomy.profiler_version}")
        w.kv("database", taxonomy.index)
        if taxonomy.host is not None:
            h = taxonomy.host
            w.kv("host reads removed", h.describe())
        else:
            w.kv("host reads removed", "not filtered (host fraction unreported)")
        w.kv("unclassified fraction", f"{taxonomy.unknown_percent:.1f}%")
        species = {k: v for k, v in taxonomy.species.items() if v > 0}
        w.kv("species detected", f"{len(species):,}")
        phyla = sorted(taxonomy.phyla.items(), key=lambda kv: -kv[1])
        if phyla:
            w.add()
            w.add(f"  {'phylum':<28}{'% of classified':>16}")
            for name, value in phyla[:8]:
                w.add(f"  {name:<28}{value:>15.1f}%")
        top = sorted(species.items(), key=lambda kv: -kv[1])[:15]
        if top:
            w.add()
            w.add(f"  {'species (measured)':<44}{'% of classified':>16}")
            for name, value in top:
                w.add(f"  {name.replace('_', ' '):<44}{value:>15.2f}%")
        eco = stage.results[0].modules.get("ecological") if stage.results else None
        if eco is not None:
            for f in eco.features:
                if f.raw_value is not None:
                    w.kv(f.feature.name.replace("_", " "), f"{f.raw_value:.2f}")

    # ---- dysbiosis anchor ---------------------------------------------- #
    w.section("general dysbiosis anchor (GMWI2)")
    if anchor is None or not anchor.valid:
        w.wrap(
            "Not computed — the taxonomic engine did not run or used an incompatible database."
            if anchor is None
            else anchor.note,
            indent=2,
        )
    else:
        w.kv("GMWI2", f"{anchor.score:+.2f}  ({anchor.band})")
        w.kv(
            "model taxa present",
            f"{anchor.n_health_taxa_present} health-associated, "
            f"{anchor.n_disease_taxa_present} disease-associated",
        )
        w.add()
        w.wrap(anchor.caution, indent=2)

    # ---- profiles -------------------------------------------------------- #
    for result in stage.results:
        profile = result.profile
        w.section(f"profile: {profile.label}  [{profile.name} v{profile.version}]")
        if profile.status == "validation_control":
            w.wrap(
                s.bold("PIPELINE VALIDATION CONTROL — NOT A CANCER SCREEN AND NOT A RISK "
                       "ESTIMATE."),
                indent=2,
            )
            w.add()

        if result.abstention.abstained:
            w.wrap(s.bold("NO SCORE — the engine abstained."), indent=2)
            w.add()
            w.add("  Why:")
            for reason in result.abstention.triggered:
                w.wrap(reason, indent=4, bullet="- ")
            w.add("  What would change this:")
            for remedy in result.abstention.remedies:
                w.wrap(remedy, indent=4, bullet="- ")
            w.add()

        w.add(f"  {'module':<14}{'score':>8}{'percentile':>12}{'features':>12}  reference")
        for name in ("taxonomic", "functional", "ecological", "phenotype"):
            module = result.modules.get(name)
            if module is None:
                continue
            score = f"{module.score:+.2f}" if module.score is not None else "—"
            pct = _ordinal(module.percentile) if module.percentile is not None else "—"
            feats = f"{module.n_measured}/{len(module.features)}"
            share = (
                f"{module.weight_in_combined:.0%} of combined"
                if module.weight_in_combined
                else "excluded from combined"
            )
            w.add(
                f"  {name:<14}{score:>8}{pct:>12}{feats:>12}  n={module.reference_n:,}  {share}"
            )
        w.add("  " + "-" * 72)
        if result.combined_percentile is not None and not result.abstention.abstained:
            w.add(
                f"  {'combined similarity':<34}"
                f"{_ordinal(result.combined_percentile):>5} percentile"
            )
        else:
            w.add(f"  {'combined similarity':<34}    —")
        if result.abstention.abstained:
            w.add(f"  {'confidence':<34}not scored")
        else:
            w.add(f"  {'confidence':<34}{result.confidence.grade}")
            for reason in result.confidence.reasons:
                w.wrap(reason, indent=6, bullet="· ")
        if result.stratum:
            w.kv("duration stratum", result.stratum)
        elif profile.duration_dependent:
            w.wrap(
                "Illness duration not supplied. This pattern differs between short-term and "
                "long-term illness and the strata must not be averaged, so both readings "
                "apply and confidence is reduced. Supply --illness-duration-years.",
                indent=2,
            )

        for check in result.cross_engine:
            w.add()
            w.add(f"  CROSS-ENGINE CHECK  {check.verdict}")
            w.add(f"    taxonomic   {check.taxonomic_detail}")
            w.add(f"    gene panel  {check.functional_detail}")
            w.wrap(check.message, indent=4)

        w.add()
        w.add(f"  {'feature':<36}{'measured':>10}{'pct':>6}{'dir':>5}{'w':>5}{'q':>5}{'s':>5}{'c':>5}{'v':>7}")
        for name in ("taxonomic", "functional", "ecological", "phenotype"):
            module = result.modules.get(name)
            if module is None:
                continue
            for f in module.features:
                if f.raw_value is None:
                    measured = "absent"
                elif f.feature.engine == "metaphlan":
                    measured = f"{f.raw_value:.3f}%"
                else:
                    measured = f"{f.raw_value:.1f}"
                pct = f"{f.percentile:.0f}" if f.percentile is not None else "—"
                v = f"{f.v:+.2f}" if f.v is not None else "—"
                arrow = "↑" if f.feature.direction == "increased" else "↓"
                w.add(
                    f"  {f.feature.name.replace('_', ' ')[:35]:<36}{measured:>10}{pct:>6}{arrow:>5}"
                    f"{f.feature.w:>5.2f}{f.feature.q:>5.2f}{f.feature.s:>5.2f}{f.feature.c:>5.2f}"
                    f"{v:>7}"
                )
        w.add()
        w.wrap(" ".join(profile.combined_split_note.split()), indent=2)
        for caveat in profile.caveats:
            w.add()
            w.wrap(" ".join(caveat.split()), indent=2, bullet="! ")

    if stage.notes:
        w.section("profile similarity — notes")
        for note in stage.notes:
            w.wrap(note, indent=2, bullet="- ")


# --------------------------------------------------------------------------- #
# json + files
# --------------------------------------------------------------------------- #


def build_results_json(
    *,
    sample: str,
    tally: TallyResult,
    selected: Sequence[PanelResult],
    qc: SampleQC | None,
    profile: CommunityProfile | None,
    diagnostics: Sequence[Diagnostic],
    run_meta: dict[str, Any],
    reference_entries: dict[str, Any],
) -> dict[str, Any]:
    return {
        "tool": "openbiota",
        "version": __version__,
        "sample": sample,
        "measures": "genetic capacity (pathway gene abundance in community DNA)",
        "does_not_measure": [
            "metabolite concentrations",
            "disease risk",
        ],
        "validated_reference_range": False,
        "qualitative_bands_are_arbitrary": True,
        "run": run_meta,
        "input": None if qc is None else qc.to_json(),
        "normalisation": {
            **tally.normalizer.to_json(),
            "formula": (
                "copies_per_genome = (target_fragments / mean_target_ref_length_aa) / "
                "(rpoB_fragments / mean_rpoB_ref_length_aa)"
            ),
            "reported_as": "copies per 100 bacterial genomes",
        },
        "search": {
            "total_hit_lines": tally.total_hit_lines,
            "fragments_with_hit": tally.fragments_with_hit,
            "hit_lines_collapsed_by_mate_pairing": (
                tally.total_hit_lines - tally.fragments_with_hit
            ),
            "hits_per_mate": dict(tally.hits_per_mate),
            "unresolved_sseqids": tally.unresolved_sseqids,
        },
        "panels": [p.to_json() for p in selected],
        "panel_bands": {
            p.panel.name: band_for(p.copies_per_100_genomes) for p in selected
        },
        "community_profile": None if profile is None else profile.to_json(),
        "diagnostics": [d.to_json() for d in diagnostics],
        "reference_entries": reference_entries,
        "standing_caveats": list(STANDING_CAVEATS),
    }


def write_outputs(
    *,
    out_dir: Path,
    summary_text: str,
    results: dict[str, Any],
) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "summary.txt"
    json_path = out_dir / "results.json"
    summary_path.write_text(summary_text, encoding="utf-8")
    json_path.write_text(json.dumps(results, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    return summary_path, json_path


def render_compact(
    *,
    sample: str,
    tally: TallyResult,
    selected: Sequence[PanelResult],
    style: Style | None = None,
) -> str:
    """One-screen answer, for `--compact`."""
    w = _Writer(style or Style(False))
    s = w.style
    w.add(s.bold(f"OpenBiota Gut Health Test — {sample}  (genetic capacity, not metabolites)"))
    w.add(
        f"rpoB denominator: {tally.normalizer.fragments:,} fragments"
        + ("  [INDETERMINATE — zero rpoB]" if tally.normalizer.indeterminate else "")
    )
    w.add()
    w.add(f"{'metabolite':<36}{'copies/100 genomes':>20}{'band':>14}{'frags':>9}  flags")
    w.rule()
    for panel in selected:
        w.add(
            f"{panel.panel.metabolite[:35]:<36}"
            f"{_fmt_copies(panel.copies_per_100_genomes):>20}"
            f"{band_for(panel.copies_per_100_genomes):>14}"
            f"{panel.accepted_fragments:>9,}  {_flags(panel)}"
        )
        if panel.copies_per_100_genomes_residue_consistent is not None:
            w.add(
                f"{'  ↳ residue-consistent subset':<36}"
                f"{_fmt_copies(panel.copies_per_100_genomes_residue_consistent):>20}"
                f"{band_for(panel.copies_per_100_genomes_residue_consistent):>14}"
                f"{(panel.residue_consistent_fragments or 0):>9,}  tighter estimate"
            )
    w.rule()
    w.wrap("Bands are arbitrary heuristics with no clinical meaning.", indent=0, bullet="* ")
    if any(f in _flags(p) for p in selected for f in ("LOW-ID", "UNSTABLE", "INDETERMINATE")):
        w.wrap(
            "Flags: UNSTABLE = too few fragments, Poisson noise dominates. LOW-ID = built on "
            "distant homologs, upper bound only. INDETERMINATE = no rpoB denominator. "
            "Run without --compact for the detail.",
            indent=0,
            bullet="* ",
        )
    return w.text()


def write_run_log(path: Path, reporter_lines: str, extra: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        f"OpenBiota Gut Health Test (openbiota) v{__version__}",
        *[f"{k}: {v}" for k, v in extra.items()],
        "-" * WIDTH,
    ]
    path.write_text("\n".join(header) + "\n" + reporter_lines, encoding="utf-8")


def format_elapsed(seconds: float) -> str:
    return human_duration(seconds)
