"""Read preprocessing: FASTQ header provenance, fastp report parsing, QC gates (spec 4.2)."""

from __future__ import annotations

import gzip
from pathlib import Path

from openbiota.preprocess import (
    MIN_USABLE_NONHOST_PAIRS,
    HeaderProvenance,
    assemble_gates,
    parse_fastp_json,
    read_header_provenance,
)

# --------------------------------------------------------------------------- #
# header provenance
# --------------------------------------------------------------------------- #


def _write_fastq(path: Path, headers: list[str], *, gz: bool = False) -> Path:
    body = "".join(f"{h}\nACGTACGTAC\n+\nIIIIIIIIII\n" for h in headers)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body)
    return path


def test_illumina_headers_yield_instrument_run_flowcell_and_lane(tmp_path):
    headers = [
        f"@LH00586:435:23LHMJLT4:6:{i}:1000:2000 1:N:0:ACGTAC+TTGACA" for i in range(1101, 1131)
    ]
    path = _write_fastq(tmp_path / "r1.fastq.gz", headers, gz=True)
    prov = read_header_provenance(path)
    assert prov.style == "illumina"
    assert prov.instrument == "LH00586"
    assert prov.run == "435"
    assert prov.flowcells == ("23LHMJLT4",)
    assert prov.lanes == (6,)
    assert prov.indexes and prov.indexes[0].startswith("ACGTAC")
    assert prov.batch_key == "23LHMJLT4:6"
    assert prov.n_examined == 30
    assert prov.to_json()["batch_key"] == "23LHMJLT4:6"


def test_two_lanes_form_one_batch_key(tmp_path):
    headers = [f"@A00123:9:HXYZ:{lane}:1101:1000:{i} 1:N:0:1" for i, lane in enumerate([1, 2] * 10)]
    prov = read_header_provenance(_write_fastq(tmp_path / "r1.fastq", headers))
    assert prov.lanes == (1, 2)
    assert prov.batch_key == "HXYZ:1,2"


def test_sra_headers_are_recognised_without_a_batch_key(tmp_path):
    headers = [f"@SRR1234567.{i} {i} length=151" for i in range(1, 21)]
    prov = read_header_provenance(_write_fastq(tmp_path / "r1.fastq", headers))
    assert prov.style == "sra"
    assert prov.sra_accession == "SRR1234567"
    assert prov.batch_key is None


def test_unknown_headers_do_not_invent_provenance(tmp_path):
    prov = read_header_provenance(_write_fastq(tmp_path / "r1.fastq", ["@read1", "@read2"]))
    assert prov.style == "unknown"
    assert prov.batch_key is None
    assert prov.example == "@read1"


def test_batch_key_needs_a_flowcell():
    assert HeaderProvenance(style="illumina").batch_key is None


# --------------------------------------------------------------------------- #
# fastp JSON
# --------------------------------------------------------------------------- #


def _payload(*, paired: bool = True, adapter_reads: int = 4_000, insert_peak: int = 105):
    payload = {
        "summary": {
            "fastp_version": "1.3.6",
            "sequencing": "paired end (151 cycles + 151 cycles)",
            "before_filtering": {
                "total_reads": 2_000_000, "total_bases": 292_000_000,
                "q20_rate": 0.99, "q30_rate": 0.959, "gc_content": 0.4625,
                "read1_mean_length": 146, "read2_mean_length": 146,
            },
            "after_filtering": {"total_reads": 1_990_000, "total_bases": 290_000_000},
        },
        "filtering_result": {"low_quality_reads": 6_000, "too_many_N_reads": 1_000, "too_short_reads": 3_000},
        "duplication": {"rate": 0.0129},
        "adapter_cutting": {"adapter_trimmed_reads": adapter_reads, "adapter_trimmed_bases": 40_000},
        "insert_size": {"peak": insert_peak, "unknown": 480_000, "histogram": [0] * 100 + [520_000]},
        "read1_before_filtering": {"total_cycles": 151, "polyg_trimmed_reads": 12},
    }
    if paired:
        payload["read2_before_filtering"] = {"total_cycles": 151}
    else:
        payload["summary"]["before_filtering"].pop("read2_mean_length")
    return payload


def test_paired_counts_are_halved_to_pairs():
    result = parse_fastp_json(_payload(), mode="assess", json_path=Path("x.json"),
                              elapsed_s=1.0, cached=False, command=("fastp",))
    assert result.paired
    assert result.pairs_in == 1_000_000
    assert result.pairs_out == 995_000
    assert result.pass_fraction == 0.995
    assert result.q30_rate == 0.959
    assert result.duplication_rate == 0.0129
    assert result.adapter_trimmed_fraction == 4_000 / 2_000_000
    assert result.insert_peak == 105
    assert result.insert_measurable_fraction == 0.52
    assert result.polyg_trimmed_reads == 12
    assert result.total_cycles == 151


def test_single_end_has_no_insert_metrics():
    result = parse_fastp_json(_payload(paired=False), mode="assess", json_path=Path("x.json"),
                              elapsed_s=1.0, cached=False, command=("fastp",))
    assert not result.paired
    assert result.pairs_in == 2_000_000
    assert result.read2_mean_length is None
    assert result.insert_peak is None
    assert result.insert_measurable_fraction is None


def test_pre_trimmed_reads_are_recognised_and_explained():
    # Mean length under the cycle count with almost no adapter found: the
    # provider already trimmed. The note must say so rather than call the
    # library short-insert.
    result = parse_fastp_json(_payload(adapter_reads=50), mode="assess", json_path=Path("x.json"),
                              elapsed_s=1.0, cached=False, command=("fastp",))
    assert result.pre_trimmed
    assert "trimmed" in result.insert_note.lower()
    payload = result.to_json()
    assert payload["pairs_in"] == 1_000_000
    assert payload["mode"] == "assess"


# --------------------------------------------------------------------------- #
# QC gates
# --------------------------------------------------------------------------- #


def _fastp(pairs_in=8_000_000, pairs_out=7_990_000, dup=0.013, q30=0.96):
    payload = _payload()
    payload["summary"]["before_filtering"]["total_reads"] = pairs_in * 2
    payload["summary"]["after_filtering"]["total_reads"] = pairs_out * 2
    payload["summary"]["before_filtering"]["q30_rate"] = q30
    payload["duplication"]["rate"] = dup
    return parse_fastp_json(payload, mode="assess", json_path=Path("x.json"),
                            elapsed_s=1.0, cached=False, command=("fastp",))


def test_gates_pass_for_a_deep_clean_sample():
    gates = assemble_gates(
        input_pairs=8_000_000, fastp=_fastp(), host_fraction=0.002, nonhost_pairs=7_980_000,
        classified_fraction=0.62, rpob_fragments=8_000,
    )
    by_name = {g.name: g for g in gates.gates}
    assert by_name["input_pairs"].status == "pass"
    assert by_name["post_trim_pairs"].status == "pass"
    assert by_name["host_fraction"].status == "pass"
    assert by_name["nonhost_pairs"].status == "pass"
    assert by_name["duplicate_rate"].status == "pass"
    assert gates.overall == "pass"
    assert gates.usable_pairs == 7_980_000
    assert gates.depth_adequate is True


def test_zero_host_reads_in_a_deep_sample_is_attributed_to_the_provider():
    """Real stool always carries some human DNA; none at all means it was
    stripped upstream, and the gate says so rather than reporting a clean 0%."""
    gates = assemble_gates(
        input_pairs=8_000_000, fastp=_fastp(), host_fraction=0.0, nonhost_pairs=8_000_000,
        classified_fraction=0.62, rpob_fragments=8_000,
    )
    gate = {g.name: g for g in gates.gates}["host_fraction"]
    assert gate.status == "not_assessable"
    assert gate.value is None
    assert "not_assessable_upstream_removed" in gate.note
    # The probe-derived state says the same thing regardless of depth.
    probed = assemble_gates(
        input_pairs=8_000_000, fastp=_fastp(), host_fraction=None, nonhost_pairs=8_000_000,
        classified_fraction=0.62, rpob_fragments=8_000, host_status="not_assessable_upstream_removed",
    )
    assert {g.name: g for g in probed.gates}["host_fraction"].status == "not_assessable"
    assert probed.overall == "pass"
    # A shallow file with no host reads is not enough evidence to claim that.
    shallow = assemble_gates(
        input_pairs=50_000, fastp=None, host_fraction=0.0, nonhost_pairs=50_000,
        classified_fraction=None, rpob_fragments=None,
    )
    assert "provider" not in {g.name: g for g in shallow.gates}["host_fraction"].note


def test_shallow_sample_fails_the_prespecified_depth_gate():
    gates = assemble_gates(
        input_pairs=400_000, fastp=None, host_fraction=0.01,
        nonhost_pairs=MIN_USABLE_NONHOST_PAIRS - 1, classified_fraction=0.6, rpob_fragments=300,
    )
    by_name = {g.name: g for g in gates.gates}
    assert by_name["nonhost_pairs"].status == "fail"
    assert by_name["post_trim_pairs"].status == "unknown"  # fastp did not run
    assert gates.overall == "fail"
    assert gates.depth_adequate is False
    assert gates.to_json()["minimum_usable_pairs"] == MIN_USABLE_NONHOST_PAIRS


def test_high_host_fraction_and_low_classification_warn():
    gates = assemble_gates(
        input_pairs=8_000_000, fastp=_fastp(), host_fraction=0.20, nonhost_pairs=6_400_000,
        classified_fraction=0.35, rpob_fragments=6_000,
    )
    by_name = {g.name: g for g in gates.gates}
    assert by_name["host_fraction"].status in ("warn", "fail")
    assert by_name["classified_fraction"].status == "warn"
    assert gates.overall in ("warn", "fail")


def test_missing_stages_are_unknown_not_passes():
    gates = assemble_gates(
        input_pairs=None, fastp=None, host_fraction=None, nonhost_pairs=None,
        classified_fraction=None, rpob_fragments=None,
    )
    assert {g.status for g in gates.gates} <= {"unknown", "pass"}
    assert all(g.status == "unknown" for g in gates.gates if g.name in
               ("input_pairs", "post_trim_pairs", "host_fraction", "nonhost_pairs"))
    assert gates.usable_pairs is None
    assert gates.depth_adequate is None


def test_provenance_is_carried_into_the_gate_table(tmp_path):
    headers = [f"@LH00586:435:23LHMJLT4:6:1101:1000:{i} 1:N:0:1" for i in range(10)]
    prov = read_header_provenance(_write_fastq(tmp_path / "r1.fastq", headers))
    gates = assemble_gates(
        input_pairs=8_000_000, fastp=_fastp(), host_fraction=0.0, nonhost_pairs=8_000_000,
        classified_fraction=0.6, rpob_fragments=8_000, provenance=prov,
    )
    assert gates.provenance is prov
    batch = next(g for g in gates.gates if g.name == "batch")
    assert "23LHMJLT4" in batch.display
    assert gates.to_json()["provenance"]["batch_key"] == "23LHMJLT4:6"
