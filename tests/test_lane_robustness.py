"""Lanes that fail, fail in the open - and the ways they used to fail quietly.

Each test here is a failure that reached a real run or the self-test:

* SingleM started with HMMER absent and died a minute in on a shell error;
* Kraken 2, told that plain FASTQ was gzipped, read nothing, exited 0 and
  the lane recorded an empty gut;
* sylph-tax read its taxonomy folder from a per-user config file holding
  one absolute path, which pointed into a repository that had been moved;
* a lane that is simply not part of this install printed as a stage that
  had finished, under a warning nobody could act on;
* tools the installer put inside the repository were found only when the
  program was started through one particular launcher.

No reference database is needed: the tools are stand-in scripts.
"""

from __future__ import annotations

import gzip
import json
import os
import stat
from pathlib import Path

import pytest

from openbiota import installcheck, projectpath
from openbiota.engines import kraken, singlem, sylph
from openbiota.expansion import lanes

REPO = Path(__file__).resolve().parents[1]
FASTQ = "@r1\nACGTACGTAC\n+\nJJJJJJJJJJ\n"


def _script(path: Path, body: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def _reads(folder: Path, *, gz: bool) -> tuple[Path, Path]:
    folder.mkdir(parents=True, exist_ok=True)
    out = []
    for mate in (1, 2):
        path = folder / (f"S_{mate}.fastq.gz" if gz else f"S_{mate}.fastq")
        path.write_bytes(gzip.compress(FASTQ.encode()) if gz else FASTQ.encode())
        out.append(path)
    return out[0], out[1]


# --------------------------------------------------------------------------- #
# Kraken 2: the compression flag is said only when true
# --------------------------------------------------------------------------- #


def test_gzip_is_recognised_by_content_not_by_name(tmp_path: Path):
    zipped, _ = _reads(tmp_path / "a", gz=True)
    plain, _ = _reads(tmp_path / "b", gz=False)
    misnamed = tmp_path / "plain_named.gz"
    misnamed.write_text(FASTQ)
    assert kraken._gzipped(zipped) and not kraken._gzipped(plain) and not kraken._gzipped(misnamed)  # noqa: SLF001
    assert kraken._read_length(misnamed) == 10, "and read as what it is"  # noqa: SLF001


def _fake_kraken(tmp_path: Path, monkeypatch, *, report_line: str = "") -> Path:
    argv = tmp_path / "kraken.argv"
    exe = _script(tmp_path / "bin" / "kraken2", f"""
printf '%s\\n' "$@" > "{argv}"
prev=""
for a in "$@"; do
  if [ "$prev" = "--report" ]; then printf '%s' "{report_line}" > "$a"; fi
  prev="$a"
done
echo "0 sequences (0.00 Mbp) processed in 0.001s" >&2
""")
    bracken = _script(tmp_path / "bin" / "bracken", "exit 1\n")
    monkeypatch.setattr(kraken, "_tools", lambda: (str(exe), str(bracken)))
    monkeypatch.setattr(kraken, "database_ready", lambda _refs, _panel="uhgg": True)
    monkeypatch.setattr(kraken, "kraken_version", lambda: "test")
    monkeypatch.setattr(kraken, "cache_file", lambda work, *_a, **_kw: work / "cache.json")
    return argv


@pytest.mark.parametrize("gz", [True, False])
def test_kraken_is_told_the_reads_are_gzipped_only_when_they_are(tmp_path: Path, monkeypatch, gz: bool):
    argv = _fake_kraken(tmp_path, monkeypatch)
    r1, r2 = _reads(tmp_path / "reads", gz=gz)
    with pytest.raises(RuntimeError):
        kraken.run_kraken(sample="S", r1=r1, r2=r2, refs_dir=tmp_path / "refs", work_dir=tmp_path / "work")
    args = argv.read_text().split("\n")
    assert ("--gzip-compressed" in args) is gz
    assert args.index("--paired") < args.index(str(r1)) < args.index(str(r2))


def test_kraken_reading_nothing_is_a_failure_not_an_empty_gut(tmp_path: Path, monkeypatch):
    _fake_kraken(tmp_path, monkeypatch)   # writes an empty report, as kraken2 does when it reads no sequence
    r1, r2 = _reads(tmp_path / "reads", gz=True)
    with pytest.raises(RuntimeError, match=r"kraken2 read no sequences from S_1\.fastq\.gz"):
        kraken.run_kraken(sample="S", r1=r1, r2=r2, refs_dir=tmp_path / "refs", work_dir=tmp_path / "work")
    assert not (tmp_path / "work" / "cache.json").exists(), "and nothing is cached as a result"


# --------------------------------------------------------------------------- #
# sylph-tax: told where its tables are
# --------------------------------------------------------------------------- #


def _fake_sylph(tmp_path: Path, monkeypatch, *, tax_exit: int = 0) -> tuple[Path, Path]:
    db_dir = tmp_path / "refs" / "sylph"
    (db_dir / "tax").mkdir(parents=True)
    (db_dir / sylph.DATABASE_FILE).write_bytes(b"db")
    with gzip.open(db_dir / "tax" / "gtdb_r232_metadata.tsv.gz", "wt") as fh:
        fh.write("GCF_000196555.1\td__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Actinomycetales;"
                 "f__Bifidobacteriaceae;g__Bifidobacterium;s__Bifidobacterium longum\n")
    row = ("S\tgtdb/GCF_000196555.1_genomic.fna.gz\t100.0\t100.0\t98.96\t5.0\tNA-NA\tHIGH\tNA-NA\t3\t3.6\t"
           "7221/9990\t98.96\t0\tNC_015067.1 Bifidobacterium longum")
    header = ("Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tTrue_cov\t"
              "ANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\t"
              "Naive_ANI\tkmers_reassigned\tContig_name")
    exe = _script(tmp_path / "bin" / "sylph", f"""
if [ "$1" = profile ]; then
  prev=""
  for a in "$@"; do
    if [ "$prev" = "-o" ]; then printf '%s\\n%s\\n' "{header}" "{row}" > "$a"; fi
    prev="$a"
  done
fi
""")
    argv = tmp_path / "tax.argv"
    tax = _script(tmp_path / "bin" / "sylph-tax", f"""
printf '%s\\n' "$@" > "{argv}"
if [ {tax_exit} -ne 0 ]; then echo "ERROR: Metadata file /somewhere/else/gtdb_r232_metadata.tsv.gz not found. Exiting."; fi
exit {tax_exit}
""")
    monkeypatch.setattr(sylph, "_tools", lambda: (exe, tax))
    return db_dir, argv


def test_sylph_tax_is_given_the_taxonomy_folder_and_no_global_config(tmp_path: Path, monkeypatch):
    db_dir, argv = _fake_sylph(tmp_path, monkeypatch)
    # A per-user config pointing into a repository that no longer exists:
    home = tmp_path / "home"
    (home / ".config" / "sylph-tax").mkdir(parents=True)
    (home / ".config" / "sylph-tax" / "config.json").write_text(json.dumps({"taxonomy_dir": "/moved/away/refs/sylph/tax"}))
    monkeypatch.setattr(Path, "home", lambda: home)
    r1, r2 = _reads(tmp_path / "reads", gz=True)
    profile = sylph.run_sylph(sample="S", r1=r1, r2=r2, db_dir=db_dir, work_dir=tmp_path / "work")
    args = argv.read_text().split("\n")
    assert args[:4] == ["--no-config", "--taxonomy-dir", str(db_dir / "tax"), "taxprof"]
    assert profile is not None
    assert [hit.species for hit in profile.hits] == ["Bifidobacterium longum"], "the lineage came from refs/, not the stale config"


def test_a_sylph_tax_failure_carries_what_it_said(tmp_path: Path, monkeypatch):
    db_dir, _ = _fake_sylph(tmp_path, monkeypatch, tax_exit=1)
    r1, r2 = _reads(tmp_path / "reads", gz=True)
    with pytest.raises(RuntimeError, match=r"sylph-tax taxprof failed \(exit 1\).*Metadata file .* not found"):
        sylph.run_sylph(sample="S", r1=r1, r2=r2, db_dir=db_dir, work_dir=tmp_path / "work")


def test_the_makefile_downloads_the_tables_without_writing_that_config():
    makefile = (REPO / "Makefile").read_text()
    assert "sylph-tax --no-config download --download-to refs/sylph/tax" in makefile


# --------------------------------------------------------------------------- #
# SingleM: the programs it starts by name
# --------------------------------------------------------------------------- #


def test_singlem_names_a_missing_program_before_it_starts(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(singlem, "_tool", lambda: tmp_path / "singlem")
    monkeypatch.setattr(singlem, "ensure_extracted", lambda _refs: True)
    monkeypatch.setattr(singlem.shutil, "which", lambda tool: None if tool == "hmmsearch" else f"/usr/bin/{tool}")
    started: list = []
    monkeypatch.setattr(singlem.subprocess, "run", lambda *a, **_kw: started.append(a))
    r1, r2 = _reads(tmp_path / "reads", gz=True)
    with pytest.raises(RuntimeError, match=r"SingleM cannot run: hmmsearch: not found on the PATH .*HMMER"):
        singlem.run_singlem(sample="S", r1=r1, r2=r2, refs_dir=tmp_path, work_dir=tmp_path / "work")
    assert started == [], "nothing was started"


def test_that_message_reads_as_an_instruction():
    from openbiota.logging_util import plain_reason

    message = "RuntimeError: SingleM cannot run: hmmsearch: not found on the PATH (hmmsearch is part of HMMER)"
    assert plain_reason(message) == "the program 'hmmsearch' is not installed - run ./install.sh to add it"


def test_singlems_programs_are_part_of_the_install():
    assert set(singlem.PATH_TOOLS) == {"orfm", "hmmsearch", "diamond"}
    install = (REPO / "install.sh").read_text()
    tools = next(line for line in install.splitlines() if line.startswith('CATALOGUE_TOOLS="'))
    assert " orfm " in tools and " hmmsearch " in tools
    assert "brew:hmmsearch|apt:hmmsearch|dnf:hmmsearch) echo hmmer" in install
    assert "hmmsearch) echo hmmer" in install.split("conda_name_for() {")[1].split("}")[0]
    assert "$(call need,hmmsearch,hmmer)" in (REPO / "Makefile").read_text()
    assert "diamond" in next(line for line in install.splitlines() if line.startswith('SYSTEM_TOOLS="'))


def test_the_readiness_check_counts_a_missing_singlem_program(tmp_path: Path, monkeypatch):
    (tmp_path / ".venv-mpa42" / "bin").mkdir(parents=True)
    (tmp_path / ".venv-mpa42" / "bin" / "metaphlan").write_text("")
    import shutil

    monkeypatch.setattr(shutil, "which", lambda tool: f"/usr/bin/{tool}")
    assert installcheck.catalogue_tools(tmp_path) is True
    monkeypatch.setattr(shutil, "which", lambda tool: None if tool == "hmmsearch" else f"/usr/bin/{tool}")
    assert installcheck.catalogue_tools(tmp_path) is False


def test_the_system_step_installs_the_lanes_binaries_too():
    """One sudo prompt, at the start; and a machine installed before a tool
    joined the list gets it on the next ./install.sh."""
    install = (REPO / "install.sh").read_text()
    done = install[install.index("done_system() {"):install.index("run_system() {")]
    assert "catalogue_tools_missing" in done and "SKIP_CATALOGUES" in done
    run = install[install.index("run_system() {"):install.index("# -- python")]
    assert "lanes_missing=$(catalogue_tools_missing)" in run
    assert 'for t in $lanes_missing; do' in run


def test_sudo_is_announced_only_when_it_is_about_to_ask():
    install = (REPO / "install.sh").read_text()
    notice = install[install.index("pkg_install() {"):install.index("ensure_micromamba() {")]
    guard, _, rest = notice.partition("needs administrator rights")
    assert "! sudo -n true 2>/dev/null" in guard.splitlines()[-2]
    assert notice.count("needs administrator rights") == 1 and "sudo -v" in rest


def test_a_python_that_is_merely_old_is_called_old():
    install = (REPO / "install.sh").read_text()
    start = install.index("ensure_tool() {")
    ensure = install[start:install.index("\n}\n", start)]
    assert "python 3.11 or newer is needed; this machine's python3 is" in ensure
    assert "from_pkg=0" in ensure, "apt would hand back the same python3; conda supplies a newer one"


# --------------------------------------------------------------------------- #
# the lane runner
# --------------------------------------------------------------------------- #


class _Recorder:
    def __init__(self) -> None:
        self.stages: list[str] = []
        self.warnings: list[str] = []
        self.recorded: list[str] = []

    def stage(self, name: str):
        import contextlib

        @contextlib.contextmanager
        def cm():
            self.stages.append(name)
            yield

        return cm()

    def warn(self, text: str) -> None:
        self.warnings.append(text)

    def record(self, text: str) -> None:
        self.recorded.append(text)


def _run_only(lane: str, tmp_path: Path, reporter: _Recorder) -> dict:
    return lanes.run_all(sample="S", r1=tmp_path / "r1", r2=tmp_path / "r2", refs_dir=tmp_path, out_dir=tmp_path,
                         threads=1, reporter=reporter, only=(lane,))["lane_status"][lane]


def test_the_optional_rescue_panel_is_recorded_not_announced(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(kraken, "available", lambda _refs, _panel="uhgg": False)
    reporter = _Recorder()
    status = _run_only("kraken_rescue", tmp_path, reporter)
    assert status["state"] == "not_assessed" and "panel not built" in status["reason"]
    assert reporter.stages == [], "a lane that is not installed is not a stage that ran"
    assert reporter.warnings == []
    assert any("panel not built" in line for line in reporter.recorded), "the log still has it"


def test_a_required_lane_that_is_not_installed_is_a_warning(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(singlem, "run_singlem", lambda **_kw: None)
    reporter = _Recorder()
    status = _run_only("singlem_globdb", tmp_path, reporter)
    assert status == {"state": "not_assessed", "reason": "tool or reference not installed"}
    assert len(reporter.warnings) == 1 and "not installed" in reporter.warnings[0]


def test_a_failed_lane_is_recorded_with_its_detail_and_warned_once(tmp_path: Path, monkeypatch):
    def boom(**_kw):
        raise RuntimeError("singlem pipe failed (exit 1): " + "x" * 400 + " hmmsearch: command not found")

    monkeypatch.setattr(singlem, "run_singlem", boom)
    reporter = _Recorder()
    status = _run_only("singlem_globdb", tmp_path, reporter)
    assert status["state"] == "failed"
    assert status["error_detail"].endswith("hmmsearch: command not found")
    assert len(reporter.warnings) == 1 and len(reporter.warnings[0]) < 320
    assert any("hmmsearch: command not found" in line for line in reporter.recorded), "the log has the end of it"


# --------------------------------------------------------------------------- #
# the project's tools, however the program is started
# --------------------------------------------------------------------------- #


def test_project_tool_folders_go_first_on_the_path_once(tmp_path: Path, monkeypatch):
    for folder in ("vendor/bin", "vendor/cargo/bin"):
        (tmp_path / folder).mkdir(parents=True)
    monkeypatch.setenv("PATH", "/usr/bin:/bin")
    added = projectpath.ensure(tmp_path)
    assert added == [str(tmp_path / "vendor/bin"), str(tmp_path / "vendor/cargo/bin")], "only folders that exist"
    assert os.environ["PATH"] == os.pathsep.join([*added, "/usr/bin", "/bin"])
    assert projectpath.ensure(tmp_path) == [], "and not a second time"
    assert os.environ["PATH"].count(str(tmp_path / "vendor/bin")) == 1


def test_the_same_folders_in_the_same_order_as_the_launcher():
    launcher = (REPO / "bin" / "openbiota").read_text()
    line = next(text for text in launcher.splitlines() if text.lstrip().startswith("export PATH"))
    positions = [line.index(f"$repo/{folder}") for folder in projectpath.TOOL_DIRS]
    assert positions == sorted(positions)


def test_every_entry_point_puts_them_on_the_path():
    cli_source = (REPO / "openbiota" / "cli.py").read_text()
    main = cli_source[cli_source.index("def main(argv"):]
    assert "projectpath.ensure(REPO_ROOT)" in main[:600]
    check = (REPO / "openbiota" / "installcheck.py").read_text()
    assert "projectpath.ensure(REPO_ROOT)" in check[check.index("def main("):]
