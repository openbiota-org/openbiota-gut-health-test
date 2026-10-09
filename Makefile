# OpenBiota Gut Health Test
#
# `make` alone runs the full screen on ./fastq.
# `make help` lists everything.

SHELL      := /bin/bash
PYTHON     ?= python3.13
VENV       := .venv
PY         := $(VENV)/bin/python
OPENBIOTA       := $(VENV)/bin/openbiota
THREADS    ?= $(shell getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)
FASTQ_DIR  ?= fastq
RESULTS    ?= results
SUBSAMPLE  ?= 100000
COHORT_SAMPLES ?= 34

# --------------------------------------------------------------------------- #
# Pinned external tool versions.
#
# Every version here is exact, never a floating minimum. MetaPhlAn is the
# reason: the 4.1 and 4.2 marker workflows produce incompatible consensus
# markers, so `metaphlan>=4.1` once installed 4.2.5 beside a 4.1.1 database and
# put a silently-wrong strain lane one command away.
# --------------------------------------------------------------------------- #
MPA3_VERSION      ?= 3.1.0
MPA4_VERSION      ?= 4.1.1
MPA4_INDEX        ?= mpa_vJun23_CHOCOPhlAnSGB_202403
STRAIN_VENV       ?= .venv-strain
STRAIN_PYTHON     ?= python3.12
# inStrain needs its own interpreter. It pins biopython<=1.74 and genuinely
# requires that API (`Bio.codonalign.codonalphabet`, removed in 1.78), and
# biopython 1.74's C source does not compile on 3.12 or 3.13
# (Bio/Align/_aligners.c uses PyUnicode_WCHAR_KIND, removed in 3.12). 3.9 is
# the newest interpreter where the pin installs cleanly. Relaxing the pin was
# tried and fails at import, so this is not a conservative bound.
INSTRAIN_VENV     ?= .venv-instrain
INSTRAIN_PYTHON   ?= python3.9
INSTRAIN_VERSION  ?= 1.10.0
ECTYPER_VERSION   ?= 2.0.0
# Whole-genome lane. Sylph is Rust and built from its tagged source; the
# GTDB species sketch is 20 GB and fetched once. Both are pinned: the lane's
# cache key includes them, so a change here invalidates every profile.
SYLPH_VERSION     ?= 1.0.0
SYLPH_DB          ?= gtdb-r232-c200-dbv2.syl2db
SYLPH_DB_URL      ?= https://faust.compbio.cs.cmu.edu/sylph-stuff/$(SYLPH_DB)
SYLPH_TAX_VERSION ?= 1.9.2
# Pinned EukDetect2 code commit for the mycobiome module (database: Zenodo 19056625).
EUKDETECT_COMMIT  ?= 8d69014727b2c5956de30b0811644b6c3a7bd4f3

# v0.8.3 extension. MICOM is pinned to the version the published replay
# fixture was produced with (BUILD_SPEC_v0.8.3 section 11.7): updating it
# requires its own versioned validation, because a solver or objective change
# moves every predicted flux.
# dbCAN snapshot for the A01 substrate panel. One dated release: the HMM,
# the substrate mapping and the CAZy set must come from the same snapshot,
# because mixing an HMM with an older mapping silently changes what every
# family is taken to mean.
DBCAN_RELEASE := V13
DBCAN_BASE := https://bcb.unl.edu/dbCAN2/download/Databases/$(DBCAN_RELEASE)
DBCAN_DIR := refs/dbcan

MICOM_VERSION     ?= 0.37.0
AGORA2_QZA        ?= agora201_refseq216_species_1.qza
AGORA2_URL        ?= https://zenodo.org/records/7739096/files/$(AGORA2_QZA)?download=1
AGORA2_MD5        ?= 6a13da2a3fd9b1bc059987d6344f32b6
# The four pinned S10 files, with the SHA-256 the specification records. A
# download that does not match is the wrong bytes, not a new revision.
S10_COMMIT        ?= a494f4e21c8510cf02618af1861f8649d92fe698
S10_RAW           ?= https://raw.githubusercontent.com/Gibbons-Lab/2024_probiotic_engraftment/$(S10_COMMIT)
S5_DATA_SHA       ?= 19ad4e8e2eea635f3c577df0d1bc04ba852249f368bac25d2519abbe71bd8c3e
ARIVALE_SHA       ?= e811818f16a1067a4e8d365a047bb0dbdf4c0fc3d185b89e840a53c3dd0edeea
EUROPEAN_MED_SHA  ?= 0299fdd9b9d1ba2938681a6b078d54ddde73add1912a6b3582c8c9dffbea8f61
HIGHFIBRE_MED_SHA ?= b29a95fc3bf1639b78bca65ea32d78ee7498753aa37ad653d366df9e2198d526
KLEBORATE_VERSION ?= 3.2.4
# The intersection of every pin in the shared venv: ECTyper accepts <= 1.85,
# MIDAS v3 accepts <= 1.83, and Kleborate would raise it past both. `pip check`
# runs after the installs so a future conflict surfaces instead of silently
# breaking one tool.
BIOPYTHON_VERSION ?= 1.83
CARGO_BIN          ?= $(HOME)/.cargo/bin
# Vendored under vendor/ because Homebrew cannot supply these cleanly here:
# `brew install blast` refuses over a libproj.a conflict with `proj`, and NCBI
# ships AMRFinderPlus binaries for Linux only.
BLAST_VERSION      ?= 2.17.0
AMRFINDER_VERSION  ?= 4.2.7
STXTYPER_VERSION   ?= 1.0.45
MIDAS_VERSION      ?= 1.0.1
SAMESTR_VERSION    ?= 1.2025.111
# Hard dependencies of the typing tools. mash is not in core homebrew.
MASH_FORMULA       ?= brewsci/bio/mash
MINIMAP2_FORMULA   ?= minimap2
STRAINGE_VERSION  ?= 1.3.9
# TRACS is pinned by commit, not tag: the v1.1.4 tag's package metadata still
# reports 1.1.3 upstream, so `pip show tracs` will say 1.1.3 even though this
# is the 1.1.4 source. Spec §11.4 asks for a locked tag/commit for exactly
# this reason -- the interfaces move faster than the version string.
TRACS_VERSION     ?= 1.1.4
TRACS_COMMIT      ?= de9282a64fba8e3d4386ce114e50d55e7fddc440
TRACS_REPORTS_VERSION ?= 1.1.3

.DEFAULT_GOAL := run
.PHONY: tools-expanded refs-expanded refs-globdb-genomes expansion-cohorts expansion-crosswalks expansion-reconcile expansion extension-refs substrate-refs help setup doctor panels build-db cohort validate depth-check all smoke run compact json test lint clean clean-all profilers strain-tools strain-system-deps strain-tools-report strain-tools-lock refs-strain refs-strain-verify refs-ctnpc instrain kleborate-optional vendor-tools genome-lane mycobiome-refs

help:
	@echo "OpenBiota Gut Health Test — targets"
	@echo
	@echo "  make setup       create $(VENV) and install the package (editable)"
	@echo "  make doctor      check dependencies, hardware and inputs"
	@echo "  make profiles    list the disease similarity profiles"
	@echo "  make panels      list the gene panels"
	@echo "  make build-db    fetch reference sets and build the DIAMOND database"
	@echo "  make cohort      build percentile reference ranges from public metagenomes"
	@echo "  make validate    measure sensitivity, specificity and accuracy"
	@echo "  make depth-check test whether the sequencing depth was sufficient"
	@echo "  make profilers   create .venv-mpa3 / .venv-mpa4 at pinned versions"
	@echo "  make strain-tools install the strain, typing and comparison tools (pinned)"
	@echo "  make strain-tools-lock  freeze the installed strain tools to a lockfile"
	@echo "  make tools-expanded     install every pinned tool for expanded detection (spec 0.8.4) and lock versions"
	@echo "  make refs-expanded      fetch every expanded reference once, verified and locked (resumable)"
	@echo "  make expansion          crosswalks + supplement reconciliation + rescue panel (needs refs-expanded)"
	@echo "  make refs-globdb-genomes fetch + unpack the 272 GB GlobDB genome archive (confirmation of GlobDB-native clusters)"
	@echo "  make expansion-cohorts  build the expanded lanes' own reference populations for percentiles"
	@echo "  make refs-strain fetch the strain/typing reference panels"
	@echo "  make refs-strain-verify   check every strain/typing reference against its lock"
	@echo "  make refs-ctnpc  delimit the CTnPc accessory region against the study's comparators"
	@echo "  make genome-lane install sylph and fetch the GTDB R232 species sketch (20 GB)"
	@echo "  make mycobiome-refs       fungal references for the mycobiome module (~90 GB, resumable, locked)"
	@echo "  make extension-refs       metabolic model library and replay fixtures for model-assisted scenarios"
	@echo "  make substrate-refs       substrate gene references for the synbiotic scenarios (one dated release)"
	@echo "  make vendor-tools         pinned binaries not packaged for this platform (blastn, makeblastdb, raxmlHPC) into vendor/"
	@echo "  make strain-system-deps   system libraries the strain tools need"
	@echo "  make strain-tools-report  what of the strain toolchain is installed, at which versions"
	@echo "  make instrain             inStrain in its own environment (optional strain comparison)"
	@echo "  make kleborate-optional   Kleborate for Klebsiella typing (optional)"
	@echo "  make expansion-crosswalks / expansion-reconcile   the two halves of make expansion"
	@echo
	@echo "  make profilers          create .venv-mpa3 / .venv-mpa4 (MetaPhlAn 3 scoring lane, 4 extended catalogue)"
	@echo "  make taxonomic-cohort   build the species reference cohort (cMD)"
	@echo "  make host-index         build the human-read filter (one-time, slow)"
	@echo "  make validate-profiles  profile AUCs on labelled cohorts"
	@echo
	@echo "  make all         every build step in order, then run"
	@echo "  make smoke       fast end-to-end test on the first $(SUBSAMPLE) read pairs"
	@echo "  make run         full screen on $(FASTQ_DIR)  [default target]"
	@echo "  make run-all     screen every sample in $(FASTQ_DIR) and compare them"
	@echo "  make compact     full screen, one-screen output"
	@echo "  make json        full screen, JSON to stdout"
	@echo "  make test        unit tests"
	@echo "  make lint        ruff, if installed"
	@echo "  make prune       list superseded reference databases (make prune-all deletes them)"
	@echo "  make clean       remove $(RESULTS)/"
	@echo "  make clean-all   remove $(RESULTS)/, refs/ and $(VENV)/"
	@echo
	@echo "  variables: THREADS=$(THREADS) FASTQ_DIR=$(FASTQ_DIR) SUBSAMPLE=$(SUBSAMPLE)"
	@echo "             COHORT_SAMPLES=$(COHORT_SAMPLES)"

$(VENV)/bin/activate: pyproject.toml
	$(PYTHON) -m venv $(VENV)
	$(PY) -m pip install --quiet --upgrade pip
	$(PY) -m pip install --quiet -e '.[dev]'
	@touch $(VENV)/bin/activate

setup: $(VENV)/bin/activate
	@$(OPENBIOTA) --version

doctor: setup
	@$(OPENBIOTA) doctor --fastq-dir $(FASTQ_DIR)

panels: setup
	@$(OPENBIOTA) panels

build-db: setup
	@$(OPENBIOTA) build-db --threads $(THREADS)

cohort: setup
	@$(OPENBIOTA) cohort --samples $(COHORT_SAMPLES) --threads $(THREADS)

validate: setup
	@$(OPENBIOTA) validate --threads $(THREADS)

depth-check: setup
	@$(OPENBIOTA) depth-check --fastq-dir $(FASTQ_DIR) --threads $(THREADS)

# MetaPhlAn 3 and 4 cannot share an environment. Version 3 is the scoring
# lane (the reference cohort was profiled with it); version 4 is the extended
# SGB catalogue, reported but never scored. Databases download on first use
# into refs/metaphlan_db and refs/metaphlan4_db (~25 GB for version 4).
profilers:
	@test -x .venv-mpa3/bin/metaphlan || ( $(PYTHON) -m venv .venv-mpa3 && .venv-mpa3/bin/pip install --quiet 'metaphlan==$(MPA3_VERSION)' )
	@test -x .venv-mpa4/bin/metaphlan || ( $(PYTHON) -m venv .venv-mpa4 && .venv-mpa4/bin/pip install --quiet 'metaphlan==$(MPA4_VERSION)' )
	@.venv-mpa3/bin/metaphlan --version
	@.venv-mpa4/bin/metaphlan --version
	@# The 4.1 and 4.2 marker workflows are not interchangeable: consensus
	@# markers built by one cannot be read by the other. Fail loudly rather
	@# than let a floating install silently break the strain lane.
	@.venv-mpa4/bin/metaphlan --version 2>&1 | grep -q '$(MPA4_VERSION)' || \
		( echo "FATAL: .venv-mpa4 is not MetaPhlAn $(MPA4_VERSION); the strain lane requires it" && exit 1 )

# --------------------------------------------------------------------------- #
# strain, typing and comparison tools
#
# Every external tool the resolution layer can use, pinned. `make strain-tools`
# is idempotent: it skips anything already installed at the pinned version.
# Tools that are absent are not a failure -- the resolution census reports them
# as reference_unavailable -- but nothing should ever be installed by hand.
# --------------------------------------------------------------------------- #

strain-tools: strain-system-deps
	@command -v $(STRAIN_PYTHON) >/dev/null 2>&1 || \
		( echo "FATAL: $(STRAIN_PYTHON) required (inStrain's biopython will not build on 3.13)" && exit 1 )
	@test -d $(STRAIN_VENV) || ( $(STRAIN_PYTHON) -m venv $(STRAIN_VENV) && \
		$(STRAIN_VENV)/bin/pip install --quiet --upgrade pip wheel )
	@# Installed one at a time on purpose: a single unavailable package must
	@# not abort the whole set. Anything that fails is reported by the
	@# resolution census as reference_unavailable.
	@for spec in 'ectyper==$(ECTYPER_VERSION)' 'strainge==$(STRAINGE_VERSION)' \
	             'midasv3==$(MIDAS_VERSION)' 'samestr==$(SAMESTR_VERSION)'; do \
		$(STRAIN_VENV)/bin/pip show "$${spec%%==*}" >/dev/null 2>&1 && continue; \
		echo "  installing $$spec"; \
		$(STRAIN_VENV)/bin/pip install --quiet "$$spec" >/dev/null 2>&1 || \
			echo "    UNAVAILABLE: $$spec"; \
	done
	@# TRACS is not on PyPI, and is pinned by commit rather than tag.
	@$(STRAIN_VENV)/bin/pip show tracs >/dev/null 2>&1 || \
		$(STRAIN_VENV)/bin/pip install --quiet \
			'git+https://github.com/gtonkinhill/tracs@$(TRACS_COMMIT)' >/dev/null 2>&1 || \
		echo "    UNAVAILABLE: tracs $(TRACS_VERSION) ($(TRACS_COMMIT))"
	@$(MAKE) --no-print-directory instrain kleborate-optional
	@$(MAKE) --no-print-directory strain-tools-report

instrain:
	@command -v $(INSTRAIN_PYTHON) >/dev/null 2>&1 || \
		( echo "    UNAVAILABLE: inStrain needs $(INSTRAIN_PYTHON)" && exit 0 )
	@test -x $(INSTRAIN_VENV)/bin/inStrain && exit 0 || true
	@$(INSTRAIN_PYTHON) -m venv $(INSTRAIN_VENV) && \
		$(INSTRAIN_VENV)/bin/pip install --quiet --upgrade pip wheel && \
		$(INSTRAIN_VENV)/bin/pip install --quiet 'inStrain==$(INSTRAIN_VERSION)' || \
		echo "    UNAVAILABLE: inStrain==$(INSTRAIN_VERSION)"

# Kleborate's kaptive dependency builds `rammappy`, a Rust + CMake extension.
# Three things have to line up, and each failed in turn while getting this to
# install:
#   1. Rust >= 1.85 -- an older toolchain cannot parse rammappy's Cargo.toml.
#   2. numba/llvmlite -- these must come from wheels. Building llvmlite from
#      source needs LLVM <= 22, and Homebrew ships 23, so a source build fails
#      with "llvmlite only officially supports 22" however LLVM_DIR is set.
#   3. biopython pinned to $(BIOPYTHON_VERSION) -- Kleborate pulls 1.88, which
#      breaks ECTyper's requirement of <= 1.85 in the same venv.
kleborate-optional:
	@test -x $(STRAIN_VENV)/bin/kleborate && exit 0 || true
	@command -v $(CARGO_BIN)/cargo >/dev/null 2>&1 || \
		( echo "    UNAVAILABLE: kleborate (needs rust >= 1.85; see rustup.rs)" && exit 0 )
	@PATH="$(CARGO_BIN):$$PATH" \
		$(STRAIN_VENV)/bin/pip install --quiet \
			--only-binary llvmlite,numba 'kleborate==$(KLEBORATE_VERSION)' \
		>/dev/null 2>&1 || echo "    UNAVAILABLE: kleborate==$(KLEBORATE_VERSION)"
	@# Re-pin biopython: Kleborate raises it past what ECTyper accepts.
	@test -x $(STRAIN_VENV)/bin/kleborate && \
		$(STRAIN_VENV)/bin/pip install --quiet 'biopython==$(BIOPYTHON_VERSION)' || true
	@$(STRAIN_VENV)/bin/pip check >/dev/null 2>&1 || \
		( echo "    WARNING: $(STRAIN_VENV) has broken requirements:" && \
		  $(STRAIN_VENV)/bin/pip check | head -4 )

# Tools that Homebrew either does not carry or cannot install without
# disturbing existing formulae. Installed under vendor/ so nothing on the
# system is relinked.
#
# BLAST+ is the case that matters: `brew install blast` refuses because it
# conflicts with `proj` (both ship libproj.a), and resolving that would mean
# unlinking a formula this project does not own. NCBI's own tarball avoids the
# question entirely.
## genome-lane: install Sylph and fetch the GTDB species sketch (20 GB)
genome-lane:
	@mkdir -p vendor refs/sylph/tax
	@test -x vendor/sylph/bin/sylph && vendor/sylph/bin/sylph --version | grep -q "$(SYLPH_VERSION)" || ( \
		echo "  building sylph $(SYLPH_VERSION) from source"; \
		cargo install --git https://github.com/bluenote-1577/sylph --tag v$(SYLPH_VERSION) \
			--root vendor/sylph sylph )
	@$(VENV)/bin/python -c "import sylph_tax" 2>/dev/null || ( \
		echo "  installing sylph-tax $(SYLPH_TAX_VERSION)"; \
		$(VENV)/bin/pip install -q "git+https://github.com/bluenote-1577/sylph-tax@v$(SYLPH_TAX_VERSION)" )
	@test -s refs/sylph/tax/gtdb_r232_metadata.tsv.gz || \
		$(VENV)/bin/sylph-tax download --download-to refs/sylph/tax
	@test -s refs/sylph/$(SYLPH_DB) || ( \
		echo "  fetching $(SYLPH_DB) (20 GB, resumable)"; \
		curl -sSL --retry 5 --retry-delay 10 -C - -o refs/sylph/$(SYLPH_DB).part "$(SYLPH_DB_URL)" && \
		mv refs/sylph/$(SYLPH_DB).part refs/sylph/$(SYLPH_DB) )
	@echo "  genome lane ready: sylph $(SYLPH_VERSION), $(SYLPH_DB), GTDB r232 taxonomy"

## mycobiome-refs: fungal references for the mycobiome module (~90 GB): EukDetect2 DB, FungiGutDB,
## CGF + RefSeq fungal genomes, strain-panel assemblies and isolate reads, decoys; then lock, sketch,
## species index and strain panels. Resumable; every file is hashed into refs/mycobiome/reference_lock.json.
mycobiome-refs: setup
	@mkdir -p refs/mycobiome/logs vendor/ncbi-datasets
	@test -x vendor/ncbi-datasets/datasets || ( \
		echo "  fetching NCBI datasets CLI"; \
		curl -sSL --retry 5 -o vendor/ncbi-datasets/datasets \
			"https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets" && chmod +x vendor/ncbi-datasets/datasets )
	@test -d vendor/eukdetect/EukDetect/.git || ( mkdir -p vendor/eukdetect && cd vendor/eukdetect && \
		git clone -q https://github.com/allind/EukDetect.git && cd EukDetect && git checkout -q $(EUKDETECT_COMMIT) )
	@test -x .venv-eukdetect/bin/eukdetect || ( \
		echo "  installing EukDetect2 into .venv-eukdetect"; \
		python3.11 -m venv .venv-eukdetect 2>/dev/null || python3 -m venv .venv-eukdetect; \
		.venv-eukdetect/bin/pip install -q --upgrade pip && \
		.venv-eukdetect/bin/pip install -q "snakemake>=9,<10" pysam biopython ete3 pyyaml six numpy && \
		.venv-eukdetect/bin/pip install -q -e vendor/eukdetect/EukDetect )
	@which bcftools >/dev/null || ( echo "  bcftools missing: brew install bcftools"; exit 1 )
	@refs/mycobiome/fetch_deposits.sh
	@refs/mycobiome/fetch_genomes.sh
	@refs/mycobiome/fetch_li2022.sh
	@$(OPENBIOTA) mycobiome references prepare --threads $(THREADS)
	@$(OPENBIOTA) mycobiome references status

## extension-refs: assets for the v0.8.3 report extension (~300 MB).
##   MICOM $(MICOM_VERSION) and an open QP solver, the AGORA2 species-model
##   library, and the four pinned public files the model replay is checked
##   against. Every download is hashed; a mismatch fails the target rather
##   than installing bytes nobody has verified.
extension-refs: setup
	@mkdir -p refs/micom/s10
	@$(VENV)/bin/python -c "import micom" 2>/dev/null || ( \
		echo "  installing micom $(MICOM_VERSION) and the osqp/highs solvers"; \
		$(VENV)/bin/pip install -q "micom==$(MICOM_VERSION)" osqp highspy )
	@test -s refs/micom/$(AGORA2_QZA) || ( \
		echo "  fetching AGORA2 species models (~264 MB, resumable)"; \
		curl -sSL --retry 5 --retry-delay 10 -C - -o refs/micom/$(AGORA2_QZA).part "$(AGORA2_URL)" && \
		mv refs/micom/$(AGORA2_QZA).part refs/micom/$(AGORA2_QZA) )
	@printf '%s  refs/micom/%s\n' "$(AGORA2_MD5)" "$(AGORA2_QZA)" | md5sum -c - 2>/dev/null \
		|| [ "$$(md5 -q refs/micom/$(AGORA2_QZA))" = "$(AGORA2_MD5)" ] \
		|| ( echo "  AGORA2 checksum mismatch; refusing to install"; exit 1 )
	@set -e; for spec in \
		"figures/source_data/S5_Data.xlsx:S5_Data.xlsx:$(S5_DATA_SHA)" \
		"data/arivale-metagenomes.csv:arivale-metagenomes.csv:$(ARIVALE_SHA)" \
		"european_medium.csv:european_medium.csv:$(EUROPEAN_MED_SHA)" \
		"high_fiber_medium.csv:high_fiber_medium.csv:$(HIGHFIBRE_MED_SHA)"; do \
		path=$${spec%%:*}; rest=$${spec#*:}; name=$${rest%%:*}; want=$${rest#*:}; \
		test -s "refs/micom/s10/$$name" || curl -sSL --retry 5 -o "refs/micom/s10/$$name" "$(S10_RAW)/$$path"; \
		got=$$(shasum -a 256 "refs/micom/s10/$$name" | cut -d' ' -f1); \
		[ "$$got" = "$$want" ] || ( echo "  $$name checksum mismatch: $$got"; exit 1 ); \
	done
	@$(OPENBIOTA) extension-assets verify --spec 0.8.3

## substrate-refs: the dbCAN snapshot behind the A01 substrate panel (~2 GB).
##   One dated release, fetched together, with a snapshot.json recording every
##   file's size and checksum, and a reconciliation against the fourteen
##   substrate definitions printed at the end.
substrate-refs: setup
	@echo "installing the dbCAN $(DBCAN_RELEASE) snapshot for the A01 substrate panel"
	@mkdir -p $(DBCAN_DIR)
	@# One release, fetched together. The substrate mapping is what the
	@# fourteen substrate views are reconciled against; a mapping from a
	@# different release would quietly redefine every family.
	@for f in fam-substrate-mapping.tsv CAZyDB.fa dbCAN-HMMdb-$(DBCAN_RELEASE).txt \
	          dbCAN_sub.hmm PUL.faa dbCAN-PUL.xlsx; do \
		test -s $(DBCAN_DIR)/$$f || ( \
			echo "  fetching $$f"; \
			curl -fsSL --retry 5 --retry-delay 10 -C - \
				-o $(DBCAN_DIR)/$$f.part "$(DBCAN_BASE)/$$f" \
			&& mv $(DBCAN_DIR)/$$f.part $(DBCAN_DIR)/$$f ) || \
			echo "  WARNING: $$f is not in this release; recording the gap"; \
	done
	@# Record exactly what was installed, so a report can name its release.
	@$(VENV)/bin/python -c "import hashlib, json, pathlib; d = pathlib.Path('$(DBCAN_DIR)'); \
		print(json.dumps({'release': '$(DBCAN_RELEASE)', 'base_url': '$(DBCAN_BASE)', \
		'files': {p.name: {'bytes': p.stat().st_size, \
		'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} \
		for p in sorted(d.iterdir()) if p.is_file() and p.suffix != '.json'}}, indent=1))" \
		> $(DBCAN_DIR)/snapshot.json
	@echo "  wrote $(DBCAN_DIR)/snapshot.json"
	@$(VENV)/bin/python -c "from pathlib import Path; \
		from openbiota.extension.substrates import load_dbcan_mapping, validate_against_dbcan; \
		import json; m = load_dbcan_mapping(Path('$(DBCAN_DIR)/fam-substrate-mapping.tsv')); \
		r = validate_against_dbcan(m); \
		print('  reconciled', r['n_families_in_mapping'], 'families; agrees:', r['agrees']); \
		[print('   ', k, v) for k, v in r['families_not_in_mapping'].items()]"


vendor-tools:
	@mkdir -p vendor
	@test -x vendor/ncbi-blast-$(BLAST_VERSION)+/bin/blastn || ( \
		echo "  installing BLAST+ $(BLAST_VERSION)"; \
		curl -sSL --max-time 900 -o vendor/blast.tar.gz \
			"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$(BLAST_VERSION)/ncbi-blast-$(BLAST_VERSION)+-x64-macosx.tar.gz" && \
		tar xzf vendor/blast.tar.gz -C vendor && rm -f vendor/blast.tar.gz )
	@# AMRFinderPlus ships Linux binaries only; build from the pinned source.
	@# The bundled stxtyper is a git submodule absent from the release tarball,
	@# so that one target fails while amrfinder itself builds fine.
	@test -x vendor/amr-amrfinder_v$(AMRFINDER_VERSION)/amrfinder || ( \
		echo "  building AMRFinderPlus $(AMRFINDER_VERSION)"; \
		curl -sSL --max-time 600 -o vendor/amr-src.tar.gz \
			"https://github.com/ncbi/amr/archive/refs/tags/amrfinder_v$(AMRFINDER_VERSION).tar.gz" && \
		tar xzf vendor/amr-src.tar.gz -C vendor && rm -f vendor/amr-src.tar.gz && \
		( cd vendor/amr-amrfinder_v$(AMRFINDER_VERSION) && $(MAKE) -j8 amrfinder amrfinder_update amrfinder_index >/dev/null 2>&1 || true ) )
	@test -d vendor/amr-amrfinder_v$(AMRFINDER_VERSION)/data/latest || ( \
		echo "  fetching the AMRFinderPlus database"; \
		cd vendor/amr-amrfinder_v$(AMRFINDER_VERSION) && \
		PATH="$(CURDIR)/vendor/ncbi-blast-$(BLAST_VERSION)+/bin:$$PATH" \
			./amrfinder_update -d ./data >/dev/null 2>&1 || \
			echo "    UNAVAILABLE: AMRFinderPlus database" )
	@# StxTyper: Linux binaries only, and AMRFinderPlus bundles it as a git
	@# submodule that release tarballs omit. Built separately from its own tag.
	@test -x vendor/stxtyper-$(STXTYPER_VERSION)/stxtyper || ( \
		echo "  building StxTyper $(STXTYPER_VERSION)"; \
		curl -sSL --max-time 600 -o vendor/stx-src.tar.gz \
			"https://github.com/ncbi/stxtyper/archive/refs/tags/v$(STXTYPER_VERSION).tar.gz" && \
		tar xzf vendor/stx-src.tar.gz -C vendor && rm -f vendor/stx-src.tar.gz && \
		( cd vendor/stxtyper-$(STXTYPER_VERSION) && $(MAKE) -j8 >/dev/null 2>&1 || true ) )
	@printf "  %-22s %s\n" "blastn" \
		"$$(vendor/ncbi-blast-$(BLAST_VERSION)+/bin/blastn -version 2>/dev/null | head -1 || echo MISSING)"
	@printf "  %-22s %s\n" "stxtyper" \
		"$$(vendor/stxtyper-$(STXTYPER_VERSION)/stxtyper --version 2>/dev/null || echo MISSING)"
	@printf "  %-22s %s\n" "amrfinder" \
		"$$(vendor/amr-amrfinder_v$(AMRFINDER_VERSION)/amrfinder --version 2>/dev/null || echo MISSING)"

strain-system-deps: vendor-tools
	@command -v brew >/dev/null 2>&1 || ( echo "Homebrew required for system deps" && exit 1 )
	@# samtools is the one hard dependency: sample2markers.py shells out to it,
	@# so without it the marker lane emits SAM and then dies at the consensus
	@# step. Everything below it is optional.
	@command -v samtools >/dev/null 2>&1 || brew install samtools
	@command -v samtools >/dev/null 2>&1 || \
		( echo "FATAL: samtools is required by the marker lane" && exit 1 )
	@# mash and minimap2 are hard dependencies of the typing tools, not
	@# optional extras: ECTyper aborts without mash, and Kleborate needs
	@# both. mash lives in brewsci/bio rather than core homebrew.
	@command -v mash >/dev/null 2>&1 || ( brew tap brewsci/bio >/dev/null 2>&1; \
		brew install $(MASH_FORMULA) >/dev/null 2>&1 ) || \
		echo "    UNAVAILABLE: mash (ECTyper and Kleborate cannot run without it)"
	@command -v minimap2 >/dev/null 2>&1 || \
		brew install $(MINIMAP2_FORMULA) >/dev/null 2>&1 || \
		echo "    UNAVAILABLE: minimap2 (Kleborate cannot run without it)"

strain-tools-report:
	@echo "strain and typing tools:"
	@printf "  %-22s %s\n" "samtools" "$$(samtools --version 2>/dev/null | head -1 || echo MISSING)"
	@printf "  %-22s %s\n" "minimap2" "$$(minimap2 --version 2>/dev/null || echo MISSING)"
	@printf "  %-22s %s\n" "mash" "$$(mash --version 2>/dev/null | head -1 || echo MISSING)"
	@for t in ectyper kleborate straingst tracs; do \
		printf "  %-22s %s\n" "$$t" "$$(test -x $(STRAIN_VENV)/bin/$$t && echo $(STRAIN_VENV)/bin/$$t || echo MISSING)"; \
	done
	@printf "  %-22s %s\n" "inStrain" "$$(test -x $(INSTRAIN_VENV)/bin/inStrain && echo $(INSTRAIN_VENV)/bin/inStrain || echo MISSING)"
	@for t in midas samestr; do \
		printf "  %-22s %s\n" "$$t" "$$(test -x $(STRAIN_VENV)/bin/$$t && echo $(STRAIN_VENV)/bin/$$t || echo MISSING)"; \
	done
	@# TRACS reports an older version string than the commit it was built
	@# from; check the pin, not the self-reported number.
	@if test -x $(STRAIN_VENV)/bin/tracs; then \
		got=$$($(STRAIN_VENV)/bin/tracs --version 2>/dev/null | tail -1 | tr -d 'tracs '); \
		printf "  %-22s reports %s, built from %s (tag v%s)\n" \
			"tracs pin" "$$got" "$(TRACS_COMMIT)" "$(TRACS_VERSION)"; \
	fi
	@echo
	@echo "Anything MISSING is reported by the resolution census as"
	@echo "reference_unavailable -- never as a negative result."

refs-strain: setup
	@$(PY) tools/fetch_strain_refs.py
	@echo
	@echo "lockfile: refs/strain_refs.lock.json"

refs-strain-verify: setup
	@$(PY) tools/fetch_strain_refs.py --verify-only

# Delimit the CTnPc accessory region by presence/absence against the study's
# own CTnPc-negative comparators. Depends only on the fetched panels and
# DIAMOND, both of which the pipeline already has.
refs-ctnpc: refs-strain
	@$(PY) tools/delimit_ctnpc.py --threads $(THREADS)

strain-tools-lock:
	@test -x $(STRAIN_VENV)/bin/pip || ( echo "run 'make strain-tools' first" && exit 1 )
	@$(STRAIN_VENV)/bin/pip freeze > requirements-strain.lock
	@echo "wrote requirements-strain.lock ($$(wc -l < requirements-strain.lock | tr -d ' ') packages)"

taxonomic-cohort: setup
	@$(OPENBIOTA) taxonomic-cohort --stability

host-index: setup
	@$(OPENBIOTA) build-host-index --threads $(THREADS)

validate-profiles: setup
	@$(OPENBIOTA) validate-profiles

profiles: setup
	@$(OPENBIOTA) profiles

all: build-db cohort validate taxonomic-cohort validate-profiles run

smoke: setup
	@$(OPENBIOTA) run --fastq-dir $(FASTQ_DIR) --out $(RESULTS) --threads $(THREADS) \
		--subsample $(SUBSAMPLE)

run: setup
	@$(OPENBIOTA) run --fastq-dir $(FASTQ_DIR) --out $(RESULTS) --threads $(THREADS)

run-all: setup
	@$(OPENBIOTA) run-all --fastq-dir $(FASTQ_DIR) --out $(RESULTS) --threads $(THREADS)

compact: setup
	@$(OPENBIOTA) run --fastq-dir $(FASTQ_DIR) --out $(RESULTS) --threads $(THREADS) --compact

json: setup
	@$(OPENBIOTA) run --fastq-dir $(FASTQ_DIR) --out $(RESULTS) --threads $(THREADS) --json --quiet

test: setup
	@$(PY) -m pytest

lint: setup
	@$(PY) -m ruff check openbiota tests 2>/dev/null || echo "ruff not installed: pip install ruff"

prune: setup
	@$(OPENBIOTA) prune

prune-all: setup
	@$(OPENBIOTA) prune --alignments --yes

clean:
	rm -rf $(RESULTS)

clean-all: clean
	rm -rf refs $(VENV) .pytest_cache .ruff_cache
	find . -name __pycache__ -type d -prune -exec rm -rf {} +



# --------------------------------------------------------------------------- #
# Spec 0.8.4 - expanded gut bacterial recognition
# --------------------------------------------------------------------------- #
MPA42_VENV   ?= .venv-mpa42
MOTUS_VENV   ?= .venv-motus
SINGLEM_VENV ?= .venv-singlem
SPADES_VERSION ?= 4.3.0

tools-expanded: setup
	@echo "expanded detection tools (pinned):"
	@# aria2 gives the reference fetcher parallel range requests (a 290 GB archive in minutes, not hours)
	@command -v aria2c >/dev/null 2>&1 || HOMEBREW_NO_AUTO_UPDATE=1 HOMEBREW_NO_INSTALL_UPGRADE=1 brew install aria2
	@command -v seqkit >/dev/null 2>&1 || HOMEBREW_NO_AUTO_UPDATE=1 HOMEBREW_NO_INSTALL_UPGRADE=1 brew install seqkit
	@# MetaPhlAn 4.2.6 from the GitHub tag (its setup.py still reports 4.2.5); Jun23's 4.1.1 env is untouched.
	@test -x $(MPA42_VENV)/bin/metaphlan || ( python3.9 -m venv $(MPA42_VENV) && $(MPA42_VENV)/bin/pip install -q --upgrade pip && \
		$(MPA42_VENV)/bin/pip install -q "metaphlan @ https://github.com/biobakery/MetaPhlAn/archive/refs/tags/4.2.6.tar.gz" )
	@# mOTUs 4.1.0 needs python >= 3.12
	@test -x $(MOTUS_VENV)/bin/motus || ( python3.12 -m venv $(MOTUS_VENV) && $(MOTUS_VENV)/bin/pip install -q --upgrade pip && \
		$(MOTUS_VENV)/bin/pip install -q "motus-tool @ https://github.com/motu-tool/mOTUs/archive/refs/tags/4.1.0.tar.gz" )
	@# SingleM 0.21.4 and its binary dependencies
	@command -v orfm >/dev/null 2>&1 || HOMEBREW_NO_AUTO_UPDATE=1 HOMEBREW_NO_INSTALL_UPGRADE=1 brew install brewsci/bio/orfm
	@test -x $$HOME/.cargo/bin/smafa || $$HOME/.cargo/bin/cargo install smafa
	@test -x $$HOME/.cargo/bin/mfqe || $$HOME/.cargo/bin/cargo install mfqe
	@test -x $(SINGLEM_VENV)/bin/singlem || ( python3.12 -m venv $(SINGLEM_VENV) && $(SINGLEM_VENV)/bin/pip install -q --upgrade pip && \
		$(SINGLEM_VENV)/bin/pip install -q "singlem==0.21.4" )
	@# Kraken2, Bracken, skani, fastANI, bwa
	@for t in kraken2 bracken skani fastani bwa; do command -v $$t >/dev/null 2>&1 || command -v fastANI >/dev/null 2>&1 || \
		HOMEBREW_NO_AUTO_UPDATE=1 HOMEBREW_NO_INSTALL_UPGRADE=1 brew install $$t; done
	@# SPAdes (macOS binary release) for targeted assembly
	@test -x vendor/SPAdes-$(SPADES_VERSION)-Darwin/bin/spades.py || ( mkdir -p vendor && \
		curl -sSL --max-time 600 -o vendor/spades.tar.gz "https://github.com/ablab/spades/releases/download/v$(SPADES_VERSION)/SPAdes-$(SPADES_VERSION)-Darwin-x86_64.tar.gz" && \
		tar xzf vendor/spades.tar.gz -C vendor && rm -f vendor/spades.tar.gz )
	@# PhyloPhlAn (inside StrainPhlAn) looks for raxmlHPC, blastn and makeblastdb by name.
	@mkdir -p vendor/bin
	@for f in mafft trimal brewsci/bio/raxml brewsci/bio/fasttree; do HOMEBREW_NO_AUTO_UPDATE=1 HOMEBREW_NO_INSTALL_UPGRADE=1 brew install $$f >/dev/null 2>&1 || true; done
	@test -e vendor/bin/raxmlHPC || ln -sf /usr/local/bin/raxmlHPC-PTHREADS-SSE3 vendor/bin/raxmlHPC
	@# relative links, so the checkout can be moved or renamed without breaking them
	@test -e vendor/bin/blastn || ln -sf "../ncbi-blast-$(BLAST_VERSION)+/bin/blastn" vendor/bin/blastn
	@test -e vendor/bin/makeblastdb || ln -sf "../ncbi-blast-$(BLAST_VERSION)+/bin/makeblastdb" vendor/bin/makeblastdb
	@$(PY) scripts/lock_tools.py

refs-expanded: setup
	@$(PY) scripts/fetch_expanded_refs.py --jobs 6
	@# The Jan26 marker index (publisher-built bowtie2 shards, ~40 GB) is installed by MetaPhlAn itself.
	@$(MPA42_VENV)/bin/metaphlan --install --index mpa_vJan26_CHOCOPhlAnSGB_202605 --db_dir refs/metaphlan4_db --nproc $(THREADS)
	@# mOTUs database: unpack the verified archive
	@test -d refs/motus/db || ( mkdir -p refs/motus/db && tar xzf refs/motus/db_mOTU.tar.gz -C refs/motus/db )
	@# mOTUs checks for the marker its own downloader writes; ours was md5-verified by the fetcher.
	@test -f refs/motus/db/db_mOTU/db_mOTU.downloaded || echo "downloaded by scripts/fetch_expanded_refs.py (Zenodo 20322482, md5 verified)" > refs/motus/db/db_mOTU/db_mOTU.downloaded
	@echo "expanded references fetched; locks in refs/expanded/locks"

expansion-crosswalks: setup
	@$(PY) -m openbiota.expansion.crosswalk

expansion-reconcile: setup
	@$(PY) -m openbiota.expansion.reconcile $(THREADS)

expansion: expansion-crosswalks expansion-reconcile
	@echo "crosswalks: refs/crosswalk/*.summary.json; reconciliation: refs/expanded/reconciliation/summary.json"

refs-globdb-genomes: setup
	@# 272 GB one-time archive of every GlobDB representative genome; deferred from refs-expanded on purpose.
	@$(PY) scripts/fetch_expanded_refs.py --only globdb_genomes --include-deferred --jobs 1
	@command -v pv >/dev/null 2>&1 || HOMEBREW_NO_AUTO_UPDATE=1 HOMEBREW_NO_INSTALL_UPGRADE=1 brew install pv
	@$(PY) scripts/unpack_globdb_genomes.py --rate-limit 150

expansion-cohorts: setup
	@# A reference population per expanded lane, profiled with that lane's own engine (see scripts/build_expansion_cohort.py).
	@$(PY) scripts/build_expansion_cohort.py --lanes globdb --pairs 3000000
