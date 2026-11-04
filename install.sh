#!/usr/bin/env bash
# install.sh - from a fresh clone to a working ./openbiota, in one command.
#
#   ./install.sh              everything a complete report needs (~15 GB)
#   ./install.sh --full       plus every optional detection lane (~1.5 TB)
#   ./install.sh --list       what the steps are and which are already done
#   ./install.sh --help       every option
#
# Re-running is always safe: each step checks what it already produced and
# skips it, so an interrupted install picks up where it stopped. Everything
# each step prints goes to install.log; the terminal shows progress and, if
# something fails, the part of the log that explains why.
#
# Works on macOS (Intel and Apple silicon) and Linux (Debian, Ubuntu, Fedora,
# and anything else through the conda fallback). Written for bash 3.2, which
# is what macOS ships, so it does not use arrays-of-arrays, associative
# arrays, mapfile or lower-casing expansions.

set -u
set -o pipefail

# --------------------------------------------------------------------------- #
# What gets installed, pinned
# --------------------------------------------------------------------------- #

MPA3_VERSION="3.1.0"                       # scoring lane; the cohort was profiled with it
MPA3_INDEX="mpa_v31_CHOCOPhlAn_201901"     # its database
MICROMAMBA_URL="https://micro.mamba.pm/api/micromamba"

# Steps in the standard install, in order. Each is a pair of functions below:
# done_<step> (0 when nothing is left to do) and run_<step>.
STANDARD_STEPS="preflight system python launcher profiler data metaphlan_db diamond_db models host_index check"
# Added by --full. Each is a Makefile target, so the pinned versions live in
# one place; these are the slow, large lanes and take most of a day.
FULL_STEPS="profiler4 vendor_tools strain_tools refs_strain pathogens substrate_refs tools_expanded refs_expanded expansion genome_lane mycobiome_refs"

# Free disk each tier needs, in GB, with headroom for the downloads' temporaries.
DISK_STANDARD_GB=20
DISK_NO_HOST_GB=10
DISK_FULL_GB=1500

# --------------------------------------------------------------------------- #
# Options
# --------------------------------------------------------------------------- #

FULL=0
SKIP_HOST_INDEX=0
NO_PATH=0
FORCE=0
ASSUME_YES=0
LIST_ONLY=0
DRY_RUN=0
OFFLINE=0
ONLY=""
SKIP=""
THREADS=""
NO_COLOR_FLAG=${NO_COLOR:-}

usage() {
    cat <<'EOF'
usage: ./install.sh [options]

  --full               also install every optional detection lane (~1.5 TB,
                       most of a day; see docs/INSTALL.md for what each adds)
  --skip-host-index    leave out the human-read filter (saves ~7 GB and an
                       hour or more; the report then cannot state the human
                       DNA fraction, which for stool is usually under 1%)
  --no-path            do not link `openbiota` into ~/.local/bin or touch
                       your shell start-up file; bin/openbiota still works
  --only STEP[,STEP]   run just these steps
  --skip STEP[,STEP]   run everything except these
  --threads N          CPU threads for the builds (default: all of them)
  --force              redo steps even when they look complete
  --yes                do not stop to ask about installing system packages
  --offline            fail fast instead of downloading anything
  --dry-run            show what would run, run nothing
  --list               show every step and whether it is done
  --no-color           plain output
  -h, --help           this

Steps, in order:
  preflight     operating system, disk, memory, network, permissions
  system        git, curl, python 3.11+, diamond, bowtie2, samtools, fastp
  python        the .venv with openbiota installed into it
  launcher      `openbiota` on your PATH (a link to bin/openbiota)
  profiler      MetaPhlAn 3 in its own .venv-mpa3
  data          the reference data this project computed and the UniProt
                responses behind the protein database (from data/, 250 MB)
  metaphlan_db  the MetaPhlAn 3 marker database (0.4 GB download, 3.5 GB built)
  diamond_db    the protein search database (built locally from the cached
                UniProt responses; ~2 min)
  models        the metabolic model library for the synbiotic scenarios
  host_index    the human-genome read filter (7 GB download + build)
  check         verify every piece, end to end

With --full, afterwards:
  profiler4 vendor_tools strain_tools refs_strain pathogens substrate_refs
  tools_expanded refs_expanded expansion genome_lane mycobiome_refs
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --full) FULL=1 ;;
        --skip-host-index) SKIP_HOST_INDEX=1 ;;
        --no-path) NO_PATH=1 ;;
        --only) shift; ONLY=${1:-}; [ -n "$ONLY" ] || { echo "--only needs a step name" >&2; exit 2; } ;;
        --only=*) ONLY=${1#--only=} ;;
        --skip) shift; SKIP=${1:-}; [ -n "$SKIP" ] || { echo "--skip needs a step name" >&2; exit 2; } ;;
        --skip=*) SKIP=${1#--skip=} ;;
        --threads) shift; THREADS=${1:-}; [ -n "$THREADS" ] || { echo "--threads needs a number" >&2; exit 2; } ;;
        --threads=*) THREADS=${1#--threads=} ;;
        --force) FORCE=1 ;;
        --yes|-y) ASSUME_YES=1 ;;
        --offline) OFFLINE=1 ;;
        --dry-run) DRY_RUN=1 ;;
        --list) LIST_ONLY=1 ;;
        --no-color) NO_COLOR_FLAG=1 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown option: $1 (try --help)" >&2; exit 2 ;;
    esac
    shift
done

case "${THREADS:-1}" in
    ''|*[!0-9]*) echo "--threads must be a whole number, not '$THREADS'" >&2; exit 2 ;;
esac

# --------------------------------------------------------------------------- #
# Where we are
# --------------------------------------------------------------------------- #

# The script's own directory is the repository, however it was invoked.
REPO=$(cd "$(dirname "$0")" && pwd -P)
cd "$REPO" || { echo "cannot enter $REPO" >&2; exit 1; }
LOG="$REPO/install.log"
LOCK="$REPO/.install.lock"
VENV="$REPO/.venv"
MPA3="$REPO/.venv-mpa3"
VENDOR_BIN="$REPO/vendor/bin"
CONDA="$REPO/vendor/conda"
mkdir -p "$VENDOR_BIN" "$REPO/refs" "$REPO/fastq" "$REPO/results"

# Tools the installer puts under vendor/ win over whatever else is on PATH,
# which is also what ./openbiota does, so the install and the run agree.
export PATH="$VENDOR_BIN:$CONDA/bin:$PATH"
export PIP_DISABLE_PIP_VERSION_CHECK=1
export HOMEBREW_NO_AUTO_UPDATE=1
export HOMEBREW_NO_INSTALL_CLEANUP=1
export HOMEBREW_NO_ENV_HINTS=1
export DEBIAN_FRONTEND=noninteractive

# --------------------------------------------------------------------------- #
# Terminal
# --------------------------------------------------------------------------- #

if [ -t 1 ] && [ -z "$NO_COLOR_FLAG" ] && [ "${TERM:-dumb}" != "dumb" ]; then
    TTY=1
    BOLD=$(printf '\033[1m'); DIM=$(printf '\033[2m'); RESET=$(printf '\033[0m')
    GREEN=$(printf '\033[32m'); YELLOW=$(printf '\033[33m'); RED=$(printf '\033[31m'); CYAN=$(printf '\033[36m')
else
    TTY=0
    BOLD=""; DIM=""; RESET=""; GREEN=""; YELLOW=""; RED=""; CYAN=""
fi
COLUMNS_NOW=$(tput cols 2>/dev/null || echo 80)
[ "$COLUMNS_NOW" -ge 40 ] 2>/dev/null || COLUMNS_NOW=80
# Block glyphs and ticks only where the terminal can show them.
case "${LC_ALL:-${LC_CTYPE:-${LANG:-}}}" in
    *UTF-8*|*utf-8*|*UTF8*|*utf8*) FULLC="█"; EMPTYC="░"; TICK="✓"; CROSS="✗" ;;
    *) FULLC="#"; EMPTYC="-"; TICK="ok"; CROSS="FAILED" ;;
esac

log() { printf '%s %s\n' "$(date '+%H:%M:%S')" "$*" >>"$LOG"; }
say() { printf '%s\n' "$*"; log "$*"; }
ok() { printf '  %s%s%s %s\n' "$GREEN" "$TICK" "$RESET" "$*"; log "OK $*"; }
warn() { printf '  %s!%s %s\n' "$YELLOW" "$RESET" "$*"; log "WARN $*"; }
note() { printf '  %s%s%s\n' "$DIM" "$*" "$RESET"; log "$*"; }

die() {
    printf '\n  %s%s %s%s\n' "$RED" "$CROSS" "$*" "$RESET" >&2
    log "FAIL $*"
    exit 1
}

# A step failed inside a command: show why, from the log, then stop.
fail_with_log() {
    printf '\n  %s%s %s%s\n' "$RED" "$CROSS" "$1" "$RESET"
    log "FAIL $1"
    if [ -s "$LOG" ]; then
        printf '\n  %sThe last lines of %s:%s\n' "$DIM" "install.log" "$RESET"
        tail -n 25 "$LOG" | sed 's/^/    /'
    fi
    printf '\n  Fix the cause and run %s./install.sh%s again; finished steps are kept.\n' "$BOLD" "$RESET"
    exit 1
}

# Overall progress: one line per step, with the step counter on the left.
STEP_N=0
STEP_TOTAL=0
STEP_T0=0
bar() {
    # bar DONE TOTAL -> a 20-cell bar
    local done=$1 total=$2 cells=20 filled i out=""
    filled=$(( total > 0 ? done * cells / total : 0 ))
    i=0
    while [ $i -lt $cells ]; do
        if [ $i -lt $filled ]; then out="${out}${FULLC}"; else out="${out}${EMPTYC}"; fi
        i=$(( i + 1 ))
    done
    printf '%s' "$out"
}
begin_step() {
    STEP_N=$(( STEP_N + 1 ))
    STEP_T0=$(date +%s)
    printf '\n%s%s%s %s%2d/%d%s  %s%s%s\n' \
        "$CYAN" "$(bar $(( STEP_N - 1 )) "$STEP_TOTAL")" "$RESET" \
        "$DIM" "$STEP_N" "$STEP_TOTAL" "$RESET" "$BOLD" "$1" "$RESET"
    log "==== step $STEP_N/$STEP_TOTAL: $1"
}
end_step() {
    local secs=$(( $(date +%s) - STEP_T0 ))
    if [ $secs -ge 90 ]; then note "done in $(( secs / 60 )) min"; fi
}

elapsed_text() {
    local s=$1
    if [ "$s" -lt 60 ]; then printf '%ds' "$s"
    elif [ "$s" -lt 3600 ]; then printf '%dm%02ds' $(( s / 60 )) $(( s % 60 ))
    else printf '%dh%02dm' $(( s / 3600 )) $(( (s % 3600) / 60 )); fi
}

# run_logged LABEL COMMAND...  - runs the command with all its output in the
# log, and keeps one live line on the terminal with the elapsed time and the
# command's latest line, so a forty-minute build never looks hung.
run_logged() {
    local label=$1; shift
    log "\$ $*"
    if [ $DRY_RUN -eq 1 ]; then
        note "would run: $*"
        return 0
    fi
    if [ $TTY -eq 0 ]; then
        printf '  %s ...\n' "$label"
        "$@" >>"$LOG" 2>&1
        return $?
    fi
    local t0 pid status last width line
    t0=$(date +%s)
    "$@" >>"$LOG" 2>&1 &
    pid=$!
    width=$(( COLUMNS_NOW - 4 ))
    while kill -0 "$pid" 2>/dev/null; do
        # A tool that redraws one line with carriage returns leaves them all
        # on the log's last line; the newest segment is the one to show.
        last=$(tail -n 1 "$LOG" 2>/dev/null | tr '\r' '\n' | tail -n 1 | tr '\t' ' ')
        line=$(printf '  %s  %s  %s' "$label" "$(elapsed_text $(( $(date +%s) - t0 )))" "$last")
        printf '\r\033[2K%s' "$(printf '%s' "$line" | cut -c1-"$width")"
        sleep 1
    done
    wait "$pid"; status=$?
    printf '\r\033[2K'
    return $status
}

# download URL DEST  - resumable, retried, with a live size readout.
download() {
    local url=$1 dest=$2 part label
    part="$dest.part"
    label="downloading $(basename "$dest")"
    [ $OFFLINE -eq 1 ] && { warn "offline: cannot fetch $url"; return 1; }
    if ! run_logged "$label" curl -fL --retry 5 --retry-delay 5 -C - \
            --connect-timeout 30 -o "$part" "$url"; then
        return 1
    fi
    mv -f "$part" "$dest"
}

confirm() {
    # confirm QUESTION -> 0 for yes. Non-interactive runs say yes with --yes
    # and otherwise refuse, so a cron job or CI never hangs on a prompt.
    [ $ASSUME_YES -eq 1 ] && return 0
    if [ ! -t 0 ]; then
        warn "$1 - not a terminal, so cannot ask; re-run with --yes to allow it"
        return 1
    fi
    printf '  %s [Y/n] ' "$1"
    local answer
    read -r answer
    case "$answer" in ''|y|Y|yes|YES|Yes) return 0 ;; *) return 1 ;; esac
}

have() { command -v "$1" >/dev/null 2>&1; }

# --------------------------------------------------------------------------- #
# Platform
# --------------------------------------------------------------------------- #

OS=$(uname -s)
ARCH=$(uname -m)
PKG=""           # brew | apt | dnf | ""
SUDO=""
case "$OS" in
    Darwin) ;;
    Linux) ;;
    *) die "unsupported operating system: $OS (macOS and Linux are supported)" ;;
esac

if [ "$(id -u)" -eq 0 ]; then
    # Root is normal inside a container and wrong on a desktop: Homebrew
    # refuses to run as root, and files the installer writes would then be
    # unreadable to the user who clones the repository.
    if [ "$OS" = "Darwin" ]; then die "do not run the installer as root on macOS (Homebrew refuses it); run it as your own user"; fi
else
    if have sudo; then SUDO="sudo"; fi
fi

detect_pkg() {
    if [ "$OS" = "Darwin" ]; then
        have brew && PKG=brew
    elif have apt-get; then PKG=apt
    elif have dnf; then PKG=dnf
    elif have brew; then PKG=brew
    fi
}
detect_pkg

cpu_count() {
    getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4
}
[ -n "$THREADS" ] || THREADS=$(cpu_count)

free_gb() {
    # free space on the filesystem holding the repository, whole GB
    df -Pk "$REPO" | awk 'NR==2 { printf "%d", $4 / 1048576 }'
}
memory_gb() {
    if [ "$OS" = "Darwin" ]; then
        sysctl -n hw.memsize 2>/dev/null | awk '{ printf "%d", $1 / 1073741824 }'
    else
        awk '/MemTotal/ { printf "%d", $2 / 1048576 }' /proc/meminfo 2>/dev/null || echo 0
    fi
}

# The newest python 3.11+ on the machine, by command name.
find_python() {
    local candidate
    for candidate in python3.13 python3.12 python3.11 python3; do
        if have "$candidate" && "$candidate" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 11) else 1)' 2>/dev/null; then
            command -v "$candidate"
            return 0
        fi
    done
    return 1
}

# --------------------------------------------------------------------------- #
# Package installation, three ways: the system's manager, then conda for
# whatever that could not supply, then a clear message.
# --------------------------------------------------------------------------- #

APT_UPDATED=0
pkg_install() {
    # pkg_install NAME... (manager package names)
    [ $OFFLINE -eq 1 ] && { warn "offline: cannot install $*"; return 1; }
    case "$PKG" in
        brew) run_logged "brew install $*" brew install "$@" ;;
        apt|dnf)
            # Ask for the password here, in the open, before the command's
            # output is sent to the log and the live status line would paint
            # over sudo's prompt.
            if [ -n "$SUDO" ] && [ $DRY_RUN -eq 0 ]; then
                printf '  %s%s needs administrator rights to install packages%s\n' "$DIM" "$PKG" "$RESET"
                sudo -v || { warn "sudo refused; the tools will come from conda instead"; return 1; }
            fi
            if [ "$PKG" = apt ]; then
                if [ $APT_UPDATED -eq 0 ]; then
                    run_logged "apt-get update" $SUDO apt-get update -qq || true
                    APT_UPDATED=1
                fi
                run_logged "apt-get install $*" $SUDO apt-get install -y -qq --no-install-recommends "$@"
            else
                run_logged "dnf install $*" $SUDO dnf install -y -q "$@"
            fi ;;
        *) return 1 ;;
    esac
}

ensure_micromamba() {
    have micromamba && return 0
    [ -x "$VENDOR_BIN/micromamba" ] && return 0
    [ $OFFLINE -eq 1 ] && return 1
    local platform
    case "$OS-$ARCH" in
        Darwin-arm64) platform=osx-arm64 ;;
        Darwin-x86_64) platform=osx-64 ;;
        Linux-x86_64) platform=linux-64 ;;
        Linux-aarch64) platform=linux-aarch64 ;;
        *) return 1 ;;
    esac
    note "fetching micromamba so the missing tools can come from conda-forge/bioconda"
    local tmp="$REPO/vendor/.micromamba.tmp"
    rm -rf "$tmp"; mkdir -p "$tmp"
    download "$MICROMAMBA_URL/$platform/latest" "$tmp/mm.tar.bz2" || return 1
    ( cd "$tmp" && tar -xjf mm.tar.bz2 bin/micromamba 2>>"$LOG" ) || return 1
    mv -f "$tmp/bin/micromamba" "$VENDOR_BIN/micromamba" && chmod +x "$VENDOR_BIN/micromamba"
    rm -rf "$tmp"
    [ -x "$VENDOR_BIN/micromamba" ]
}

conda_install() {
    # conda_install PKG... into vendor/conda, then expose the binaries.
    ensure_micromamba || return 1
    local mm
    mm=$(command -v micromamba || echo "$VENDOR_BIN/micromamba")
    export MAMBA_ROOT_PREFIX="$REPO/vendor/.mamba"
    if [ ! -d "$CONDA" ]; then
        run_logged "conda: creating vendor/conda" "$mm" create -y -q -p "$CONDA" -c conda-forge -c bioconda "$@" || return 1
    else
        run_logged "conda: installing $*" "$mm" install -y -q -p "$CONDA" -c conda-forge -c bioconda "$@" || return 1
    fi
}

# Map a tool to its package name per manager. Everything the standard install
# needs is in Homebrew core, Debian/Ubuntu universe and bioconda.
pkg_name_for() {
    local tool=$1
    case "$PKG:$tool" in
        apt:diamond) echo diamond-aligner ;;
        apt:python) echo "python3 python3-venv python3-pip" ;;
        dnf:python) echo "python3 python3-pip" ;;
        brew:python) echo "python@3.13" ;;
        *:python) echo python3 ;;
        *) echo "$tool" ;;
    esac
}
conda_name_for() {
    case "$1" in
        python) echo "python=3.12 pip" ;;
        *) echo "$1" ;;
    esac
}

ensure_tool() {
    # ensure_tool NAME  - installed by whatever route works; 0 when present.
    local tool=$1
    if [ "$tool" = python ]; then
        find_python >/dev/null 2>&1 && return 0
    else
        have "$tool" && return 0
    fi
    note "$tool is not installed"
    if [ -n "$PKG" ]; then
        # shellcheck disable=SC2046
        pkg_install $(pkg_name_for "$tool") || warn "$PKG could not install $tool; trying conda"
    fi
    if [ "$tool" = python ]; then find_python >/dev/null 2>&1 && return 0; else have "$tool" && return 0; fi
    # shellcheck disable=SC2046
    conda_install $(conda_name_for "$tool") || true
    hash -r
    if [ "$tool" = python ]; then find_python >/dev/null 2>&1 && return 0; else have "$tool" && return 0; fi
    return 1
}

# --------------------------------------------------------------------------- #
# Steps
# --------------------------------------------------------------------------- #

# -- preflight ---------------------------------------------------------------
done_preflight() { return 1; }   # always runs; it is the check
run_preflight() {
    ok "$OS $ARCH, $(cpu_count) CPUs, $(memory_gb) GB memory, using $THREADS threads"
    if [ "$OS" = "Darwin" ]; then
        if ! xcode-select -p >/dev/null 2>&1; then
            warn "the Xcode command line tools are not installed; macOS needs them for git and compilers"
            if [ $DRY_RUN -eq 0 ] && confirm "install them now? (macOS opens a dialog)"; then
                xcode-select --install 2>>"$LOG" || true
                die "finish the Xcode command line tools dialog, then run ./install.sh again"
            else
                die "run 'xcode-select --install', then run ./install.sh again"
            fi
        fi
        if [ -z "$PKG" ]; then
            warn "Homebrew is not installed; it is how macOS gets diamond, bowtie2, samtools and fastp"
            if [ $DRY_RUN -eq 0 ] && [ $OFFLINE -eq 0 ] && confirm "install Homebrew now? (its installer asks for your password)"; then
                /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" </dev/tty || die "the Homebrew installer did not finish"
                # Homebrew's own instruction, applied here so this shell sees it.
                for prefix in /opt/homebrew /usr/local; do
                    [ -x "$prefix/bin/brew" ] && eval "$("$prefix/bin/brew" shellenv)" && break
                done
                detect_pkg
            fi
            [ -n "$PKG" ] || die "install Homebrew from https://brew.sh, then run ./install.sh again"
        fi
        ok "Homebrew at $(brew --prefix)"
    else
        if [ -z "$PKG" ]; then
            warn "no apt-get, dnf or brew found; the tools will come from conda-forge/bioconda instead"
        elif [ "$PKG" != brew ] && [ "$(id -u)" -ne 0 ] && [ -z "$SUDO" ]; then
            warn "no sudo: cannot use $PKG; the tools will come from conda-forge/bioconda instead"
            PKG=""
        else
            ok "packages via $PKG"
        fi
    fi

    local need free
    if [ $FULL -eq 1 ]; then need=$DISK_FULL_GB
    elif [ $SKIP_HOST_INDEX -eq 1 ]; then need=$DISK_NO_HOST_GB
    else need=$DISK_STANDARD_GB; fi
    free=$(free_gb)
    if [ "$free" -lt "$need" ]; then
        die "only $free GB free on the drive holding this folder; the install needs about $need GB. Free some space, or clone onto a larger drive."
    fi
    ok "$free GB free (need about $need GB)"

    local mem; mem=$(memory_gb)
    if [ "$mem" -gt 0 ] && [ "$mem" -lt 8 ]; then
        warn "$mem GB memory is tight: the human-genome index build needs about 8 GB; consider --skip-host-index"
    fi

    if [ $OFFLINE -eq 0 ] && [ $DRY_RUN -eq 0 ]; then
        if ! online; then
            die "cannot reach the internet (pypi.org); the install downloads references. Check the connection, or use --offline to only verify what is here."
        fi
        ok "internet reachable"
    fi
    if [ ! -w "$REPO" ]; then die "$REPO is not writable"; fi
    case "$REPO" in
        *" "*|*"	"*)
            # bowtie2's wrapper scripts, MetaPhlAn and several of the full-tier
            # tools split paths on whitespace. Say so now, not forty minutes in.
            die "the folder path contains a space ($REPO); the bioinformatics tools cannot cope with that. Move or re-clone the project into a path without spaces, e.g. ~/openbiota-gut-health-test" ;;
    esac
}

# Reachability without assuming curl is installed yet.
online() {
    if have curl; then
        curl -fsI --connect-timeout 15 https://pypi.org/simple/ >/dev/null 2>&1 && return 0
        return 1
    fi
    local py
    py=$(find_python 2>/dev/null || command -v python3 2>/dev/null) || return 0   # cannot tell; let the download say
    "$py" -c 'import urllib.request; urllib.request.urlopen("https://pypi.org/simple/", timeout=15)' >/dev/null 2>&1
}

# -- system ------------------------------------------------------------------
SYSTEM_TOOLS="git curl python diamond bowtie2 samtools fastp"
done_system() {
    local t py
    for t in $SYSTEM_TOOLS; do
        if [ "$t" = python ]; then
            py=$(find_python) || return 1
            "$py" -c 'import venv, ensurepip' >/dev/null 2>&1 || return 1
        else have "$t" || return 1; fi
    done
    return 0
}
run_system() {
    local t missing=""
    for t in $SYSTEM_TOOLS; do
        if [ "$t" = python ]; then find_python >/dev/null 2>&1 || missing="$missing $t"
        else have "$t" || missing="$missing $t"; fi
    done
    if [ -z "$missing" ]; then ok "every tool is already installed"; else
        note "to install:$missing"
        if [ -n "$PKG" ] && [ "$PKG" != brew ] && [ $DRY_RUN -eq 0 ]; then
            confirm "install them with $PKG (uses sudo)?" || die "cannot continue without:$missing"
        fi
        for t in $missing; do
            ensure_tool "$t" || die "could not install $t by any route. Install it yourself (see docs/INSTALL.md), then run ./install.sh again."
        done
    fi
    # Debian and Ubuntu ship python3 without the venv module; a python that
    # cannot make an environment is, for our purposes, not installed.
    local py; py=$(find_python)
    if ! "$py" -c 'import venv, ensurepip' >/dev/null 2>&1; then
        note "$py cannot create virtual environments (no venv/ensurepip module)"
        if [ "$PKG" = apt ]; then
            local minor; minor=$("$py" -c 'import sys; print("%d.%d" % sys.version_info[:2])')
            pkg_install "python$minor-venv" || pkg_install python3-venv || true
        fi
        if ! "$py" -c 'import venv, ensurepip' >/dev/null 2>&1; then
            conda_install python=3.12 pip || true
            hash -r
            py=$(find_python) || die "no usable python"
            "$py" -c 'import venv, ensurepip' >/dev/null 2>&1 \
                || die "$py has no venv module. On Debian/Ubuntu: sudo apt-get install python3-venv; then run ./install.sh again."
        fi
    fi
    # Some of these are easy to have at the wrong version; say what we got.
    ok "python   $("$(find_python)" --version 2>&1 | sed 's/Python //')  ($(find_python))"
    ok "diamond  $(diamond version 2>/dev/null | head -1 | sed 's/diamond version //')"
    ok "bowtie2  $(bowtie2 --version 2>/dev/null | head -1 | sed 's/.*version //')"
    ok "samtools $(samtools --version 2>/dev/null | head -1 | sed 's/samtools //')"
    ok "fastp    $(fastp --version 2>&1 | head -1 | sed 's/fastp //')"
    # pigz is optional and halves gzip time; make is only needed for --full.
    have pigz || { [ -n "$PKG" ] && pkg_install pigz >/dev/null 2>&1 || true; }
    if [ $FULL -eq 1 ] && ! have make; then
        ensure_tool make || die "--full needs make"
    fi
}

# -- python ------------------------------------------------------------------
done_python() {
    [ -x "$VENV/bin/openbiota" ] && "$VENV/bin/python" -c 'import openbiota, reportlab, yaml' >/dev/null 2>&1
}
run_python() {
    local py; py=$(find_python) || die "no python 3.11+ (the system step should have installed one)"
    if [ ! -x "$VENV/bin/python" ] || ! "$VENV/bin/python" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 11) else 1)' 2>/dev/null; then
        [ -d "$VENV" ] && warn "replacing an unusable .venv"
        rm -rf "$VENV"
        run_logged "creating .venv with $py" "$py" -m venv "$VENV" || fail_with_log "python -m venv failed (on Debian/Ubuntu: sudo apt-get install python3-venv)"
    fi
    run_logged "upgrading pip" "$VENV/bin/python" -m pip install --quiet --upgrade pip || fail_with_log "pip could not upgrade itself"
    # --prefer-binary: a slightly older wheel beats a newer source release
    # that would compile for half an hour (or fail for want of a compiler).
    run_logged "installing openbiota and its python dependencies" \
        "$VENV/bin/python" -m pip install --quiet --prefer-binary -e ".[dev]" || fail_with_log "pip install of openbiota failed"
    done_python || fail_with_log ".venv/bin/openbiota does not import after installation"
    ok "openbiota $("$VENV/bin/openbiota" --version 2>&1 | head -1 | sed 's/openbiota //') in .venv"
}

# -- launcher ----------------------------------------------------------------
# bin/openbiota always works from the repository. This step also makes plain
# `openbiota` work from any new terminal: a link in ~/.local/bin, and that
# folder on the PATH if it is not already, added to the shell's start-up file
# inside a marked block so it is easy to find and remove.
LOCAL_BIN="${XDG_BIN_HOME:-${HOME:-/nonexistent}/.local/bin}"
PATH_MARK="# added by openbiota install.sh"

shell_rc() {
    # The start-up file of the user's login shell. zsh is macOS's default,
    # bash is Linux's; anything else is told what to add and left alone.
    case "$(basename "${SHELL:-}")" in
        zsh) echo "$HOME/.zshrc" ;;
        bash) if [ "$OS" = "Darwin" ]; then echo "$HOME/.bash_profile"; else echo "$HOME/.bashrc"; fi ;;
        *) echo "" ;;
    esac
}
path_has_local_bin() {
    case ":$PATH:" in *":$LOCAL_BIN:"*) return 0 ;; esac
    local rc; rc=$(shell_rc)
    [ -n "$rc" ] && [ -f "$rc" ] && grep -qF "$PATH_MARK" "$rc"
}
done_launcher() {
    [ -x "$REPO/bin/openbiota" ] || return 1
    [ $NO_PATH -eq 1 ] && return 0
    [ "$(readlink "$LOCAL_BIN/openbiota" 2>/dev/null)" = "$REPO/bin/openbiota" ] && path_has_local_bin
}
run_launcher() {
    chmod +x "$REPO/bin/openbiota" 2>/dev/null || true
    [ -x "$REPO/bin/openbiota" ] || die "bin/openbiota is missing from this checkout"
    ok "bin/openbiota"
    if [ $NO_PATH -eq 1 ]; then note "--no-path: not linking it onto your PATH"; return 0; fi
    if [ $DRY_RUN -eq 1 ]; then note "would link $LOCAL_BIN/openbiota -> bin/openbiota"; return 0; fi
    if [ -z "${HOME:-}" ] || [ ! -w "$HOME" ]; then
        warn "no writable home folder; use bin/openbiota (or pass --no-path to silence this)"
        return 0
    fi
    mkdir -p "$LOCAL_BIN" 2>/dev/null || { warn "cannot create $LOCAL_BIN; use bin/openbiota"; return 0; }
    if [ -e "$LOCAL_BIN/openbiota" ] && [ ! -L "$LOCAL_BIN/openbiota" ]; then
        warn "$LOCAL_BIN/openbiota exists and is not a link; leaving it alone. bin/openbiota still works."
        return 0
    fi
    ln -sfn "$REPO/bin/openbiota" "$LOCAL_BIN/openbiota" || { warn "could not link into $LOCAL_BIN; use bin/openbiota"; return 0; }
    ok "openbiota -> $LOCAL_BIN/openbiota"
    if path_has_local_bin; then
        return 0
    fi
    local rc; rc=$(shell_rc)
    if [ -z "$rc" ]; then
        warn "add $LOCAL_BIN to your PATH to use plain 'openbiota'; bin/openbiota works regardless"
        return 0
    fi
    {
        printf '\n%s\n' "$PATH_MARK"
        # shellcheck disable=SC2016
        printf 'export PATH="%s:$PATH"\n' "$LOCAL_BIN"
    } >>"$rc" || { warn "could not write to $rc; use bin/openbiota"; return 0; }
    export PATH="$LOCAL_BIN:$PATH"
    NEW_TERMINAL_NOTE=1
    ok "added $LOCAL_BIN to the PATH in $(basename "$rc") (takes effect in new terminals)"
}
NEW_TERMINAL_NOTE=0

# -- profiler (MetaPhlAn 3) --------------------------------------------------
done_profiler() {
    [ -x "$MPA3/bin/metaphlan" ] && "$MPA3/bin/metaphlan" --version 2>&1 | grep -q "$MPA3_VERSION"
}
run_profiler() {
    local py; py=$(find_python)
    if [ -d "$MPA3" ] && ! done_profiler; then
        warn "replacing .venv-mpa3: it is not MetaPhlAn $MPA3_VERSION"
        rm -rf "$MPA3"
    fi
    if [ ! -x "$MPA3/bin/python" ]; then
        run_logged "creating .venv-mpa3" "$py" -m venv "$MPA3" || fail_with_log "python -m venv failed for .venv-mpa3"
        run_logged "upgrading pip" "$MPA3/bin/python" -m pip install --quiet --upgrade pip wheel || fail_with_log "pip upgrade failed in .venv-mpa3"
    fi
    run_logged "installing MetaPhlAn $MPA3_VERSION" \
        "$MPA3/bin/python" -m pip install --quiet --prefer-binary "metaphlan==$MPA3_VERSION" || fail_with_log "pip install metaphlan==$MPA3_VERSION failed"
    done_profiler || fail_with_log ".venv-mpa3/bin/metaphlan is not version $MPA3_VERSION"
    ok "MetaPhlAn $MPA3_VERSION in .venv-mpa3"
}

# -- data (the archive in data/) ---------------------------------------------
# Asked of the bundles' own manifests: every item each one carries is under
# refs/. Needs only a python with the standard library, so this step has no
# dependency on the .venv.
any_python() {
    if [ -x "$VENV/bin/python" ]; then echo "$VENV/bin/python"; else find_python; fi
}
done_data() {
    local py; py=$(any_python) || return 1
    "$py" scripts/derived_refs.py check --quiet >/dev/null 2>&1
}
run_data() {
    [ -s "$REPO/data/derived-refs.manifest.json" ] \
        || die "data/derived-refs.manifest.json is missing from this checkout; re-clone the repository (data/ is tracked in git)"
    local py
    py=$(any_python) || die "no python 3.11+ to unpack the archives with"
    run_logged "verifying the archives in data/" "$py" scripts/derived_refs.py verify \
        || fail_with_log "an archive in data/ does not match its manifest; the checkout may be incomplete (git lfs is not involved; a plain re-clone fixes a truncated download)"
    run_logged "unpacking data/derived-refs" "$py" scripts/derived_refs.py restore derived-refs \
        || fail_with_log "could not unpack data/derived-refs.tar.gz"
    "$py" scripts/derived_refs.py check derived-refs --quiet >>"$LOG" 2>&1 \
        || fail_with_log "refs/ is still missing files after unpacking"
    ok "reference ranges, scoring cohort, age model, crosswalks and panels in refs/"
    # The UniProt responses make the next-but-one step a two-minute local
    # build instead of six hours of fetching; without them it still works.
    if run_logged "unpacking data/uniprot-cache" "$py" scripts/derived_refs.py restore uniprot-cache; then
        ok "UniProt responses for the protein database in refs/cache/"
    else
        warn "the UniProt cache did not unpack; build-db will fetch from UniProt instead (slow)"
    fi
}

# ready NAME - the pipeline's own answer to "is this piece installed?"
ready() {
    [ -x "$VENV/bin/python" ] && "$VENV/bin/python" -m openbiota.installcheck "$1" >/dev/null 2>&1
}

# -- metaphlan_db ------------------------------------------------------------
done_metaphlan_db() { ready metaphlan_db; }
run_metaphlan_db() {
    [ $OFFLINE -eq 1 ] && die "offline: the MetaPhlAn database is not installed and cannot be downloaded"
    mkdir -p "$REPO/refs/metaphlan_db"
    note "0.4 GB from the MetaPhlAn servers, then indexed with bowtie2 to 3.5 GB; the index build is the slow part (20-40 min)"
    run_logged "installing $MPA3_INDEX" "$MPA3/bin/metaphlan" --install --index "$MPA3_INDEX" \
        --bowtie2db "$REPO/refs/metaphlan_db" --nproc "$THREADS" \
        || fail_with_log "MetaPhlAn could not install its database (a dropped connection is the usual cause; re-run to resume)"
    done_metaphlan_db || fail_with_log "the MetaPhlAn database is not complete after installing"
    ok "MetaPhlAn database $MPA3_INDEX"
}

# -- diamond_db --------------------------------------------------------------
done_diamond_db() { ready diamond_db; }
run_diamond_db() {
    [ $OFFLINE -eq 1 ] && [ -z "$(ls "$REPO"/refs/cache/*.tsv.gz 2>/dev/null)" ] && die "offline: the protein database is not built and the UniProt cache is not unpacked"
    if [ -n "$(ls "$REPO"/refs/cache/*.tsv.gz 2>/dev/null)" ]; then
        note "building the DIAMOND index from the cached UniProt responses"
    else
        note "no UniProt cache in refs/cache: fetching every panel's reference proteins from UniProt (hours)"
    fi
    run_logged "openbiota build-db" "$VENV/bin/openbiota" build-db --threads "$THREADS" \
        || fail_with_log "openbiota build-db failed"
    done_diamond_db || fail_with_log "the DIAMOND database is not complete after building"
    ok "protein search database built"
}

# -- models (MICOM + AGORA2) ---------------------------------------------------
done_models() { ready models; }
run_models() {
    [ $OFFLINE -eq 1 ] && die "offline: the model library is not installed"
    have make || ensure_tool make || die "this step uses make"
    run_logged "make extension-refs" make --no-print-directory extension-refs THREADS="$THREADS" \
        || fail_with_log "make extension-refs failed"
    done_models || fail_with_log "the model library is incomplete after installing"
    ok "MICOM and the AGORA2 model library"
}

# -- host_index --------------------------------------------------------------
done_host_index() {
    [ $SKIP_HOST_INDEX -eq 1 ] && return 0
    ready host_index
}
run_host_index() {
    if [ $SKIP_HOST_INDEX -eq 1 ]; then note "skipped (--skip-host-index)"; return 0; fi
    [ $OFFLINE -eq 1 ] && die "offline: the host index is not built and needs the GRCh38 download"
    note "downloads the human genome (0.9 GB) and indexes it; about an hour on a desktop, longer on a laptop"
    run_logged "openbiota build-host-index" "$VENV/bin/openbiota" build-host-index --threads "$THREADS" \
        || fail_with_log "openbiota build-host-index failed"
    done_host_index || fail_with_log "the host index is not complete after building"
    ok "human-read filter index"
}

# -- check -------------------------------------------------------------------
done_check() { return 1; }   # always runs
run_check() {
    local problems=0
    [ -x "$REPO/bin/openbiota" ] || { warn "bin/openbiota is not executable"; problems=1; }
    run_logged "bin/openbiota --version" "$REPO/bin/openbiota" --version || { warn "bin/openbiota does not start"; problems=1; }
    run_logged "openbiota doctor" "$VENV/bin/openbiota" doctor || { warn "openbiota doctor reported a problem (see install.log)"; problems=1; }
    run_logged "importing the report modules" "$VENV/bin/python" -c \
        'import openbiota.cli, openbiota.pdfreport, openbiota.solutions, openbiota.engines.metaphlan' \
        || { warn "a python import failed (see install.log)"; problems=1; }
    local name
    for name in python profiler data metaphlan_db diamond_db models; do
        ready "$name" || { warn "$name is missing"; problems=1; }
    done
    if [ $SKIP_HOST_INDEX -eq 0 ]; then ready host_index || { warn "host index missing"; problems=1; }; fi
    "$VENV/bin/python" -m openbiota.installcheck --all >>"$LOG" 2>&1 || true
    if [ $problems -ne 0 ]; then
        fail_with_log "the installation is not complete; the warnings above say what is missing"
    fi
    ok "every piece is in place"
}

# -- the --full steps: each is a Makefile target --------------------------------
make_step() {
    # make_step TARGET  - run a Makefile target as a step, logged
    have make || ensure_tool make || die "--full needs make"
    run_logged "make $1" make --no-print-directory "$1" THREADS="$THREADS" || fail_with_log "make $1 failed"
}
done_profiler4()      { [ -x "$REPO/.venv-mpa4/bin/metaphlan" ]; }
run_profiler4()       { make_step profilers; ok "MetaPhlAn 4 in .venv-mpa4"; }
done_vendor_tools()   { [ -x "$REPO/vendor/bin/blastn" ] || [ -n "$(ls -d "$REPO"/vendor/ncbi-blast-*/bin/blastn 2>/dev/null)" ]; }
run_vendor_tools()    { make_step vendor-tools; ok "vendored BLAST+, AMRFinderPlus and StxTyper"; }
done_strain_tools()   { [ -x "$REPO/.venv-strain/bin/ectyper" ]; }
run_strain_tools()    { make_step strain-tools; ok "strain and typing tools"; }
done_refs_strain()    { [ -s "$REPO/refs/strain_refs.lock.json" ] && [ -d "$REPO/refs/panels/probiotics" ]; }
run_refs_strain()     { make_step refs-strain; ok "strain reference panels"; }
done_pathogens()      { [ -s "$REPO/refs/pathogens/bundle/manifest.json" ]; }
run_pathogens() {
    note "377 reference genomes from NCBI, then k-mer ownership and a competitive index; this is the heaviest build"
    run_logged "building the pathogen bundle" "$VENV/bin/python" tools/build_pathogen_bundle.py --threads "$THREADS" \
        || fail_with_log "tools/build_pathogen_bundle.py failed"
    ok "pathogen reference bundle"
}
done_substrate_refs() { [ -s "$REPO/refs/dbcan/snapshot.json" ]; }
run_substrate_refs()  { make_step substrate-refs; ok "dbCAN substrate references"; }
done_tools_expanded() { [ -x "$REPO/.venv-mpa42/bin/metaphlan" ] && [ -x "$REPO/.venv-motus/bin/motus" ] && [ -x "$REPO/.venv-singlem/bin/singlem" ]; }
run_tools_expanded()  { make_step tools-expanded; ok "expanded detection tools"; }
done_refs_expanded()  { [ -d "$REPO/refs/motus/db" ] && [ -d "$REPO/refs/kraken2" ] && [ -d "$REPO/refs/metaphlan4_db" ]; }
run_refs_expanded()   { make_step refs-expanded; ok "expanded references"; }
done_expansion()      { [ -s "$REPO/refs/expanded/reconciliation/summary.json" ]; }
run_expansion()       { make_step expansion; ok "crosswalks and reconciliation"; }
done_genome_lane()    { [ -x "$REPO/vendor/sylph/bin/sylph" ] && [ -n "$(ls "$REPO"/refs/sylph/*.syl2db 2>/dev/null)" ]; }
run_genome_lane()     { make_step genome-lane; ok "sylph and the GTDB sketch"; }
done_mycobiome_refs() { [ -s "$REPO/refs/mycobiome/reference_lock.json" ]; }
run_mycobiome_refs()  { make_step mycobiome-refs; ok "mycobiome references"; }

# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #

all_steps() {
    if [ $FULL -eq 1 ]; then echo "$STANDARD_STEPS $FULL_STEPS"; else echo "$STANDARD_STEPS"; fi
}

in_list() {
    # in_list NEEDLE "a,b,c"
    local n=$1 list=$2 item
    local old_ifs=$IFS; IFS=','
    for item in $list; do
        IFS=$old_ifs
        [ "$item" = "$n" ] && return 0
        IFS=','
    done
    IFS=$old_ifs
    return 1
}

known_step() {
    local s
    for s in $STANDARD_STEPS $FULL_STEPS; do [ "$s" = "$1" ] && return 0; done
    return 1
}

# Validated here, in the main shell, because a `die` inside a command
# substitution only ends the subshell and the script would carry on with
# the error text as its step list.
validate_step_names() {
    local s list=$1
    local old_ifs=$IFS; IFS=','
    for s in $list; do
        IFS=$old_ifs
        known_step "$s" || die "unknown step: $s (see --help)"
        IFS=','
    done
    IFS=$old_ifs
}
[ -z "$ONLY" ] || validate_step_names "$ONLY"
[ -z "$SKIP" ] || validate_step_names "$SKIP"

selected_steps() {
    local s out=""
    if [ -n "$ONLY" ]; then
        local old_ifs=$IFS; IFS=','
        for s in $ONLY; do
            IFS=$old_ifs
            out="$out $s"
            IFS=','
        done
        IFS=$old_ifs
        echo "$out"
        return
    fi
    for s in $(all_steps); do
        if [ -n "$SKIP" ] && in_list "$s" "$SKIP"; then continue; fi
        out="$out $s"
    done
    echo "$out"
}

# Only one installer at a time: two of them building the same index would
# corrupt it, and that is easy to do by accident from two terminals.
if [ -e "$LOCK" ] && [ $LIST_ONLY -eq 0 ]; then
    other=$(cat "$LOCK" 2>/dev/null || echo "?")
    if [ -n "$other" ] && kill -0 "$other" 2>/dev/null; then
        die "another ./install.sh (pid $other) is running; wait for it, or remove .install.lock if it is stale"
    fi
    rm -f "$LOCK"
fi
if [ $LIST_ONLY -eq 0 ]; then
    echo $$ >"$LOCK"
    cleanup() {
        rm -f "$LOCK"
        [ $TTY -eq 1 ] && printf '\033[2K\r'
    }
    trap cleanup EXIT
    trap 'printf "\n"; warn "interrupted; run ./install.sh again to resume"; exit 130' INT TERM
fi

if [ $LIST_ONLY -eq 1 ]; then
    printf '%sInstall steps%s (%s)\n\n' "$BOLD" "$RESET" "$([ $FULL -eq 1 ] && echo full || echo standard)"
    for s in $(all_steps); do
        if [ "$s" = preflight ] || [ "$s" = check ]; then state="  runs every time"
        elif "done_$s" 2>/dev/null; then state="${GREEN}  done${RESET}"
        else state="${YELLOW}  to do${RESET}"; fi
        printf '  %-16s%s\n' "$s" "$state"
    done
    exit 0
fi

: >>"$LOG"
log "===== install.sh started: $(date) in $REPO ($OS $ARCH) full=$FULL skip_host=$SKIP_HOST_INDEX threads=$THREADS"

STEPS=$(selected_steps)
STEP_TOTAL=0
for step in $STEPS; do STEP_TOTAL=$(( STEP_TOTAL + 1 )); done
T_START=$(date +%s)

printf '%sOpenBiota installer%s  %s%s%s\n' "$BOLD" "$RESET" "$DIM" "$([ $FULL -eq 1 ] && echo "full install" || echo "standard install")" "$RESET"
printf '%sLog: %s%s\n' "$DIM" "install.log" "$RESET"

for step in $STEPS; do
    label=$step
    case "$step" in
        preflight) label="Checking this machine" ;;
        system) label="System tools" ;;
        python) label="Python environment" ;;
        launcher) label="The openbiota command" ;;
        profiler) label="MetaPhlAn 3 profiler" ;;
        data) label="Reference data from the repository" ;;
        metaphlan_db) label="MetaPhlAn marker database" ;;
        diamond_db) label="Protein search database" ;;
        models) label="Metabolic model library" ;;
        host_index) label="Human-read filter index" ;;
        check) label="Final check" ;;
        profiler4) label="MetaPhlAn 4 (full)" ;;
        vendor_tools) label="Vendored tools (full)" ;;
        strain_tools) label="Strain and typing tools (full)" ;;
        refs_strain) label="Strain references (full)" ;;
        pathogens) label="Pathogen reference bundle (full)" ;;
        substrate_refs) label="Substrate references (full)" ;;
        tools_expanded) label="Expanded detection tools (full)" ;;
        refs_expanded) label="Expanded references (full)" ;;
        expansion) label="Crosswalks and reconciliation (full)" ;;
        genome_lane) label="Whole-genome lane (full)" ;;
        mycobiome_refs) label="Mycobiome references (full)" ;;
    esac
    begin_step "$label"
    if [ $FORCE -eq 0 ] && "done_$step"; then
        ok "already done"
        continue
    fi
    "run_$step"
    end_step
done

total=$(( $(date +%s) - T_START ))
printf '\n%s%s%s %s%d/%d%s\n' "$CYAN" "$(bar "$STEP_TOTAL" "$STEP_TOTAL")" "$RESET" "$DIM" "$STEP_TOTAL" "$STEP_TOTAL" "$RESET"
printf '\n%s%sInstalled%s in %s.\n' "$GREEN" "$BOLD" "$RESET" "$(elapsed_text $total)"
cat <<EOF

Two steps left, and they are yours:

  1. Copy your FASTQ files into a folder under ${BOLD}fastq/${RESET}, one folder per sample:

       fastq/my-sample/my-sample_1.fastq.gz
       fastq/my-sample/my-sample_2.fastq.gz

     The files are usually named something like NAME_1.fastq.gz and
     NAME_2.fastq.gz (or NAME_R1 / NAME_R2). The sample is named by the file.

  2. Run the screen, from this folder:

       ${BOLD}bin/openbiota run${RESET}

     With several samples in fastq/, name one: ${BOLD}bin/openbiota run --sample my-sample${RESET}
     The report is written to ${BOLD}results/my-sample/my-sample_report.pdf${RESET}
EOF
if [ $NEW_TERMINAL_NOTE -eq 1 ]; then
    printf '\n     In any %snew%s terminal window the same command is just %sopenbiota run%s.\n' "$BOLD" "$RESET" "$BOLD" "$RESET"
elif [ $NO_PATH -eq 0 ] && have openbiota; then
    printf '\n     %sopenbiota run%s (without bin/) works too.\n' "$BOLD" "$RESET"
fi
printf '\n'
