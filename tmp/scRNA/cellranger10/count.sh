#!/usr/bin/env bash
# One 10x GEX sample; run from the directory that should contain the output ID.
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: bash $0 ID SAMPLE_PREFIX TRANSCRIPTOME FASTQ_DIR" >&2
    echo "Environment: CELLRANGER_BIN=cellranger10 LOCALCORES=20 [LOCALMEM=GB]" >&2
    exit 2
fi

id=$1
sample=$2
reference=$3
fastqs=$4
bin=${CELLRANGER_BIN:-cellranger10}
cores=${LOCALCORES:-20}
[[ $id =~ ^[a-zA-Z0-9_-]+$ ]] || { echo 'Invalid output ID' >&2; exit 2; }
[[ -n $sample && -d $reference && -d $fastqs ]] || {
    echo 'Require a sample prefix and existing reference/FASTQ directories' >&2
    exit 2
}
[[ $cores =~ ^[1-9][0-9]*$ ]] || { echo 'LOCALCORES must be a positive integer' >&2; exit 2; }
if [[ -n ${LOCALMEM:-} && ! $LOCALMEM =~ ^[1-9][0-9]*$ ]]; then
    echo 'LOCALMEM must be a positive integer (GB)' >&2
    exit 2
fi
command -v "$bin" >/dev/null || { echo "Executable not found: $bin" >&2; exit 127; }
version=$("$bin" --version)
printf 'Executable: %s\nVersion: %s\n' "$(command -v "$bin")" "$version"
[[ $version =~ ^cellranger[[:space:]]+(cellranger-)?10\.[0-9]+\.[0-9]+ ]] || {
    echo 'This script requires Cell Ranger 10.x; check CELLRANGER_BIN' >&2
    exit 2
}

cmd=("$bin" count "--id=$id" "--sample=$sample"
    "--transcriptome=$reference" "--fastqs=$fastqs"
    --nosecondary --create-bam=true "--localcores=$cores")
[[ -z ${LOCALMEM:-} ]] || cmd+=("--localmem=$LOCALMEM")
printf 'Started: %s\nCommand:' "$(date -Is)"
printf ' %q' "${cmd[@]}"
printf '\n'
finish() {
    local status=$?
    printf 'Finished: %s; exit_code=%s\n' "$(date -Is)" "$status"
}
trap finish EXIT
"${cmd[@]}"
