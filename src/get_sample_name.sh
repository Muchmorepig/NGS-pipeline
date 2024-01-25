#!/bin/bash
get_unique_sample_names() {
    local indir=$(realpath "$1")
    local suffix=$2

    find "${indir}" -maxdepth 1 -name "*${suffix}" -printf '%f\n' |
        sed -e "s/${suffix}$//" | sort -u
}
