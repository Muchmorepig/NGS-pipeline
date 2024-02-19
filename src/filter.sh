#!/bin/bash

filter_bam() {
    local dir=$1
    local th=$2
    local max_jobs=4

    # Check if log directory exists, if not, create it
    if [ ! -d "${dir}/log" ]; then
        echo -e "  Created directory: ${dir}/log"
        mkdir -p "${dir}/log"
    fi

    local job_count=0
    for bam in "${dir}/"*.bam; do
        local base=$(basename "${bam}" .sorted.bam)
        # Run sambamba markdup followed by samtools view and index in a subshell in the background
        (
            sambamba markdup -t "$th" "${bam}" "${dir}/${base}.sorted.markdup.bam" \
                2>"${dir}/log/${base}.sambamba" &&
                samtools view -@ "$th" -bF 1804 -q 20 "${dir}/${base}.sorted.markdup.bam" -o "${dir}/${base}.flt.bam" &&
                samtools index -@ "$th" "${dir}/${base}.flt.bam" &&
                echo >&2 "${base} done" &&
                rm -v "${dir}/${base}.sorted.markdup.bam" "${dir}/${base}.sorted.markdup.bam.bai" &&
                rm -v "${dir}/${base}.sorted.bam" "${dir}/${base}.sorted.bam.bai"
                
        ) &

        # Increment and check if we need to wait for jobs to finish
        ((job_count++))
        if ((job_count >= max_jobs)); then
            wait -n
            ((job_count--))
        fi
    done

    wait
}
