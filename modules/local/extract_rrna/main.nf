process EXTRACT_RRNA {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/emboss:6.6.0--h0f19ade_14'

    input:
    tuple val(meta), path(sequence), path(ssu_tab), path(lsu_tab)

    output:
    tuple val(meta), path("${meta.id}_SSU.fasta"), emit: ssu_fasta, optional: true
    tuple val(meta), path("${meta.id}_LSU.fasta"), emit: lsu_fasta, optional: true
    tuple val(meta), path("${meta.id}_ITS.fasta"), emit: its_fasta, optional: true
    tuple val(meta), path("${meta.id}_rrna.tsv"),  emit: rrna_log, optional: true
    path "versions.yml",                           emit: versions

    shell = ['/bin/bash', '-euo', 'pipefail']

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    write_versions() {
        local emboss_ver
        emboss_ver=\$(seqret -version 2>&1 | grep -oE 'EMBOSS:[[:space:]]*[0-9.]+' | sed 's/EMBOSS:[[:space:]]*//' || echo "unknown")
        printf '"${task.process}":\\n    emboss: %s\\n' "\$emboss_ver" > versions.yml
    }

    skip_sample() {
        echo "WARNING: \$1 -- skipping sample ${meta.id}" >&2
        write_versions
        exit 0
    }

    # -------------------------------------------------------------------------
    # 0. Handle gzip-compressed input.
    # -------------------------------------------------------------------------
    WORK_FASTA="${sequence}"
    if od -A n -t x1 -N 2 "${sequence}" | grep -q "1f 8b"; then
        gunzip -c "${sequence}" > decompressed_input.fasta
        WORK_FASTA="decompressed_input.fasta"
    fi

    # -------------------------------------------------------------------------
    # 1. parse_tab / parse_normalised helpers
    # -------------------------------------------------------------------------
    parse_tab() {
        awk '!/^#/ && NF > 0 { print \$1, \$8, \$9, \$10, \$15 }' "\$1"
    }

    parse_normalised() {
        local tabfile=\$1
        awk '!/^#/ && NF > 0 {
            start = (\$8 < \$9) ? \$8 : \$9
            end   = (\$8 < \$9) ? \$9 : \$8
            print \$1, start, end, \$10, \$15
        }' "\$tabfile"
    }

    # -------------------------------------------------------------------------
    # 1b. Check both tab files have at least one real hit.
    # -------------------------------------------------------------------------
    SSU_COUNT=\$(parse_tab "${ssu_tab}" | wc -l)
    LSU_COUNT=\$(parse_tab "${lsu_tab}" | wc -l)

    [ "\$SSU_COUNT" -eq 0 ] && skip_sample "No SSU (18S) hits in ${ssu_tab}"
    [ "\$LSU_COUNT" -eq 0 ] && skip_sample "No LSU (28S) hits in ${lsu_tab}"

    ALL_SSU=\$(parse_normalised "${ssu_tab}")
    ALL_LSU=\$(parse_normalised "${lsu_tab}")

    [ -z "\$ALL_SSU" ] && skip_sample "No SSU (18S) hits in ${ssu_tab}"
    [ -z "\$ALL_LSU" ] && skip_sample "No LSU (28S) hits in ${lsu_tab}"

    # -------------------------------------------------------------------------
    # 3. Pair formation, filtering, and selection.
    # -------------------------------------------------------------------------
    ALL_PAIRS=\$(awk -v sample="${meta.id}" '
    BEGIN {
        MAX_SPAN = 10000
        OFS      = "\\t"
    }

    NR == FNR {
        ssu_contig[NR] = \$1
        ssu_start[NR]  = \$2
        ssu_end[NR]    = \$3
        ssu_strand[NR] = \$4
        ssu_score[NR]  = \$5
        n_ssu          = NR
        next
    }

    {
        idx = NR - n_ssu
        lsu_contig[idx] = \$1
        lsu_start[idx]  = \$2
        lsu_end[idx]    = \$3
        lsu_strand[idx] = \$4
        lsu_score[idx]  = \$5
        n_lsu           = idx
    }

    END {
        n_pairs = 0

        for (i = 1; i <= n_ssu; i++) {

            closest_dist = -1
            closest_j    = -1

            for (j = 1; j <= n_lsu; j++) {
                if (ssu_contig[i] != lsu_contig[j]) continue

                # LSU must come AFTER SSU (biologically: SSU-ITS-LSU order)
                if (ssu_strand[i] == "+") {
                    if (lsu_start[j] <= ssu_end[i]) continue
                    dist = lsu_start[j] - ssu_end[i]
                } else {
                    if (lsu_end[j] >= ssu_start[i]) continue
                    dist = ssu_start[i] - lsu_end[j]
                }

                if (closest_dist < 0 || dist < closest_dist) {
                    closest_dist = dist
                    closest_j    = j
                }
            }

            if (closest_j < 0) continue

            j = closest_j

            locus_start = (ssu_start[i] < lsu_start[j]) ? ssu_start[i] : lsu_start[j]
            locus_end   = (ssu_end[i]   > lsu_end[j])   ? ssu_end[i]   : lsu_end[j]
            span        = locus_end - locus_start + 1

            if (span > MAX_SPAN) continue

            if (ssu_end[i] < lsu_start[j]) {
                its_start = ssu_end[i]   + 1
                its_end   = lsu_start[j] - 1
            } else {
                its_start = lsu_end[j]   + 1
                its_end   = ssu_start[i] - 1
            }

            if (its_end < its_start) continue

            combined    = ssu_score[i] + lsu_score[j]
            ssu_len     = ssu_end[i]   - ssu_start[i] + 1
            lsu_len     = lsu_end[j]   - lsu_start[j] + 1
            its_len     = its_end      - its_start     + 1

            if (ssu_strand[i] == "+") {
                strand_span = lsu_end[j]   - ssu_start[i] + 1
            } else {
                strand_span = ssu_end[i]   - lsu_start[j] + 1
            }

            n_pairs++
            pair_contig[n_pairs]      = ssu_contig[i]
            pair_ssu_start[n_pairs]   = ssu_start[i]
            pair_ssu_end[n_pairs]     = ssu_end[i]
            pair_ssu_strand[n_pairs]  = ssu_strand[i]
            pair_ssu_len[n_pairs]     = ssu_len
            pair_lsu_start[n_pairs]   = lsu_start[j]
            pair_lsu_end[n_pairs]     = lsu_end[j]
            pair_lsu_strand[n_pairs]  = lsu_strand[j]
            pair_lsu_len[n_pairs]     = lsu_len
            pair_its_start[n_pairs]   = its_start
            pair_its_end[n_pairs]     = its_end
            pair_its_strand[n_pairs]  = ssu_strand[i]
            pair_its_len[n_pairs]     = its_len
            pair_span[n_pairs]        = span
            pair_score[n_pairs]       = combined
            pair_strand_span[n_pairs] = strand_span
        }

        if (n_pairs == 0) exit 1

        best_score = -1
        for (k = 1; k <= n_pairs; k++) {
            if (pair_score[k] > best_score) best_score = pair_score[k]
        }

        best_strand_span = -1
        best_k           = -1
        for (k = 1; k <= n_pairs; k++) {
            if (pair_score[k] != best_score) continue
            if (pair_strand_span[k] > best_strand_span) {
                best_strand_span = pair_strand_span[k]
                best_k           = k
            }
        }

        k = best_k
        print "true", sample, \
            pair_contig[k], \
            pair_ssu_start[k], pair_ssu_end[k], pair_ssu_strand[k], pair_ssu_len[k], \
            pair_lsu_start[k], pair_lsu_end[k], pair_lsu_strand[k], pair_lsu_len[k], \
            pair_its_start[k], pair_its_end[k], pair_its_strand[k], pair_its_len[k], \
            pair_span[k], pair_score[k], pair_strand_span[k]

        for (k = 1; k <= n_pairs; k++) {
            if (k == best_k) continue
            print "false", sample, \
                pair_contig[k], \
                pair_ssu_start[k], pair_ssu_end[k], pair_ssu_strand[k], pair_ssu_len[k], \
                pair_lsu_start[k], pair_lsu_end[k], pair_lsu_strand[k], pair_lsu_len[k], \
                pair_its_start[k], pair_its_end[k], pair_its_strand[k], pair_its_len[k], \
                pair_span[k], pair_score[k], pair_strand_span[k]
        }
    }
    ' <(echo "\$ALL_SSU") <(echo "\$ALL_LSU")) || skip_sample "No valid SSU/LSU pair found within 10kb span"

    # -------------------------------------------------------------------------
    # 4. Write TSV (header + all pairs)
    # -------------------------------------------------------------------------
    printf 'is_best\tsample\tcontig\tssu_start\tssu_end\tssu_strand\tssu_length\tlsu_start\tlsu_end\tlsu_strand\tlsu_length\tits_start\tits_end\tits_strand\tits_length\tlocus_span_bp\tcombined_score\tstrand_aware_span_bp\n' \
        > "${prefix}_rrna.tsv"

    echo "\$ALL_PAIRS" >> "${prefix}_rrna.tsv"

    # -------------------------------------------------------------------------
    # 5. Extract best pair coordinates from first row of ALL_PAIRS
    # -------------------------------------------------------------------------
    BEST_ROW=\$(echo "\$ALL_PAIRS" | head -1)

    CHOSEN_CONTIG=\$(echo "\$BEST_ROW" | awk '{ print \$3 }')
    SSU_START=\$(    echo "\$BEST_ROW" | awk '{ print \$4 }')
    SSU_END=\$(      echo "\$BEST_ROW" | awk '{ print \$5 }')
    SSU_STRAND=\$(   echo "\$BEST_ROW" | awk '{ print \$6 }')
    LSU_START=\$(    echo "\$BEST_ROW" | awk '{ print \$8 }')
    LSU_END=\$(      echo "\$BEST_ROW" | awk '{ print \$9 }')
    ITS_START=\$(    echo "\$BEST_ROW" | awk '{ print \$12 }')
    ITS_END=\$(      echo "\$BEST_ROW" | awk '{ print \$13 }')
    LOCUS_STRAND="\$SSU_STRAND"
    LOCUS_SPAN=\$(   echo "\$BEST_ROW" | awk '{ print \$16 }')

    echo "Selected contig : \$CHOSEN_CONTIG"
    echo "Best SSU        : \$CHOSEN_CONTIG:\$SSU_START-\$SSU_END (\$SSU_STRAND)"
    echo "Best LSU        : \$CHOSEN_CONTIG:\$LSU_START-\$LSU_END"
    echo "ITS region      : \$CHOSEN_CONTIG:\$ITS_START-\$ITS_END"
    echo "Locus span      : \$LOCUS_SPAN bp"

    # -------------------------------------------------------------------------
    # 6. Extract sequences with EMBOSS seqret + descseq.
    # -------------------------------------------------------------------------
    STRAND_FLAG=""
    [ "\$LOCUS_STRAND" = "-" ] && STRAND_FLAG="r"

    extractSequence() {
        local name=\$1
        local start=\$2
        local end=\$3
        local flag=\$4
        local outfile=\$5
        local coords="\${CHOSEN_CONTIG}[\${start}:\${end}:\${flag}]"
        seqret "\${WORK_FASTA}:\${coords}" -filter \\
            | descseq -filter \\
                -name "\$name" \\
                -desc "\${CHOSEN_CONTIG}:\${start}-\${end}(\${LOCUS_STRAND})" \\
            > "\$outfile"
    }

    extractSequence "${prefix}_SSU" "\$SSU_START" "\$SSU_END" "\$STRAND_FLAG" "${prefix}_SSU.fasta"
    extractSequence "${prefix}_ITS" "\$ITS_START" "\$ITS_END" "\$STRAND_FLAG" "${prefix}_ITS.fasta"
    extractSequence "${prefix}_LSU" "\$LSU_START" "\$LSU_END" "\$STRAND_FLAG" "${prefix}_LSU.fasta"

    # -------------------------------------------------------------------------
    # 7. Versions
    # -------------------------------------------------------------------------
    write_versions
    """
}