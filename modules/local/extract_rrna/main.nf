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
        emboss_ver=\$(seqret -version 2>&1 | grep -oP 'EMBOSS:\\s*\\K[\\d.]+' || echo "unknown")
        printf '"${task.process}":\\n    emboss: %s\\n' "\$emboss_ver" > versions.yml
    }

    skip_sample() {
        echo "WARNING: \$1 -- skipping sample ${meta.id}" >&2
        write_versions
        exit 0
    }

    # -------------------------------------------------------------------------
    # 0. Handle gzip-compressed input.
    #    We use 'od' to read the magic bytes instead of 'file', which is not
    #    present in the EMBOSS container. Gzip files start with bytes 1f 8b.
    # -------------------------------------------------------------------------
    WORK_FASTA="${sequence}"
    if od -A n -t x1 -N 2 "${sequence}" | grep -q "1f 8b"; then
        gunzip -c "${sequence}" > decompressed_input.fasta
        WORK_FASTA="decompressed_input.fasta"
    fi

    # -------------------------------------------------------------------------
    # 1. parse_tab: strip comment lines and emit:
    #    contig  seq_from  seq_to  strand  score
    # -------------------------------------------------------------------------
    parse_tab() {
        awk '!/^#/ && NF > 0 { print \$1, \$8, \$9, \$10, \$16 }' "\$1"
    }

    # -------------------------------------------------------------------------
    # 1b. Check both tab files have at least one real hit.
    # -------------------------------------------------------------------------
    SSU_COUNT=\$(parse_tab "${ssu_tab}" | wc -l)
    LSU_COUNT=\$(parse_tab "${lsu_tab}" | wc -l)

    [ "\$SSU_COUNT" -eq 0 ] && skip_sample "No SSU (18S) hits in ${ssu_tab}"
    [ "\$LSU_COUNT" -eq 0 ] && skip_sample "No LSU (28S) hits in ${lsu_tab}"

# -------------------------------------------------------------------------
    # 2. Parse all hits from both files into two flat lists.
    #    Fields: contig  start  end  strand  score
    #    Coordinates are normalised (start < end) as before.
    # -------------------------------------------------------------------------
    parse_normalised() {
        local tabfile=\$1
        awk '!/^#/ && NF > 0 {
            start = (\$8 < \$9) ? \$8 : \$9
            end   = (\$8 < \$9) ? \$9 : \$8
            print \$1, start, end, \$10, \$16
        }' "\$tabfile"
    }

    ALL_SSU=\$(parse_normalised "${ssu_tab}")
    ALL_LSU=\$(parse_normalised "${lsu_tab}")

    [ -z "\$ALL_SSU" ] && skip_sample "No SSU (18S) hits in ${ssu_tab}"
    [ -z "\$ALL_LSU" ] && skip_sample "No LSU (28S) hits in ${lsu_tab}"

    # -------------------------------------------------------------------------
    # 3. For each SSU hit, find its closest LSU hit on the same contig.
    #    This defines the "pair" — the two features most likely from the same
    #    rDNA repeat unit.
    #
    #    Then apply a span filter: the total region from the start of 18S to
    #    the end of 28S must be <= 10,000 bp. A typical rDNA unit is ~8kb so
    #    anything larger is likely a spurious cross-repeat pairing.
    #
    #    Among all pairs that pass the filter, pick the one with the highest
    #    combined SSU+LSU bit score. Ties are broken by SSU file order.
    # -------------------------------------------------------------------------
    BEST_PAIR=\$(awk '
    BEGIN {
        best_score = -1
        MAX_SPAN   = 10000
    }

    # First file: load all SSU hits
    NR == FNR {
        ssu_contig[NR] = \$1
        ssu_start[NR]  = \$2
        ssu_end[NR]    = \$3
        ssu_strand[NR] = \$4
        ssu_score[NR]  = \$5
        n_ssu          = NR
        next
    }

    # Second file: load all LSU hits into an array
    {
        lsu_contig[NR - n_ssu]  = \$1
        lsu_start[NR - n_ssu]   = \$2
        lsu_end[NR - n_ssu]     = \$3
        lsu_strand[NR - n_ssu]  = \$4
        lsu_score[NR - n_ssu]   = \$5
        n_lsu                   = NR - n_ssu
    }

    END {
        # For every SSU hit, find the closest LSU hit on the same contig
        for (i = 1; i <= n_ssu; i++) {

            closest_dist  = -1
            closest_j     = -1

            for (j = 1; j <= n_lsu; j++) {

                if (ssu_contig[i] != lsu_contig[j]) continue

                # Genomic gap between the two features
                if      (lsu_start[j] > ssu_end[i])    dist = lsu_start[j] - ssu_end[i]
                else if (ssu_start[i] > lsu_end[j])    dist = ssu_start[i] - lsu_end[j]
                else                                    dist = 0

                if (closest_dist < 0 || dist < closest_dist) {
                    closest_dist = dist
                    closest_j    = j
                }
            }

            # No LSU found on the same contig for this SSU — skip
            if (closest_j < 0) continue

            j = closest_j

            # Compute total span of the rDNA unit
            locus_start = (ssu_start[i] < lsu_start[j]) ? ssu_start[i] : lsu_start[j]
            locus_end   = (ssu_end[i]   > lsu_end[j])   ? ssu_end[i]   : lsu_end[j]
            span        = locus_end - locus_start + 1

            # Apply span filter
            if (span > MAX_SPAN) continue

            combined = ssu_score[i] + lsu_score[j]

            # Keep the best pair using the tiebreaker chain:
            #   1. highest combined score
            #   2. if tied: longest span
            #   3. if still tied: SSU file order (implicit, i iterates in order)
            if (best_score < 0 || \
                combined > best_score || \
                (combined == best_score && span > best_span)) {
                best_score      = combined
                best_span       = span
                best_contig     = ssu_contig[i]
                best_ssu_start  = ssu_start[i]
                best_ssu_end    = ssu_end[i]
                best_ssu_strand = ssu_strand[i]
                best_lsu_start  = lsu_start[j]
                best_lsu_end    = lsu_end[j]
                best_lsu_strand = lsu_strand[j]
            }
        }

        if (best_score < 0) exit 1

        print best_contig, \
              best_ssu_start, best_ssu_end, best_ssu_strand, \
              best_lsu_start, best_lsu_end, best_lsu_strand, \
              best_span
    }
    ' <(echo "\$ALL_SSU") <(echo "\$ALL_LSU")) || skip_sample "No valid SSU/LSU pair found within 10kb span"

    CHOSEN_CONTIG=\$(echo "\$BEST_PAIR" | awk '{ print \$1 }')
    SSU_START=\$(    echo "\$BEST_PAIR" | awk '{ print \$2 }')
    SSU_END=\$(      echo "\$BEST_PAIR" | awk '{ print \$3 }')
    SSU_STRAND=\$(   echo "\$BEST_PAIR" | awk '{ print \$4 }')
    LSU_START=\$(    echo "\$BEST_PAIR" | awk '{ print \$5 }')
    LSU_END=\$(      echo "\$BEST_PAIR" | awk '{ print \$6 }')
    LSU_STRAND=\$(   echo "\$BEST_PAIR" | awk '{ print \$7 }')
    LOCUS_SPAN=\$(   echo "\$BEST_PAIR" | awk '{ print \$8 }')

    echo "Selected contig : \$CHOSEN_CONTIG"
    echo "Best SSU        : \$CHOSEN_CONTIG:\$SSU_START-\$SSU_END (\$SSU_STRAND)"
    echo "Best LSU        : \$CHOSEN_CONTIG:\$LSU_START-\$LSU_END (\$LSU_STRAND)"
    echo "Locus span      : \$LOCUS_SPAN bp"

    # -------------------------------------------------------------------------
    # 5. Compute ITS coordinates.
    # -------------------------------------------------------------------------
    if [ "\$SSU_END" -lt "\$LSU_START" ]; then
        ITS_START=\$(( SSU_END   + 1 ))
        ITS_END=\$(( LSU_START   - 1 ))
        LOCUS_STRAND="\$SSU_STRAND"
    else
        ITS_START=\$(( LSU_END   + 1 ))
        ITS_END=\$(( SSU_START   - 1 ))
        LOCUS_STRAND="\$SSU_STRAND"
    fi

    if [ "\$ITS_END" -lt "\$ITS_START" ]; then
        skip_sample "ITS coordinates invalid (start=\$ITS_START end=\$ITS_END) -- SSU and LSU may overlap"
    fi

    echo "ITS region: \$CHOSEN_CONTIG:\$ITS_START-\$ITS_END"

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
    # 8. Write extraction log as a TSV file.
    #    Columns: sample, contig, ssu_start, ssu_end, ssu_strand,
    #             lsu_start, lsu_end, lsu_strand, its_start, its_end
    # -------------------------------------------------------------------------
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "sample" \
        "contig" \
        "ssu_start" "ssu_end" "ssu_strand" \
        "lsu_start" "lsu_end" "lsu_strand" \
        "its_start" "its_end" "its_strand" \
        "locus_span_bp" \
        > "${prefix}_rrna.tsv"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${meta.id}" \
        "\$CHOSEN_CONTIG" \
        "\$SSU_START" "\$SSU_END" "\$SSU_STRAND" \
        "\$LSU_START" "\$LSU_END" "\$LSU_STRAND" \
        "\$ITS_START" "\$ITS_END" "\$LOCUS_STRAND" \
        "\$LOCUS_SPAN" \
        >> "${prefix}_rrna.tsv"

    # -------------------------------------------------------------------------
    # 9. Versions
    # -------------------------------------------------------------------------
    write_versions
    """
}