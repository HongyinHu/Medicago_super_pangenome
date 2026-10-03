SP=genome_474

cat > header.txt <<'EOF'
TRASH_chr	TRASH_start	TRASH_end	TRASH_id	monomer_len	array_len	peak_chr	peak_start	peak_end	peak_id	peak_score	peak_strand	peak_signal	peak_pvalue	peak_qvalue
EOF

cat header.txt 05_centromere_define/${SP}.TRASH_overlap_CENH3peak.tsv \
  > 05_centromere_define/${SP}.TRASH_overlap_CENH3peak.with_header.tsv
