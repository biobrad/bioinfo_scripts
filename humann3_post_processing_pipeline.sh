#!/bin/bash
# =============================================================================
# HUMAnN3 post-processing pipeline
# -----------------------------------------------------------------------------
# Inputs:  per-sample HUMAnN3 output directories or files containing
#          *_genefamilies.tsv, *_pathabundance.tsv, *_pathcoverage.tsv
# Outputs: joined, CPM-normalized, named, and split stratified/unstratified
#          tables suitable for downstream R analysis.
#
# Notes:
#   * --update-snames drops the trailing "_Abundance" / "_Abundance-RPKs" so
#     R sample names match metadata cleanly.
#   * Both CPM-normalized and raw RPK outputs are split: CPM for MaAsLin2 and
#     visualization; RPK (rounded to integers) for ALDEx2 / compositional DA.
# =============================================================================
 
set -euo pipefail
 
# ---- EDIT THESE PATHS --------------------------------------------------------
INPUT_DIR="humann3_output"       # contains per-sample HUMAnN3 outputs
WORK_DIR="humann3_processed"     # all derived tables land here
# -----------------------------------------------------------------------------
 
mkdir -p "$WORK_DIR"
 
echo "[1/5] Joining per-sample tables ..."
humann_join_tables -i "$INPUT_DIR" -o "$WORK_DIR/all_genefamilies.tsv"  --file_name genefamilies  -s
humann_join_tables -i "$INPUT_DIR" -o "$WORK_DIR/all_pathabundance.tsv" --file_name pathabundance -s
humann_join_tables -i "$INPUT_DIR" -o "$WORK_DIR/all_pathcoverage.tsv"  --file_name pathcoverage  -s
 
echo "[2/5] Normalizing to CPM ..."
humann_renorm_table -i "$WORK_DIR/all_pathabundance.tsv" \
                    -o "$WORK_DIR/all_pathabundance_cpm.tsv" \
                    --units cpm --update-snames
humann_renorm_table -i "$WORK_DIR/all_genefamilies.tsv" \
                    -o "$WORK_DIR/all_genefamilies_cpm.tsv" \
                    --units cpm --update-snames
 
echo "[3/5] Regrouping UniRef90 -> KO and EC ..."
humann_regroup_table -i "$WORK_DIR/all_genefamilies.tsv" \
                     -o "$WORK_DIR/all_ko.tsv" --groups uniref90_ko
humann_regroup_table -i "$WORK_DIR/all_genefamilies.tsv" \
                     -o "$WORK_DIR/all_ec.tsv" --groups uniref90_level4ec
 
humann_renorm_table -i "$WORK_DIR/all_ko.tsv" \
                    -o "$WORK_DIR/all_ko_cpm.tsv" --units cpm --update-snames
humann_renorm_table -i "$WORK_DIR/all_ec.tsv" \
                    -o "$WORK_DIR/all_ec_cpm.tsv" --units cpm --update-snames
 
echo "[4/5] Adding human-readable names ..."
humann_rename_table -i "$WORK_DIR/all_pathabundance_cpm.tsv" \
                    -o "$WORK_DIR/all_pathabundance_cpm_named.tsv" --names metacyc-pwy
humann_rename_table -i "$WORK_DIR/all_ko_cpm.tsv" \
                    -o "$WORK_DIR/all_ko_cpm_named.tsv" --names kegg-orthology
humann_rename_table -i "$WORK_DIR/all_ec_cpm.tsv" \
                    -o "$WORK_DIR/all_ec_cpm_named.tsv" --names ec
 
echo "[5/5] Splitting stratified vs unstratified ..."
# CPM-normalized, named (for MaAsLin2 + viz)
humann_split_stratified_table -i "$WORK_DIR/all_pathabundance_cpm_named.tsv" -o "$WORK_DIR/path_split_cpm/"
humann_split_stratified_table -i "$WORK_DIR/all_ko_cpm_named.tsv"            -o "$WORK_DIR/ko_split_cpm/"
 
# Raw RPK (for ALDEx2 / compositional methods that prefer count-like input)
humann_split_stratified_table -i "$WORK_DIR/all_pathabundance.tsv" -o "$WORK_DIR/path_split_rpk/"
humann_split_stratified_table -i "$WORK_DIR/all_ko.tsv"            -o "$WORK_DIR/ko_split_rpk/"
 
cat <<EOF
 
Post-processing complete. Inputs for the R Markdown:
 
  Community-level (CPM, named):
    $WORK_DIR/path_split_cpm/all_pathabundance_cpm_named_unstratified.tsv
    $WORK_DIR/ko_split_cpm/all_ko_cpm_named_unstratified.tsv
 
  Stratified (CPM, named) - for "who's driving this pathway":
    $WORK_DIR/path_split_cpm/all_pathabundance_cpm_named_stratified.tsv
 
  Community-level (RPK, for ALDEx2):
    $WORK_DIR/path_split_rpk/all_pathabundance_unstratified.tsv
    $WORK_DIR/ko_split_rpk/all_ko_unstratified.tsv
 
EOF