#!/bin/bash
# =============================================================================
# HUMAnN3 post-processing pipeline (v2 - handles missing utility_mapping)
# -----------------------------------------------------------------------------
# Auto-detects whether the full utility_mapping database is installed.
#  - If yes: regroups UniRef90 -> KO and EC (uniref90_ko, uniref90_level4ec)
#  - If no:  falls back to uniref90_rxn (MetaCyc reactions, ships by default)
#            and warns the user how to get KO/EC.
# =============================================================================

set -euo pipefail

# ---- EDIT THESE PATHS --------------------------------------------------------
INPUT_DIR="humann3_output"        # per-sample HUMAnN3 output directories
WORK_DIR="humann3_processed"      # all derived tables land here
# -----------------------------------------------------------------------------

mkdir -p "$WORK_DIR"

# ---- Detect available regroup options ---------------------------------------
REGROUP_OPTS=$(humann_regroup_table -h 2>&1 | grep -oE 'uniref90_[a-z0-9]+' | sort -u || true)
HAS_KO=0
HAS_EC=0
echo "$REGROUP_OPTS" | grep -qx "uniref90_ko"        && HAS_KO=1
echo "$REGROUP_OPTS" | grep -qx "uniref90_level4ec"  && HAS_EC=1

if [[ $HAS_KO -eq 1 && $HAS_EC -eq 1 ]]; then
  echo ">>> Full utility_mapping detected (KO + EC available)."
else
  cat <<'EOF' >&2
>>> WARNING: full utility_mapping is NOT installed.
    Only uniref90_rxn (MetaCyc reactions) will be produced.
    To enable KO / EC / Pfam / eggNOG / GO regrouping, run:

      humann_databases --download utility_mapping full <PATH>

    then re-run this script.
EOF
fi

# ---- 1. Join ----------------------------------------------------------------
echo "[1/5] Joining per-sample tables ..."
humann_join_tables -i "$INPUT_DIR" -o "$WORK_DIR/all_genefamilies.tsv"  --file_name genefamilies  -s
humann_join_tables -i "$INPUT_DIR" -o "$WORK_DIR/all_pathabundance.tsv" --file_name pathabundance -s
humann_join_tables -i "$INPUT_DIR" -o "$WORK_DIR/all_pathcoverage.tsv"  --file_name pathcoverage  -s

# ---- 2. Normalize -----------------------------------------------------------
echo "[2/5] Normalizing to CPM ..."
humann_renorm_table -i "$WORK_DIR/all_pathabundance.tsv" \
                    -o "$WORK_DIR/all_pathabundance_cpm.tsv" \
                    --units cpm --update-snames
humann_renorm_table -i "$WORK_DIR/all_genefamilies.tsv" \
                    -o "$WORK_DIR/all_genefamilies_cpm.tsv" \
                    --units cpm --update-snames

# ---- 3. Regroup (KO + EC if available, otherwise RXN) -----------------------
echo "[3/5] Regrouping gene families ..."

if [[ $HAS_KO -eq 1 ]]; then
  humann_regroup_table -i "$WORK_DIR/all_genefamilies.tsv" \
                       -o "$WORK_DIR/all_ko.tsv" --groups uniref90_ko
  humann_renorm_table  -i "$WORK_DIR/all_ko.tsv" \
                       -o "$WORK_DIR/all_ko_cpm.tsv" --units cpm --update-snames
fi

if [[ $HAS_EC -eq 1 ]]; then
  humann_regroup_table -i "$WORK_DIR/all_genefamilies.tsv" \
                       -o "$WORK_DIR/all_ec.tsv" --groups uniref90_level4ec
  humann_renorm_table  -i "$WORK_DIR/all_ec.tsv" \
                       -o "$WORK_DIR/all_ec_cpm.tsv" --units cpm --update-snames
fi

# Always also produce the rxn regrouping (built-in; useful as MetaCyc QC layer)
humann_regroup_table -i "$WORK_DIR/all_genefamilies.tsv" \
                     -o "$WORK_DIR/all_rxn.tsv" --groups uniref90_rxn
humann_renorm_table  -i "$WORK_DIR/all_rxn.tsv" \
                     -o "$WORK_DIR/all_rxn_cpm.tsv" --units cpm --update-snames

# ---- 4. Rename --------------------------------------------------------------
echo "[4/5] Adding human-readable names ..."
humann_rename_table -i "$WORK_DIR/all_pathabundance_cpm.tsv" \
                    -o "$WORK_DIR/all_pathabundance_cpm_named.tsv" --names metacyc-pwy

if [[ $HAS_KO -eq 1 ]]; then
  humann_rename_table -i "$WORK_DIR/all_ko_cpm.tsv" \
                      -o "$WORK_DIR/all_ko_cpm_named.tsv" --names kegg-orthology
fi
if [[ $HAS_EC -eq 1 ]]; then
  humann_rename_table -i "$WORK_DIR/all_ec_cpm.tsv" \
                      -o "$WORK_DIR/all_ec_cpm_named.tsv" --names ec
fi
humann_rename_table -i "$WORK_DIR/all_rxn_cpm.tsv" \
                    -o "$WORK_DIR/all_rxn_cpm_named.tsv" --names metacyc-rxn

# ---- 5. Split stratified / unstratified -------------------------------------
echo "[5/5] Splitting stratified vs unstratified ..."

humann_split_stratified_table -i "$WORK_DIR/all_pathabundance_cpm_named.tsv" -o "$WORK_DIR/path_split_cpm/"
humann_split_stratified_table -i "$WORK_DIR/all_pathabundance.tsv"           -o "$WORK_DIR/path_split_rpk/"
humann_split_stratified_table -i "$WORK_DIR/all_rxn_cpm_named.tsv"           -o "$WORK_DIR/rxn_split_cpm/"
humann_split_stratified_table -i "$WORK_DIR/all_rxn.tsv"                     -o "$WORK_DIR/rxn_split_rpk/"

if [[ $HAS_KO -eq 1 ]]; then
  humann_split_stratified_table -i "$WORK_DIR/all_ko_cpm_named.tsv" -o "$WORK_DIR/ko_split_cpm/"
  humann_split_stratified_table -i "$WORK_DIR/all_ko.tsv"           -o "$WORK_DIR/ko_split_rpk/"
fi
if [[ $HAS_EC -eq 1 ]]; then
  humann_split_stratified_table -i "$WORK_DIR/all_ec_cpm_named.tsv" -o "$WORK_DIR/ec_split_cpm/"
  humann_split_stratified_table -i "$WORK_DIR/all_ec.tsv"           -o "$WORK_DIR/ec_split_rpk/"
fi

cat <<EOF

Post-processing complete. Outputs in: $WORK_DIR
  Pathways (CPM, named, unstratified): path_split_cpm/all_pathabundance_cpm_named_unstratified.tsv
  Pathways (RPK, unstratified):        path_split_rpk/all_pathabundance_unstratified.tsv
  Pathways (CPM, stratified):          path_split_cpm/all_pathabundance_cpm_named_stratified.tsv
  MetaCyc reactions (CPM, named):      rxn_split_cpm/all_rxn_cpm_named_unstratified.tsv
EOF

if [[ $HAS_KO -eq 1 ]]; then
  echo "  KOs (CPM, named):                    ko_split_cpm/all_ko_cpm_named_unstratified.tsv"
  echo "  KOs (RPK):                           ko_split_rpk/all_ko_unstratified.tsv"
fi
if [[ $HAS_EC -eq 1 ]]; then
  echo "  ECs (CPM, named):                    ec_split_cpm/all_ec_cpm_named_unstratified.tsv"
fi
