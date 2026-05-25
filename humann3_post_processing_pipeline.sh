
#!/bin/bash
# =============================================================================
# HUMAnN3 post-processing pipeline (v3)
# -----------------------------------------------------------------------------
# Change from v2: split stratified BEFORE renaming. humann_rename_table
# silently matches nothing when run on stratified output because the row
# keys are "PATHWAY|taxon" and the mapping file only has bare IDs.
# =============================================================================
 
set -euo pipefail
 
# ---- EDIT THESE PATHS --------------------------------------------------------
INPUT_DIR="humann3_output"
WORK_DIR="humann3_processed"
# -----------------------------------------------------------------------------
 
mkdir -p "$WORK_DIR"
 
# ---- Detect available regroup options ---------------------------------------
REGROUP_OPTS=$(humann_regroup_table -h 2>&1 | grep -oE 'uniref90_[a-z0-9]+' | sort -u || true)
HAS_KO=0; HAS_EC=0
echo "$REGROUP_OPTS" | grep -qx "uniref90_ko"       && HAS_KO=1
echo "$REGROUP_OPTS" | grep -qx "uniref90_level4ec" && HAS_EC=1
 
if [[ $HAS_KO -eq 1 && $HAS_EC -eq 1 ]]; then
  echo ">>> Full utility_mapping detected (KO + EC available)."
else
  echo ">>> WARNING: full utility_mapping is NOT installed." >&2
fi
 
# ---- 1. Join ----------------------------------------------------------------
echo "[1/6] Joining per-sample tables ..."
humann_join_tables -i "$INPUT_DIR" -o "$WORK_DIR/all_genefamilies.tsv"  --file_name genefamilies  -s
humann_join_tables -i "$INPUT_DIR" -o "$WORK_DIR/all_pathabundance.tsv" --file_name pathabundance -s
humann_join_tables -i "$INPUT_DIR" -o "$WORK_DIR/all_pathcoverage.tsv"  --file_name pathcoverage  -s
 
# ---- 2. Normalize -----------------------------------------------------------
echo "[2/6] Normalizing to CPM ..."
humann_renorm_table -i "$WORK_DIR/all_pathabundance.tsv" \
                    -o "$WORK_DIR/all_pathabundance_cpm.tsv" \
                    --units cpm --update-snames
humann_renorm_table -i "$WORK_DIR/all_genefamilies.tsv" \
                    -o "$WORK_DIR/all_genefamilies_cpm.tsv" \
                    --units cpm --update-snames
 
# ---- 3. Regroup gene families ----------------------------------------------
echo "[3/6] Regrouping gene families ..."
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
humann_regroup_table -i "$WORK_DIR/all_genefamilies.tsv" \
                     -o "$WORK_DIR/all_rxn.tsv" --groups uniref90_rxn
humann_renorm_table  -i "$WORK_DIR/all_rxn.tsv" \
                     -o "$WORK_DIR/all_rxn_cpm.tsv" --units cpm --update-snames
 
# ---- 4. Split stratified / unstratified FIRST -------------------------------
echo "[4/6] Splitting stratified vs unstratified ..."
humann_split_stratified_table -i "$WORK_DIR/all_pathabundance_cpm.tsv" -o "$WORK_DIR/path_split_cpm/"
humann_split_stratified_table -i "$WORK_DIR/all_pathabundance.tsv"     -o "$WORK_DIR/path_split_rpk/"
humann_split_stratified_table -i "$WORK_DIR/all_rxn_cpm.tsv"           -o "$WORK_DIR/rxn_split_cpm/"
humann_split_stratified_table -i "$WORK_DIR/all_rxn.tsv"               -o "$WORK_DIR/rxn_split_rpk/"
if [[ $HAS_KO -eq 1 ]]; then
  humann_split_stratified_table -i "$WORK_DIR/all_ko_cpm.tsv" -o "$WORK_DIR/ko_split_cpm/"
  humann_split_stratified_table -i "$WORK_DIR/all_ko.tsv"     -o "$WORK_DIR/ko_split_rpk/"
fi
if [[ $HAS_EC -eq 1 ]]; then
  humann_split_stratified_table -i "$WORK_DIR/all_ec_cpm.tsv" -o "$WORK_DIR/ec_split_cpm/"
  humann_split_stratified_table -i "$WORK_DIR/all_ec.tsv"     -o "$WORK_DIR/ec_split_rpk/"
fi
 
# ---- 5. Conditionally rename (only when IDs are bare) ----------------------
# Modern HUMAnN3 pathway IDs already carry ": <name>" suffixes, but KO and EC
# IDs come bare. Re-renaming a named file silently appends ": NO_NAME" and
# destroys the existing names, so we check each file first.
echo "[5/6] Renaming files that need it ..."
 
needs_rename () {
  # Returns 0 (true) if <50% of non-special rows carry ": " in the ID column.
  local f=$1
  local named total
  total=$(tail -n +2 "$f" | grep -cv "^UNMAPPED\|^UNINTEGRATED" || true)
  [[ $total -eq 0 ]] && return 1
  named=$(tail -n +2 "$f" | grep -v "^UNMAPPED\|^UNINTEGRATED" \
          | cut -f1 | grep -c ": " || true)
  awk -v n=$named -v t=$total 'BEGIN{exit !(n/t < 0.5)}'
}
 
rename_if_needed () {
  local in=$1 out=$2 mapping=$3 label=$4
  if needs_rename "$in"; then
    echo "  $label: renaming ($(basename "$in"))"
    humann_rename_table -i "$in" -o "$out" --names "$mapping"
  else
    echo "  $label: already named, copying through"
    cp "$in" "$out"
  fi
}
 
rename_if_needed \
  "$WORK_DIR/path_split_cpm/all_pathabundance_cpm_unstratified.tsv" \
  "$WORK_DIR/path_split_cpm/all_pathabundance_cpm_named_unstratified.tsv" \
  metacyc-pwy "pathways"
 
rename_if_needed \
  "$WORK_DIR/rxn_split_cpm/all_rxn_cpm_unstratified.tsv" \
  "$WORK_DIR/rxn_split_cpm/all_rxn_cpm_named_unstratified.tsv" \
  metacyc-rxn "reactions"
 
if [[ $HAS_KO -eq 1 ]]; then
  rename_if_needed \
    "$WORK_DIR/ko_split_cpm/all_ko_cpm_unstratified.tsv" \
    "$WORK_DIR/ko_split_cpm/all_ko_cpm_named_unstratified.tsv" \
    kegg-orthology "KOs"
fi
if [[ $HAS_EC -eq 1 ]]; then
  rename_if_needed \
    "$WORK_DIR/ec_split_cpm/all_ec_cpm_unstratified.tsv" \
    "$WORK_DIR/ec_split_cpm/all_ec_cpm_named_unstratified.tsv" \
    ec "ECs"
fi
 
# ---- 6. Sanity-check name presence -----------------------------------------
echo "[6/6] Verifying that all output files carry names ..."
check_named () {
  local f=$1
  [[ -f "$f" ]] || return 0
  local with_name total
  total=$(tail -n +2 "$f" | grep -cv "^UNMAPPED\|^UNINTEGRATED" || true)
  with_name=$(tail -n +2 "$f" | grep -v "^UNMAPPED\|^UNINTEGRATED" \
              | cut -f1 | grep -c ": " || true)
  if [[ $total -gt 0 ]]; then
    printf "  %-65s named: %d / %d (%.0f%%)\n" \
      "$(basename "$f")" "$with_name" "$total" \
      "$(awk -v a=$with_name -v b=$total 'BEGIN{print (a/b)*100}')"
  fi
}
check_named "$WORK_DIR/path_split_cpm/all_pathabundance_cpm_named_unstratified.tsv"
[[ $HAS_KO -eq 1 ]] && check_named "$WORK_DIR/ko_split_cpm/all_ko_cpm_named_unstratified.tsv"
[[ $HAS_EC -eq 1 ]] && check_named "$WORK_DIR/ec_split_cpm/all_ec_cpm_named_unstratified.tsv"
 
cat <<EOF
 
Post-processing complete. Use these files in the Rmd:
  proc_dir <- "$WORK_DIR"
  path_cpm_unstrat_file <- file.path(proc_dir, "path_split_cpm", "all_pathabundance_cpm_named_unstratified.tsv")
  path_cpm_strat_file   <- file.path(proc_dir, "path_split_cpm", "all_pathabundance_cpm_stratified.tsv")
  ko_cpm_unstrat_file   <- file.path(proc_dir, "ko_split_cpm",   "all_ko_cpm_named_unstratified.tsv")
  path_rpk_unstrat_file <- file.path(proc_dir, "path_split_rpk", "all_pathabundance_unstratified.tsv")
  ko_rpk_unstrat_file   <- file.path(proc_dir, "ko_split_rpk",   "all_ko_unstratified.tsv")
EOF
