#!/bin/bash
# Download the composite logs in fetch_list.txt and the Volve well picks into raw/ (skips complete files).
E="$(cd "$(dirname "$0")" && pwd)"
DB=$(ls -d /c/Users/qq175/AppData/Local/Microsoft/WinGet/Packages/Databricks.DatabricksCLI_*/ | head -1)databricks.exe
V=dbfs:/Volumes/equinor_asa_volve_data_village/public/volve
get() {
  for try in 1 2 3 4; do
    "$DB" fs cp "$1" "$2" --overwrite >/dev/null 2>&1 && break
    echo "retry $try $2"; sleep $((2**try))
  done
  [ -f "$2" ] && echo "ok $(stat -c %s "$2") $2" || echo "FAIL $2"
}
mkdir -p "$E/raw/picks"
get "$V/Geophysical_Interpretations/Wells/Well_picks_Volve_v1.dat" "$E/raw/picks/Well_picks_Volve_v1.dat"
get "$V/Geophysical_Interpretations/Wells/README.txt" "$E/raw/picks/README.txt"
while IFS='|' read -r well file; do
  [ -z "$well" ] && continue
  d="$E/raw/composite/${well// /_}"; mkdir -p "$d"
  [ -s "$d/$file" ] || get "$V/Well_logs_pr_WELL/$well/04.COMPOSITE/$file" "$d/$file"
done < "$E/fetch_list.txt"
