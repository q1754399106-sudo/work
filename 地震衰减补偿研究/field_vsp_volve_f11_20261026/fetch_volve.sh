#!/bin/bash
# Download the files in fetch_list.txt from the Volve Databricks volume into raw/ (skips files already complete).
E="$(cd "$(dirname "$0")" && pwd)"
DB=$(ls -d /c/Users/qq175/AppData/Local/Microsoft/WinGet/Packages/Databricks.DatabricksCLI_*/ | head -1)databricks.exe
V=dbfs:/Volumes/equinor_asa_volve_data_village/public/volve
W="$V/Well_logs_pr_WELL/15_9-F-11 T2"
while IFS='|' read -r src dst; do
  [ -z "$src" ] && continue
  if [[ "$src" == @* ]]; then url="$V/${src#@}"; else url="$W/$src"; fi
  out="$E/raw/$dst"; mkdir -p "$(dirname "$out")"
  for try in 1 2 3 4; do
    "$DB" fs cp "$url" "$out" --overwrite >/dev/null 2>&1 && break
    echo "retry $try $dst"; sleep $((2**try))
  done
  [ -f "$out" ] && echo "ok $(stat -c %s "$out") $dst" || echo "FAIL $dst"
done < "$E/fetch_list.txt"
