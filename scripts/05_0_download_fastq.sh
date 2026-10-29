#!/usr/bin/env bash
# Downloads FASTQ files from the ENA FTP server, several at a time, and checks each file against its MD5 sum.
# Runs on the machine that does the quantification (05_1 header), with bash, wget or curl, md5sum and xargs.
#
# Usage:  bash scripts/05_0_download_fastq.sh LIST.tsv OUT_DIR [N_PARALLEL]
#   LIST.tsv    tab-separated, with a header holding the columns ftp (host and path of one file, as in the
#               fastq_ftp field of the ENA file report) and md5, for example
#               ncbi_samples/PRJNA1160017_skin_stomach_fastq.tsv
#   OUT_DIR     where the files go (created if missing)
#   N_PARALLEL  files downloaded at the same time (default 4)
# Run it again after an interruption: complete files with a matching MD5 are skipped and partial files are resumed.
# A file whose MD5 does not match is downloaded again from the start, up to three times.
# Output: OUT_DIR/<file name>, OUT_DIR/download_log.tsv (file, status), failed files listed at the end.

set -euo pipefail
if [ $# -lt 2 ]; then
  sed -n '5,10p' "$0"
  exit 1
fi
list=$1
out=$2
jobs=${3:-4}
mkdir -p "$out"
log="$out/download_log.tsv"
printf 'file\tstatus\n' > "$log"

if command -v wget > /dev/null; then
  export DL=wget
elif command -v curl > /dev/null; then
  export DL=curl
else
  echo "Needs wget or curl" >&2
  exit 1
fi
export OUT=$out LOG=$log

# One file: skip if complete, else download (resuming a partial file) and check the MD5
fetch() {
  local url=$1 md5=$2 file
  file="$OUT/$(basename "$url")"
  if [ -s "$file" ] && [ "$(md5sum "$file" | cut -d' ' -f1)" = "$md5" ]; then
    printf '%s\tpresent\n' "$(basename "$file")" >> "$LOG"
    return 0
  fi
  for attempt in 1 2 3; do
    if [ "$DL" = wget ]; then
      wget -q -c --tries=5 --waitretry=30 -P "$OUT" "ftp://$url" || true
    else
      curl -sS -C - --retry 5 --retry-delay 30 -o "$file" "ftp://$url" || true
    fi
    if [ -s "$file" ] && [ "$(md5sum "$file" | cut -d' ' -f1)" = "$md5" ]; then
      printf '%s\tok\n' "$(basename "$file")" >> "$LOG"
      echo "done   $(basename "$file")"
      return 0
    fi
    echo "retry  $(basename "$file") (attempt $attempt, MD5 mismatch or download error)" >&2
    rm -f "$file"
  done
  printf '%s\tfailed\n' "$(basename "$file")" >> "$LOG"
  echo "FAILED $(basename "$file")" >&2
}
export -f fetch

# url and md5 columns, found by their header names
awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
            { print $col["ftp"], $col["md5"] }' "$list" |
  xargs -P "$jobs" -n 2 bash -c 'fetch "$@"' _

n_failed=$(awk -F'\t' '$2 == "failed"' "$log" | wc -l)
echo "$(($(wc -l < "$log") - 1)) files, $n_failed failed (see $log)"
[ "$n_failed" -eq 0 ]
