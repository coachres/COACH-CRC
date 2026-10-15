#!/bin/bash
# =============================================================================
# Analysis pipeline: ascatNgs (CNV)  +  BRASS (SV)
# =============================================================================

set -euo pipefail

# =========================== Path config ===========================
CANCERIT_TOOLS="/opt/bioinfo/CancerIT-BGI/v20210714"            # CancerIT tool root (contains ascatNgs / BRASS)
CANCERIT_DB="$CANCERIT_TOOLS/database/GRCh38_hla_decoy_ebv" # ascat / brass / vagrent database directory
R_BIN="/opt/R/4.1.2/bin"                                # R runtime
REF_FA="/opt/bioinfo/database/Homo_sapiens_assembly38.fasta"   # reference genome fasta
ASCATNGS_DIR="/path/to/ascatNgs_result"                     # ascatNgs output root
BRASS_DIR="/path/to/BRASS_output"                           # BRASS output directory

# =========================== Environment variables  ===========================
export CGP_PERLLIBS="$CANCERIT_TOOLS/lib/perl5"
export PERL5LIB="$CANCERIT_TOOLS/lib/perl5"
export PATH="$R_BIN:$CANCERIT_TOOLS/bin:$CANCERIT_TOOLS/biobambam2/bin:$PATH"
export R_LIBS_USER="$CANCERIT_TOOLS/R-libs"
export LD_LIBRARY_PATH="$CANCERIT_TOOLS/lib:$LD_LIBRARY_PATH"

# =========================== Input arguments ===========================
# Usage:
#   bash runBRASS.sh <normal_bam> <tumor_bam> <sample_name> [normal_name]
# Arguments:
#     $1  normal_bam    normal control BAM
#     $2  tumor_bam     tumor BAM
#     $3  sample_name   sample name (used for ascatNgs output dir and samplestatistics filename)
#     $4  normal_name   normal sample name (optional, defaults to sample_name; used for pair dir name)
usage() {
    echo "Usage: bash $0 <normal_bam> <tumor_bam> <sample_name> [normal_name]" >&2
    exit 1
}
[ "$#" -lt 3 ] && usage

NBAM=$1
TBAM=$2
TNAME=$3
NNAME=${4:-$TNAME}

echo "hostname: $(hostname)"
echo "start:    $(date)"

# =========================== Step 1: ascatNgs (CNV) ===========================
PAIRNAME="${NNAME}-VS-${TNAME}"
ASCAT_OUT="$ASCATNGS_DIR/$PAIRNAME"
mkdir -p "$ASCAT_OUT"

ascat.pl -g L -q 20 -rs 'human' -ra GRCh38 -pr WGS -c 4 -force -nobigwig \
    -o "$ASCAT_OUT" \
    -t "$TBAM" \
    -n "$NBAM" \
    -r "$REF_FA" \
    -sg "$CANCERIT_DB/ascat/SnpGcCorrections.tsv" \
    -l  "$CANCERIT_DB/ascat/hg38_chrY.loci"

if [ -d "$ASCAT_OUT/logs" ]; then
    ( cd "$ASCAT_OUT" && tar -czf logs.tar.gz logs && rm -rf logs )
fi

echo "ascatNgs done: $(date)"

# =========================== Step 2: BRASS (SV) ===========================
# BRASS depends on the samplestatistics file produced by ascatNgs
if [ ! -e "$NBAM.bas" ]; then
    bam_stats --input "$NBAM" --output "$NBAM.bas" --ref-file "$REF_FA" --num_threads 24
fi
echo "nbam stats done: $(date)"

if [ ! -e "$TBAM.bas" ]; then
    bam_stats --input "$TBAM" --output "$TBAM.bas" --ref-file "$REF_FA" --num_threads 24
fi
echo "tbam stats done: $(date)"

brass.pl -j 4 -c 24 -s 'human' -as GRCh38 -pr WGS \
    -t "$TBAM" -n "$NBAM" -g "$REF_FA" -o "$BRASS_DIR" \
    -d "$CANCERIT_DB/brass/HiDepth.bed.gz" \
    -f "$CANCERIT_DB/brass/brass_np.groups.gz" \
    -g_cache "$CANCERIT_DB/vagrent/vagrent.human.GRCh38.homo_sapiens_core_104_38.cache.gz" \
    -vi "$CANCERIT_DB/brass/viral.genomic.fa.2bit" \
    -mi "$CANCERIT_DB/brass/all_ncbi_bacteria" \
    -b "$CANCERIT_DB/brass/500bp_windows.gc.bed.gz" \
    -ct "$CANCERIT_DB/brass/CentTelo.tsv" \
    -cb "$CANCERIT_DB/brass/cytoband.txt" \
    -ss "$ASCAT_OUT/$TNAME.samplestatistics.txt"

echo "BRASS done: $(date)"
