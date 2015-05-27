#!/bin/bash

set -x
set -e

# get path to this file
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
# include helper functions
. ${DIR}/_helpers.sh
# put path to vt binary into handy variable
VT=${DIR}/../../vt

# ---------------------------------------------------------------------------
# Test 01: block decomposition of even-length blocks
# ---------------------------------------------------------------------------

${VT} \
    decompose_blocksub \
    ${DIR}/01_IN_even_length.vcf \
    >${DIR}/01_OUT_even_length.vcf.stdout \
    2> >(strip_stderr > ${DIR}/01_OUT_even_length.vcf.stderr)

# ---------------------------------------------------------------------------
# Test 02: block decomposition with alignment
# ---------------------------------------------------------------------------

${VT} \
    decompose_blocksub \
    -a \
    ${DIR}/02_IN_uneven_length.vcf \
    >${DIR}/02_OUT_uneven_length.vcf.stdout \
    2> >(strip_stderr > ${DIR}/02_OUT_uneven_length.vcf.stderr)
