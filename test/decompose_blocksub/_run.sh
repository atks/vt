#!/bin/bash

# get path to this file
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
# include helper functions
. ${DIR}/_helpers.sh
# put path to vt binary into handy variable
VT=${DIR}/../../vt
# immediately quit on errors
set -e

# create temporary directory and ensure cleanup on termination
export TMPDIR=$(mktemp -d)
trap "set -x; rm -rf ${TMPDIR}" EXIT KILL TERM INT HUP

echo "Tests for vt decompose_blocksub" >&2
echo "-------------------------------" >&2


# ---------------------------------------------------------------------------
# Test 01: block decomposition of even-length blocks
# ---------------------------------------------------------------------------

set -x
mkdir -p ${TMPDIR}/01

# call program
${VT} \
    decompose_blocksub \
    ${DIR}/01_IN_even_length.vcf \
    >${TMPDIR}/01/01_OUT_even_length.vcf.stdout \
    2> >(strip_stderr > ${TMPDIR}/01/01_OUT_even_length.vcf.stderr)

# compare results
diff ${DIR}/01_OUT_even_length.vcf.stdout ${TMPDIR}/01/01_OUT_even_length.vcf.stdout
diff ${DIR}/01_OUT_even_length.vcf.stderr ${TMPDIR}/01/01_OUT_even_length.vcf.stderr

set +x

# ---------------------------------------------------------------------------
# Test 02: block decomposition with alignment
# ---------------------------------------------------------------------------

set -x
mkdir -p ${TMPDIR}/02

# call program
${VT} \
    decompose_blocksub -a \
    ${DIR}/02_IN_uneven_length.vcf \
    >${TMPDIR}/02/02_OUT_uneven_length.vcf.stdout \
    2> >(strip_stderr > ${TMPDIR}/02/02_OUT_uneven_length.vcf.stderr) \

# compare results
diff ${DIR}/02_OUT_uneven_length.vcf.stdout ${TMPDIR}/02/02_OUT_uneven_length.vcf.stdout
diff ${DIR}/02_OUT_uneven_length.vcf.stderr ${TMPDIR}/02/02_OUT_uneven_length.vcf.stderr

set +x
