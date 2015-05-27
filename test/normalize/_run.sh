#!/bin/bash

# get path to this file
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
# include helper functions
. ${DIR}/../common/_helpers.sh
# put path to vt binary into handy variable
VT=${DIR}/../../vt
# put path to vt binary into handy variable
REF=${DIR}/../ref/20.fa.gz
STRIP_VAR=${DIR}/../common/strip_var
# immediately quit on errors
set -e

# create temporary directory and ensure cleanup on termination
export TMPDIR=${DIR}/tmp
mkdir -p ${TMPDIR}
trap "set -x; rm -rf ${TMPDIR}" EXIT KILL TERM INT HUP

echo "----------------------" >&2
echo "Tests for vt normalize" >&2
echo "----------------------" >&2

# ---------------------------------------------------------------------------
# Test 01: normalization of indels from Mills et al. chromosome 20
# ---------------------------------------------------------------------------

# call program

set -x

${VT} \
    normalize \
    ${DIR}/01_IN.vcf \
    -r ${REF} \
    -o ${TMPDIR}/01_OUT.vcf \
    2>&1 | ${STRIP_VAR} > ${TMPDIR}/01_OUT.stderr

# compare results
diff ${DIR}/01_OUT.vcf ${TMPDIR}/01_OUT.vcf
diff ${DIR}/01_OUT.stderr ${TMPDIR}/01_OUT.stderr
   
set +x
