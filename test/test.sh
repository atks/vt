#!/bin/bash

# Removes the time from vt stderr output.
strip_stderr()
{
    sed 's/Time elapsed.*/Time elapsed <stripped>/g' | sed 's/file.*/file <stripped>/g'
}

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
VT=${DIR}/../vt
REF=${DIR}/ref/20.fa.gz
TMPDIRS="";

NO_TESTS=0
PASSED_TESTS=0

echo "++++++++++++++++++++++" >&2
echo "Tests for vt normalize" >&2
echo "++++++++++++++++++++++" >&2

# create temporary directory and ensure cleanup on termination
CMDDIR=${DIR}/normalize
TMPDIR=${CMDDIR}/tmp
mkdir -p ${TMPDIR}
TMPDIRS+=" $TMPDIR";

#-----------------------
echo "testing normalize"
#-----------------------
if [ "$1" == "debug" ]; then
    set -x
fi

${OUT}${VT} \
            normalize \
            ${CMDDIR}/01_IN.vcf \
            -r ${REF} \
            -o ${TMPDIR}/01_OUT.vcf \
            2>&1 | strip_stderr > ${TMPDIR}/01_OUT.stderr

OUT=`diff ${CMDDIR}/01_OUT.vcf ${TMPDIR}/01_OUT.vcf`
ERR=`diff ${CMDDIR}/01_OUT.stderr ${TMPDIR}/01_OUT.stderr`

set +x

((NO_TESTS++))

echo -n "             output VCF file :"
if [ "$OUT" == "" ]; then
    echo " ok"
    ((PASSED_TESTS++))
else
    echo " NOT OK!!!"
fi

echo -n "             output logs     :"
if [ "$ERR" == "" ]; then
    echo " ok"
else
    echo " NOT OK!!!"
fi

echo "+++++++++++++++++++++++++++++++" >&2
echo "Tests for vt decompose_blocksub" >&2
echo "+++++++++++++++++++++++++++++++" >&2

# create temporary directory and ensure cleanup on termination
CMDDIR=${DIR}/decompose_blocksub
TMPDIR=${CMDDIR}/tmp
mkdir -p ${TMPDIR}
TMPDIRS+=" $TMPDIR";

#------------------------------------------------------ 
echo "testing decompose_blocksub of even-length blocks"
#------------------------------------------------------ 

if [ "$1" == "debug" ]; then
    set -x
fi

${VT} \
    decompose_blocksub \
    ${CMDDIR}/01_IN_even_length.vcf \
    -o ${TMPDIR}/01_OUT_even_length.vcf \
    2>&1 | strip_stderr > ${TMPDIR}/01_OUT_even_length.stderr

OUT=`diff ${CMDDIR}/01_OUT_even_length.vcf ${TMPDIR}/01_OUT_even_length.vcf`
ERR=`diff ${CMDDIR}/01_OUT_even_length.stderr ${TMPDIR}/01_OUT_even_length.stderr`

set +x

((NO_TESTS++))

echo -n "             output VCF file :"
if [ "$OUT" == "" ]; then
    echo " ok"
    ((PASSED_TESTS++))
else
    echo " NOT OK!!!"
fi

echo -n "             output logs     :"
if [ "$ERR" == "" ]; then
    echo " ok"
else
    echo " NOT OK!!!"
fi

#----------------------------------------------- 
echo "testing decompose_blocksub with alignment"
#----------------------------------------------- 

if [ "$1" == "debug" ]; then
    set -x
fi

${VT} \
    decompose_blocksub -a \
    ${CMDDIR}/02_IN_uneven_length.vcf \
    -o ${TMPDIR}/02_OUT_uneven_length.vcf \
    2>&1 | strip_stderr > ${TMPDIR}/02_OUT_uneven_length.stderr \

OUT=`diff ${CMDDIR}/02_OUT_uneven_length.vcf ${TMPDIR}/02_OUT_uneven_length.vcf`
ERR=`diff ${CMDDIR}/02_OUT_uneven_length.stderr ${TMPDIR}/02_OUT_uneven_length.stderr`

set +x

((NO_TESTS++))

echo -n "             output VCF file :"
if [ "$OUT" == "" ]; then
    echo " ok"
    ((PASSED_TESTS++))
else
    echo " NOT OK!!!"
fi

echo -n "             output logs     :"
if [ "$ERR" == "" ]; then
    echo " ok"
else
    echo " NOT OK!!!"
fi

#------------------------------------------------------------- 
echo "testing decompose_blocksub of phased even-length blocks"
#------------------------------------------------------------- 

if [ "$1" == "debug" ]; then
    set -x
fi

${VT} \
    decompose_blocksub \
    -p \
    ${CMDDIR}/03_IN_phased_even_length.vcf \
    -o ${TMPDIR}/03_OUT_phased_even_length.vcf \
    2>&1 | strip_stderr > ${TMPDIR}/03_OUT_phased_even_length.stderr

OUT=`diff ${CMDDIR}/03_OUT_phased_even_length.vcf ${TMPDIR}/03_OUT_phased_even_length.vcf`
ERR=`diff ${CMDDIR}/03_OUT_phased_even_length.stderr ${TMPDIR}/03_OUT_phased_even_length.stderr`

set +x

((NO_TESTS++))

echo -n "             output VCF file :"
if [ "$OUT" == "" ]; then
    echo " ok"
    ((PASSED_TESTS++))
else
    echo " NOT OK!!!"
fi

echo -n "             output logs     :"
if [ "$ERR" == "" ]; then
    echo " ok"
else
    echo " NOT OK!!!"
fi

if [ "$1" != "debug" ]; then
    trap "rm -rf ${TMPDIRS}" EXIT KILL TERM INT HUP
fi

echo "++++++++++++++++++++++" >&2
echo "Tests for vt decompose" >&2
echo "++++++++++++++++++++++" >&2

# create temporary directory and ensure cleanup on termination
CMDDIR=${DIR}/decompose
TMPDIR=${CMDDIR}/tmp
mkdir -p ${TMPDIR}
TMPDIRS+=" $TMPDIR";

#-----------------------
echo "testing decompose for a triallelic variant"
#----------------------- 

if [ "$1" == "debug" ]; then
    set -x
fi

${VT} \
    decompose \
    ${CMDDIR}/01_IN_multi.vcf \
    -o ${TMPDIR}/01_OUT_multi.vcf \
    2>&1 | strip_stderr > ${TMPDIR}/01_OUT_multi.stderr

OUT=`diff ${CMDDIR}/01_OUT_multi.vcf ${TMPDIR}/01_OUT_multi.vcf`
ERR=`diff ${CMDDIR}/01_OUT_multi.stderr ${TMPDIR}/01_OUT_multi.stderr`

set +x

((NO_TESTS++))

echo -n "             output VCF file :"
if [ "$OUT" == "" ]; then
    echo " ok"
    ((PASSED_TESTS++))
else
    echo " NOT OK!!!"
fi

echo -n "             output logs     :"
if [ "$ERR" == "" ]; then
    echo " ok"
else
    echo " NOT OK!!!"
fi

if [ "$1" != "debug" ]; then
    trap "rm -rf ${TMPDIRS}" EXIT KILL TERM INT HUP
fi

echo
echo -n Passed tests :
echo -n " "
echo -n ${PASSED_TESTS}  
echo -n " / "
echo ${NO_TESTS}  
