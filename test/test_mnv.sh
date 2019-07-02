#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
VT=${DIR}/../vt

. ${DIR}/ssshtest

run simple_mnv_mnv_mode_off ${VT} decompose_blocksub ${DIR}/decompose_blocksub/mnv/simple_mnv.vcf
assert_exit_code 0
assert_in_stdout "1	159030	.	T	C"
assert_in_stdout "1	159031	.	A	G"

run simple_mnv ${VT} decompose_blocksub -m ${DIR}/decompose_blocksub/mnv/simple_mnv.vcf
assert_exit_code 0
assert_in_stdout "1	159030	.	TA	CG"

run 1bp_dist_mnv ${VT} decompose_blocksub -m ${DIR}/decompose_blocksub/mnv/1bp_dist_mnv.vcf
assert_exit_code 0
assert_in_stdout "1	159030	.	TAA	GAT"

run 2bp_dist_mnv ${VT} decompose_blocksub -m ${DIR}/decompose_blocksub/mnv/2bp_dist_mnv.vcf
assert_exit_code 0
assert_in_stdout "1	159030	.	T	G"
assert_in_stdout "1	159033	.	C	T"

run complex_mnv ${VT} decompose_blocksub -m ${DIR}/decompose_blocksub/mnv/complex_mnv.vcf
assert_exit_code 0
assert_in_stdout "1	159030	.	T	G"
assert_in_stdout "1	159033	.	CC	TG"