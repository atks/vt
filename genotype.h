#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <regex.h>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "program.h"
#include "hts_utils.h"
#include "bam_ordered_reader.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "variant_manip.h"
#include "utils.h"
#include "lhmm.h"
#include "log_tool.h"

void genotype(int argc, char ** argv);


#endif