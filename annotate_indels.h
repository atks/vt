#ifndef ANNOTATE_INDELS_H
#define ANNOTATE_INDELS_H

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include <string>
#include <iostream>
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "hts_utils.h"
#include "utils.h"
#include "interval_tree.h"
#include "variant_manip.h"
#include "program.h"
#include "variant_manip.h"
#include "log_tool.h"

void annotate_indels(int argc, char ** argv);

#endif