/* The MIT License
   Copyright (c) 2014 Adrian Tan <atks@umich.edu>
   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:
   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include "genotyping_record.h"

/**
 * Constructor.
 * @v - VCF record.
 */
GenotypingRecord::GenotypingRecord(bcf_hdr_t *h, bcf1_t *v, int32_t vtype)
{
    clear();
    this->h = h;
    this->v = v;
    rid = bcf_get_rid(v);
    pos1 = bcf_get_pos1(v);
    this->vtype = vtype;

    if (vtype==VT_SNP && bcf_get_n_allele(v)==2)
    {
        rid = bcf_get_rid(v);
        beg1 = bcf_get_pos1(v);
        end1 = beg1;
    }
    else if (vtype==VT_INDEL && bcf_get_n_allele(v)==2)
    {
        rid = bcf_get_rid(v);
        char** alleles = bcf_get_allele(v);
        dlen = strlen(alleles[1])-strlen(alleles[0]);
        len = abs(dlen);

        int32_t *flanks_pos1 = NULL;
        int32_t n = 0;
                
        if (bcf_get_info_int32(h, v, "FLANKS", &flanks_pos1, &n)>0)
        {
           this->beg1 = flanks_pos1[0];
           this->end1 = flanks_pos1[1];
           free(flanks_pos1);
        }
        else
        {
            this->beg1 = bcf_get_pos1(v) - 3;
            this->end1 = bcf_get_end_pos1(v) + 3;
        }

        if (dlen>0)
        {
            indel.append(&alleles[1][1]);
        }
        else
        {
            indel.append(&alleles[0][1]);
        }
    }
    else if (vtype==VT_VNTR)
    {
        rid = bcf_get_rid(v);
        beg1 = bcf_get_pos1(v) - 1;
        end1 = bcf_get_end_pos1(v) + 1;
        
        char *motif = NULL;
        int32_t n = 0;
        
        if (bcf_get_info_string(h, v, "MOTIF", &motif, &n)>0)
        {
           this->motif.assign(motif);
           free(motif);
        }
    }
}

/**
 * Clears this record.
 */
void GenotypingRecord::clear()
{
    v =NULL;
    vtype = -1;

    no_nonref = 0;

    quals.clear();
    map_quals.clear();
    strands.clear();
    alleles.clear();
    cycles.clear();
    no_mismatches.clear();

    allele_depth_fwd.resize(2,0);
    allele_depth_rev.resize(2,0);
    depth = 0;
    depth_fwd = 0;
    depth_rev = 0;
    base_qualities_sum = 0;
}

/**
 * Destructor.
 */
GenotypingRecord::~GenotypingRecord()
{
    if (v) bcf_destroy(v);
}
