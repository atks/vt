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

#ifndef VARIANT_FILTER_H
#define VARIANT_FILTER_H

#include <vector>
#include <map>
#include "Rmath/Rmath.h"
#include "hts_utils.h"

/**
 * Class for filtering variants in discover
 */
class VariantFilter
{
    float snp_p;
    float del_p;
    float ins_p;

    float snp_e;
    float del_e;
    float ins_e;

    //bia
    float reference_bias;

    //general cutoffs
    float lr_cutoff;

    //snp filters
    uint32_t snp_baseq_cutoff;
    uint32_t snp_e_cutoff;
    float snp_f_cutoff;
    float snp_desired_type_I_error;
    float snp_desired_type_II_error;
    bool snp_adaptive_cutoff;

    //deletion filters
    uint32_t deletion_e_cutoff;
    float deletion_f_cutoff;
    float deletion_desired_type_I_error;
    float deletion_desired_type_II_error;
    bool deletion_adaptive_cutoff;

    //insertion filters
    uint32_t insertion_e_cutoff;
    float insertion_f_cutoff;
    float insertion_desired_type_I_error;
    float insertion_desired_type_II_error;
    bool insertion_adaptive_cutoff;

    //soft clip filters
    float sclip_mq_cutoff;
    uint32_t sclip_u_cutoff;

    public:

    /**
     * Constructor.
     */
    VariantFilter();

    /**
     * Filters a SNP.
     */
    bool filter_snp(uint32_t evidence_no, uint32_t read_no);

    /**
     * Filters a deletion.
     */
    bool filter_del(uint32_t evidence_no, uint32_t read_no);

    /**
     * Filters an insertion.
     */
    bool filter_ins(uint32_t evidence_no, uint32_t read_no);

    /**
     * Sync variables.
     */
    void sync();

    /**
     * Setters for reference bias.
     */
    void set_reference_bias(float reference_bias);

    float get_reference_bias();

    /**
     * Setters for general filters.
     */
    void set_lr_cutoff(float lr_cutoff);

    float get_lr_cutoff();

    /**
     * Setters for SNP filters.
     */
    void set_snp_baseq_cutoff(uint32_t snp_baseq_cutoff);

    void set_snp_e_cutoff(uint32_t snp_e_cutoff);

    void set_snp_f_cutoff(float snp_f_cutoff);

    void set_snp_desired_type_I_error(float snp_desired_type_I_error);

    void set_snp_desired_type_II_error(float snp_desired_type_II_error);

    uint32_t get_snp_baseq_cutoff();

    uint32_t get_snp_e_cutoff();

    float get_snp_f_cutoff();

    float get_snp_desired_type_I_error();

    float get_snp_desired_type_II_error();

    /**
     * Setters for deletion filters.
     */
    void set_deletion_e_cutoff(uint32_t deletion_e_cutoff);

    void set_deletion_f_cutoff(float deletion_f_cutoff);

    void set_deletion_desired_type_I_error(float deletion_desired_type_I_error);

    void set_deletion_desired_type_II_error(float deletion_desired_type_II_error);

    uint32_t get_deletion_e_cutoff();

    float get_deletion_f_cutoff();

    float get_deletion_desired_type_I_error();

    float get_deletion_desired_type_II_error();

    /**
     * Setters and getters for insertion filters.
     */
    void set_insertion_e_cutoff(uint32_t insertion_e_cutoff);

    void set_insertion_f_cutoff(float insertion_f_cutoff);

    void set_insertion_desired_type_I_error(float insertion_desired_type_I_error);

    void set_insertion_desired_type_II_error(float insertion_desired_type_II_error);

    uint32_t get_insertion_e_cutoff();

    float get_insertion_f_cutoff();

    float get_insertion_desired_type_I_error();

    float get_insertion_desired_type_II_error();

    /**
     * Setters and getters for soft clips filters.
     */
    void set_sclip_mq_cutoff(float sclip_mq_cutoff);

    void set_sclip_u_cutoff(float sclip_u_cutoff);

    float get_sclip_mq_cutoff();

    uint32_t get_sclip_u_cutoff();
};

#endif