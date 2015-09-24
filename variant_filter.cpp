/* The MIT License

   Copyright (c) 2015 Adrian Tan <atks@umich.edu>

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

#include "variant_filter.h"

VariantFilter::VariantFilter()
{
    snp_p = 0.5-reference_bias;
    del_p = 0.5-reference_bias;
    ins_p = 0.5-reference_bias;

    snp_e = 0.003;
    del_e = 0.00023;
    ins_e = 0.00017;

    snp_adaptive_cutoff = false;
    deletion_adaptive_cutoff = false;
    insertion_adaptive_cutoff = false;
    
    lr_cutoff = -1;
}

/**
 * Filters a SNP.
 */
bool VariantFilter::filter_snp(uint32_t evidence_no, uint32_t read_no)
{
    if (snp_adaptive_cutoff)
    {
        return evidence_no>=snp_e_cutoff && pbinom(evidence_no, read_no, snp_p, 1, 0)>=snp_desired_type_II_error;
    }
    else
    {
        return evidence_no>=snp_e_cutoff && (evidence_no/(float)read_no)>=snp_f_cutoff;
    }
}

/**
 * Filters a deletion.
 *
 * Returns true if variant is to be discarded.
 */
bool VariantFilter::filter_del(uint32_t evidence_no, uint32_t read_no)
{
    if (deletion_adaptive_cutoff)
    {
        return evidence_no>=deletion_e_cutoff && pbinom(evidence_no, read_no, del_p, 1, 0)>=deletion_desired_type_II_error;
    }
    else
    {
        return evidence_no>=deletion_e_cutoff && (evidence_no/(float)read_no)>=deletion_f_cutoff;
    }
}

/**
 * Filters an insertion.
 */
bool VariantFilter::filter_ins(uint32_t evidence_no, uint32_t read_no)
{
    if (insertion_adaptive_cutoff)
    {
        return evidence_no>=insertion_e_cutoff && pbinom(evidence_no, read_no, ins_p, 1, 0)>=insertion_desired_type_II_error;
    }
    else
    {
        return evidence_no>=insertion_e_cutoff && (evidence_no/(float)read_no)>=insertion_f_cutoff;
    }
}

/**
 * Sync variables.
 */
void VariantFilter::sync()
{
    snp_adaptive_cutoff = snp_desired_type_I_error!=0 || snp_desired_type_II_error!=0;
    deletion_adaptive_cutoff = deletion_desired_type_I_error!=0 || deletion_desired_type_II_error!=0;
    insertion_adaptive_cutoff = insertion_desired_type_I_error!=0 || insertion_desired_type_II_error!=0;
    
    if (lr_cutoff!=-1)
    {
        snp_adaptive_cutoff = true;
        deletion_adaptive_cutoff = true;
        insertion_adaptive_cutoff = true;
    }
}

/**
 * Setters for reference bias.
 */
void VariantFilter::set_reference_bias(float reference_bias)
{
    this->reference_bias = reference_bias;
    snp_p = 0.5-reference_bias;
    del_p = 0.5-reference_bias;
    ins_p = 0.5-reference_bias;
}

float VariantFilter::get_reference_bias()
{
    return this->reference_bias;
}

/**
 * Setters for general filters.
 */
void VariantFilter::set_lr_cutoff(float lr_cutoff)
{
     this->lr_cutoff = lr_cutoff;
}

float VariantFilter::get_lr_cutoff()
{
    return lr_cutoff;
}

/**
 * Setters for SNP filters.
 */
void VariantFilter::set_snp_baseq_cutoff(uint32_t snp_baseq_cutoff)
{
    this->snp_baseq_cutoff = snp_baseq_cutoff;
}

void VariantFilter::set_snp_e_cutoff(uint32_t snp_e_cutoff)
{
    this->snp_e_cutoff = snp_e_cutoff;
}

void VariantFilter::set_snp_f_cutoff(float snp_f_cutoff)
{
    this->snp_f_cutoff = snp_f_cutoff;
}

void VariantFilter::set_snp_desired_type_I_error(float snp_desired_type_I_error)
{
    this->snp_desired_type_I_error = snp_desired_type_I_error;
}

void VariantFilter::set_snp_desired_type_II_error(float snp_desired_type_II_error)
{
    this->snp_desired_type_II_error = snp_desired_type_II_error;
}

uint32_t VariantFilter::get_snp_baseq_cutoff()
{
    return snp_baseq_cutoff;
}

uint32_t VariantFilter::get_snp_e_cutoff()
{
    return snp_e_cutoff;
}

float VariantFilter::get_snp_f_cutoff()
{
    return snp_f_cutoff;
}

float VariantFilter::get_snp_desired_type_I_error()
{
    return snp_desired_type_I_error;
}

float VariantFilter::get_snp_desired_type_II_error()
{
    return snp_desired_type_II_error;
}

/**
 * Setters for deletion filters.
 */
void VariantFilter::set_deletion_e_cutoff(uint32_t deletion_e_cutoff)
{
    this->deletion_e_cutoff = deletion_e_cutoff;
}

void VariantFilter::set_deletion_f_cutoff(float deletion_f_cutoff)
{
    this->deletion_f_cutoff = deletion_f_cutoff;
}

void VariantFilter::set_deletion_desired_type_I_error(float deletion_desired_type_I_error)
{
    this->deletion_desired_type_I_error = deletion_desired_type_I_error;
}

void VariantFilter::set_deletion_desired_type_II_error(float deletion_desired_type_II_error)
{
    this->deletion_desired_type_II_error = deletion_desired_type_II_error;
}

uint32_t VariantFilter::get_deletion_e_cutoff()
{
    return deletion_e_cutoff;
}

float VariantFilter::get_deletion_f_cutoff()
{
    return deletion_f_cutoff;
}

float VariantFilter::get_deletion_desired_type_I_error()
{
    return deletion_desired_type_I_error;
}

float VariantFilter::get_deletion_desired_type_II_error()
{
    return deletion_desired_type_II_error;
}

/**
 * Setters and getters for insertion filters.
 */
void VariantFilter::set_insertion_e_cutoff(uint32_t insertion_e_cutoff)
{
    this->insertion_e_cutoff = insertion_e_cutoff;
}

void VariantFilter::set_insertion_f_cutoff(float insertion_f_cutoff)
{
    this->insertion_f_cutoff = insertion_f_cutoff;
}

void VariantFilter::set_insertion_desired_type_I_error(float insertion_desired_type_I_error)
{
    this->insertion_desired_type_I_error = insertion_desired_type_I_error;
}

void VariantFilter::set_insertion_desired_type_II_error(float insertion_desired_type_II_error)
{
    this->insertion_desired_type_II_error = insertion_desired_type_II_error;
}

uint32_t VariantFilter::get_insertion_e_cutoff()
{
    return insertion_e_cutoff;
}

float VariantFilter::get_insertion_f_cutoff()
{
    return insertion_f_cutoff;
}

float VariantFilter::get_insertion_desired_type_I_error()
{
    return insertion_desired_type_I_error;
}

float VariantFilter::get_insertion_desired_type_II_error()
{
    return insertion_desired_type_II_error;
}

/**
 * Setters and getters for soft clips filters.
 */
void VariantFilter::set_sclip_mq_cutoff(float sclip_mq_cutoff)
{
    this->sclip_mq_cutoff = sclip_mq_cutoff;
}

void VariantFilter::set_sclip_u_cutoff(float sclip_u_cutoff)
{
    this->sclip_u_cutoff = sclip_u_cutoff;
}

float VariantFilter::get_sclip_mq_cutoff()
{
    return sclip_mq_cutoff;
}

uint32_t VariantFilter::get_sclip_u_cutoff()
{
    return sclip_u_cutoff;
}