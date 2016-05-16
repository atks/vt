/* The MIT License
   Copyright (c) 2016 Hyun Min Kang <atks@umich.edu>
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

#ifndef JOINT_GENOTYPING_UNIT_H
#define JOINT_GENOTYPING_UNIT_H

class JointGenotypingUnit
{
    public:
  
    bcf_hdr_t *h;
    bcf1_t *v;
    int32_t rid;
    int32_t pos1; //position of variant
    //[beg1,end1] is the required overlapping of the variant against the aligned read necessary to make a genotype call.
    //for SNPs, beg1=end1=pos1
    //
    //for Indels, this refers to the flanking positions
    //   insertion 
    //        if T/TG - beg1=pos1, end1=pos1+1
    //        if T/GT - beg1=pos1-1, end1=pos1
    //   deletion  
    //        if TG/T - beg1=pos1, end1=pos1+length(REF)
    //        if TG/G - beg1=pos1-1, end1=pos1+length(REF)-1
    int32_t beg1;
    int32_t end1;
    int32_t vtype;

    //indel specific record
    int32_t dlen;
    uint32_t len;
    std::string indel;

    //vntr specific record
    std::string motif;

    //vntr specific record
    //std::vector<float> counts;

    // sample level information
    int32_t nsamples;
    kstring_t alleles;

    uint8_t* pls;
    uint8_t* ads;
    Estimator* est;

    // sufficient statistics for computing INFO field
    float bqr_num, bqr_den;
    float mqr_num, mqr_den;
    float cyr_num, cyr_den;
    float str_num, str_den;
    float nmr_num, nmr_den;
    float ior_num, ior_den;
    float nm0_num, nm0_den;
    float nm1_num, nm1_den;
    float abe_num, abe_den;
    float abz_num, abz_den;
    float ns_nref, dp_sum, max_gq;

    int32_t tmp_dp_q20;
    int32_t tmp_dp_ra;
    int32_t tmp_bq_s1, tmp_bq_s2;
    int32_t tmp_mq_s1, tmp_mq_s2;
    float tmp_cy_s1, tmp_cy_s2;
    int32_t tmp_st_s1, tmp_st_s2;
    int32_t tmp_al_s1, tmp_bq_al, tmp_mq_al;
    float  tmp_cy_al;
    int32_t tmp_st_al, tmp_nm_al;
    int32_t tmp_nm_s1, tmp_nm_s2;
    double tmp_oth_exp_q20, tmp_oth_obs_q20;
    double tmp_pls[3];
    double tmp_ads[3];

    // temporary information to be cleared out per-sample basis
    

    /**
     * Constructor.
     * @v - VCF record.
     */
    JointGenotypingRecord(bcf_hdr_t*h, bcf1_t *v, int32_t vtype, int32_t nsamples);

    /**
     * Clears this record.
     */
    void clear();
    void clearTemp();
    bcf1_t* flush_variant(bcf_hdr_t* hdr);
    void flush_sample( int32_t sampleIndex );
    void add_allele( double contam, int32_t allele, uint8_t mapq, bool fwd, uint32_t q, uint32_t cycle, uint32_t nm );
    void process_read(AugmentedBAMRecord& as, int32_t sampleIndex, double contam);    

    /**
     * Destructor.
     */
    ~JointGenotypingRecord();

    float compute_correlation(int32_t n, int32_t xy, int32_t x1, int32_t x2, int32_t y1, int32_t y2, float buffer) {
      if ( n == 0 ) return 0;
      float xsd = x2/(float)n - (x1/(float)n)*(x1/(float)n);
      float ysd = y2/(float)n - (y1/(float)n)*(y1/(float)n);
      return ( ( xy/(float)n - x1 * y1 / (float) n / (float) n ) / sqrt( xsd * ysd + buffer ) );
    }

    float compute_correlation_f(int32_t n, float xy, float x1, float x2, float y1, float y2, float buffer) {
      if ( n == 0 ) return 0;
      float xsd = x2/(float)n - (x1/(float)n)*(x1/(float)n);
      float ysd = y2/(float)n - (y1/(float)n)*(y1/(float)n);
      return ( ( xy/(float)n - x1 * y1 / (float) n / (float) n ) / sqrt( xsd * ysd + buffer ) );
    }
    
};

#endif
