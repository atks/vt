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

#include "joint_genotyping_buffered_reader.h"

/**
 * Constructor.
 */
JointGenotypingBufferedReader::JointGenotypingBufferedReader(std::string filename, std::vector<GenomeInterval>& intervals, std::string out_vcf_file_name, int32_t nsamples)
{
    vm = new VariantManip();  

    // read input BCF files and create genotyping records for every variants
    odr = new BCFOrderedReader(filename, intervals);
    
    bcf1_t* v = bcf_init();
    while( odr->read(v) ) {
      int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
      //if ( (rand() % 10000) == 0 )
      //notice("foo");
      int32_t pos1 = bcf_get_pos1(v);
      
      if ( ( pos1 < intervals[0].start1 ) || ( pos1 > intervals[0].end1 ) ) continue;
	
      if ( ( vtype != VT_VNTR ) && ( bcf_get_n_allele(v) == 2 ) ) { 
	JointGenotypingRecord* jgr = new JointGenotypingRecord(odr->hdr, v, vtype, nsamples);
	gRecords.push_back(jgr);
      }
    }
    bcf_destroy(v);

    //odr->close();
    //delete odr;

    //////////////////////////
    //options initialization//
    //////////////////////////
    output_annotations = false;

    ////////////////////////
    //stats initialization//
    ////////////////////////
    no_snps_genotyped = 0;
    no_indels_genotyped = 0;
    no_vntrs_genotyped = 0;

    lastFirst = 0;
    //currentSampleIndex = -1;

    ////////////////////////
    //tools initialization//
    ////////////////////////
    sample_names.resize(nsamples);
    sample_contams.resize(nsamples);
}

void JointGenotypingBufferedReader::flush_sample(int32_t sampleIndex) {
  /*
  if ( sampleIndex < 0 )
    sampleIndex = currentSampleIndex;
  */
  
  for(int32_t i=0; i < (int)gRecords.size(); ++i) {
    gRecords[i]->flush_sample( sampleIndex );
  }
}

void JointGenotypingBufferedReader::set_sample(int32_t sampleIndex, const char* sampleName, double contam) {
  lastFirst = 0;
  sample_names[sampleIndex] = sampleName;
  sample_contams[sampleIndex] = contam;
  //currentSampleIndex = sampleIndex;
}

/**
 * Collects sufficient statistics from read for variants to be genotyped.
 *
 * The VCF records in the buffer must never occur before
 */
int32_t JointGenotypingBufferedReader::process_read(bam_hdr_t *h, bam1_t *s, int32_t sampleIndex)
{
    //wrap bam1_t in AugmentBAMRecord
    as.initialize(h, s);

    uint32_t tid = bam_get_tid(s);
    uint32_t beg1 = as.beg1;
    uint32_t end1 = as.end1;

    int32_t nvisited = 0;

    //collect statistics for variant records that are in the buffer and overlap with the read
    JointGenotypingRecord* jgr;
    for(int32_t i = lastFirst; i < (int)gRecords.size(); ++i) {
      jgr = gRecords[i];
      
      //same chromosome
      if (tid == jgr->rid) {
	if (end1 < jgr->beg1) // read ends earlier than the last record to visit -- no need to investigate
	  return nvisited;
	else if (beg1 > jgr->end1) { // read begins later than the last record to visit -- advance the last record
	  ++lastFirst;
	  continue;
	}
	else if ((beg1 <= jgr->pos1) && (jgr->pos1 <= end1)) { // variant position is covered by the read
	  ++nvisited;
	  //jgr->process_read(as, currentSampleIndex, sample_contams[currentSampleIndex]);
	  jgr->process_read(as, sampleIndex, sample_contams[sampleIndex]);
	  //	  if (beg1 <= jgr->beg1 && jgr->end1 <= end1) {
	  //                    std::cerr << "COMPLETE";
	  //                }
	  //                else
	  //                {
	  //drop
	  //                    std::cerr << "PARTIAL";
          //      }
          //  }
          //  else
          //  {
          //      //other types of overlap, just ignore
          //  }
	  //            std::cerr << "\n";
        }
      }
      else if ( tid < jgr->rid )
	return nvisited;
      else if ( tid > jgr->rid ) {
	++lastFirst;	
	continue;
      }
      else
	abort(); 
    }
    return nvisited;
    
    //this means end of file
    //bcf_destroy(v);
}

/**
 * Flush records.
 */
bcf1_t* JointGenotypingBufferedReader::flush_variant(int32_t variantIndex, bcf_hdr_t* hdr) {
  return gRecords[variantIndex]->flush_variant(hdr);
}

void JointGenotypingBufferedReader::write_header(BCFOrderedWriter* odw) {
  // contig and sample names must be added beforehand
  
  bcf_hdr_append(odw->hdr, "##INFO=<ID=AVGDP,Number=1,Type=Float,Description=\"Average Depth per Sample\">\n");	
  bcf_hdr_append(odw->hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Number Allele Counts\">\n");
  //bcf_hdr_append(odw->hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Reads\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequency from Best-guess Genotypes\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=GC,Number=G,Type=Integer,Description=\"Genotype Counts\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=GN,Number=1,Type=Integer,Description=\"Total Number of Genotypes\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=HWEAF,Number=A,Type=Float,Description=\"Genotype likelihood based Allele Frequency assuming HWE\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=HWDGF,Number=G,Type=Float,Description=\"Genotype likelihood based Genotype Frequency ignoring HWE\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=IBC,Number=1,Type=Float,Description=\"Inbreeding Coefficients calculated from genotype likelihoods\">\n");	
  bcf_hdr_append(odw->hdr, "##INFO=<ID=HWE_SLP,Number=1,Type=Float,Description=\"Signed log p-values testing  statistics based Hardy Weinberg ln(Likelihood Ratio)\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=ABE,Number=1,Type=Float,Description=\"Expected allele Balance towards Reference Allele on Heterozygous Sites\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=ABZ,Number=1,Type=Float,Description=\"Average Z-scores of Allele Balance towards Reference Allele on Heterozygous Sites\">\n");	
  bcf_hdr_append(odw->hdr, "##INFO=<ID=NS_NREF,Number=1,Type=Integer,Description=\"Number of samples with non-reference reads\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=BQZ,Number=1,Type=Float,Description=\"Correlation between base quality and alleles\">\n");
  //bcf_hdr_append(odw->hdr, "##INFO=<ID=MQZ,Number=1,Type=Float,Description=\"Correlation between mapping quality and alleles\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=CYZ,Number=1,Type=Float,Description=\"Correlation between cycle and alleles\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=STZ,Number=1,Type=Float,Description=\"Correlation between strand and alleles\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=NMZ,Number=1,Type=Float,Description=\"Correlation between mismatch counts per read and alleles\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=IOR,Number=1,Type=Float,Description=\"Inflated rate of observing of other alleles in log10 scale\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=NM0,Number=1,Type=Float,Description=\"Average number of mismatches in the reads with ref alleles\">\n");	
  bcf_hdr_append(odw->hdr, "##INFO=<ID=NM1,Number=1,Type=Float,Description=\"Average number of mismatches in the reads with non-ref alleles\">\n");
  bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
  bcf_hdr_append(odw->hdr, "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depth\">\n");
  bcf_hdr_append(odw->hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Allele Depth, including unidentifiable alleles\">\n");
  bcf_hdr_append(odw->hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scale Genotype Likelihoods\">\n");

  bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_snp,Description=\"Overlaps with snp\">");
  bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_indel,Description=\"Overlaps with indel\">");
  bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_vntr,Description=\"Overlaps with VNTR\">");  

  odw->write_hdr();
}
