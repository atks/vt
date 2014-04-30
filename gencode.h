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

#ifndef GENCODE_H
#define GENCODE_H

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
#include "htslib/tbx.h"
#include "hts_utils.h"
#include "utils.h"
#include "interval_tree.h"
#include "variant_manip.h"
#include "genome_interval.h"
#include "tbx_ordered_reader.h"

#define GC_FT_EXON 0
#define GC_FT_CDS  1
#define GC_FT_START_CODON 2
#define GC_FT_STOP_CODON 3

#define NT_N 0
#define NT_A 1
#define NT_C 2
#define NT_G 4
#define NT_T 8

#define ALA 0
#define ARG 1
#define ASN 2
#define ASP 3
#define CYS 4
#define GLN 5
#define GLU 6
#define GLY 7
#define HIS 8
#define ILE 9
#define LEU 10
#define LYS 11
#define MET 12
#define PHE 13
#define PRO 14
#define SER 15
#define THR 16
#define TRP 17
#define TYR 18
#define VAL 19

class GENCODERecord : public Interval
{
    public:
    std::string gene;
    int32_t feature;
    std::string chrom;
    char strand;
    int32_t frame;
    int32_t *syn; // synonymous array for sequence, defined only if CDS.
    int32_t exonNo;
    bool fivePrimeConservedEssentialSpliceSite;
    bool threePrimeConservedEssentialSpliceSite;
    bool containsStartCodon;
    bool containsStopCodon;
    int32_t level;

    GENCODERecord(std::string& chrom, int32_t start, int32_t end, char strand,
                  std::string& gene, int32_t feature, int32_t frame, int32_t exonNo,
                  bool fivePrimeConservedEssentialSpliceSite, bool threePrimeConservedEssentialSpliceSite,
                  bool containsStartCodon, bool containsStopCodon,
                  int32_t level);

    /**
     * Checks if base at position position is synonymous.
     */
    bool is_synonymous(int32_t pos1, char base);

    /**
     * Prints this GENCODE record to STDERR.
     */
    void print();

    /**
     * Converts feature to string.
     */
    void feature2string(int32_t feature, kstring_t *s);

    private:
};

KHASH_MAP_INIT_STR(aadict, int32_t)

class GENCODE
{
    public:
    std::string gencode_gtf_file;
    std::string ref_fasta_file;
    faidx_t *fai;
    std::map<std::string, IntervalTree*> CHROM;
    std::stringstream token;
    khash_t(aadict) *codon2syn;

    /**
     * Constructs and initialize a GENCODE object.
     */
    GENCODE(std::string& gencode_gtf_file, std::string& ref_fasta_file, std::vector<GenomeInterval>& intervals);

    /**
     * Constructs a GENCODE object.
     */
    GENCODE(std::string& gencode_gtf_file, std::string& ref_fasta_file);

    /**
     * Initialize a vector of intervals.
     */
    void initialize(std::vector<GenomeInterval>& intervals);

    /**
     * Initialize a chromosome in the GENCODE tree.
     */
    void initialize(std::string& chrom);

    /**
     * Gets overlapping intervals with chrom:start1-end1.
     */
    void search(std::string& chrom, int32_t start1, int32_t end1, std::vector<Interval*>& intervals);

    /**
     * Splits a line into a map - PERL style.
     */
    void split_gtf_attribute_field(std::map<std::string, std::string>& map, std::string& str);

    /**
     * Generate array for ease of checking synonymous, non synonymous SNPs.
     */
    void fill_synonymous(GENCODERecord *g);

    /**
     * Parses a string to a GENCODERecord
     */
    GENCODERecord* parse_gencode(kstring_t *s)
    {
//        std::vector<std::string> fields;
//        std::map<std::string, std::string> attrib_map;
//                
//        split(fields, "\t", s.s);
//
//        std::string chrom = fields[0]=="M" ? std::string("MT") : fields[0];
//        std::string& feature = fields[2];
//        int32_t gencode_feature = -1;
//
//        if (feature!="exon" && feature!="CDS" && feature!="start_codon" && feature!="stop_codon")
//        {
//            continue;
//        }
//
//        //process fields
//        std::string& attrib = fields[8];
//        split_gtf_attribute_field(attrib_map, attrib);
//
//        int32_t start1; str2int32(fields[3], start1);
//        int32_t end1; str2int32(fields[4], end1);
//        char strand = fields[6].at(0);
//        std::string& gene = attrib_map["gene_name"];
//
//        int32_t frame;
//        if (!str2int32(fields[7], frame)) frame = -1;
//
//        int32_t level;
//        if (!str2int32(attrib_map["level"], level))
//        {
//            level = -1;
//        }
//
//        bool fivePrimeConservedEssentialSpliceSite = false;
//        bool threePrimeConservedEssentialSpliceSite = false;
//        bool containsStartCodon = false;
//        bool containsStopCodon = false;
//
//        int32_t exon_no = -1;
//
//        if (feature=="CDS")
//        {
//            gencode_feature = GC_FT_CDS;
//        }
//
//        if (feature=="exon")
//        {
//            gencode_feature = GC_FT_EXON;
//
//            if(attrib_map.find("exon_number")!=attrib_map.end() && !str2int32(attrib_map["exon_number"], exon_no))
//            {
//                exon_no = -1;
//            }
//
//            int32_t ref_len1 = 0;
//            int32_t ref_len2 = 0;
//
//            char *dnc1 = faidx_fetch_seq(fai, chrom.c_str(), start1-3, start1-2, &ref_len1);
//            char *dnc2 = faidx_fetch_seq(fai, chrom.c_str(), end1, end1+1, &ref_len2);
//
//            if(strand=='+')
//            {
//                if (!strcmp(dnc1,"AG"))
//                {
//                    fivePrimeConservedEssentialSpliceSite = true;
//                }
//
//                if (!strcmp(dnc2,"GT"))
//                {
//                    threePrimeConservedEssentialSpliceSite = true;
//                }
//            }
//
//            if(strand=='-')
//            {
//                if (!strcmp(dnc2,"CT"))
//                {
//                    fivePrimeConservedEssentialSpliceSite = true;
//                }
//
//                if (!strcmp(dnc1,"AC"))
//                {
//                    threePrimeConservedEssentialSpliceSite = true;
//                }
//            }
//
//            if (ref_len1) free(dnc1);
//            if (ref_len2) free(dnc2);
//        }
//
//        if (feature=="stop_codon")
//        {
//            gencode_feature = GC_FT_START_CODON;
//
//            CHROM[chrom]->search(start1, end1, overlaps);
//
//            for (uint32_t i=0; i<overlaps.size(); ++i)
//            {
//                GENCODERecord* record = (GENCODERecord*)overlaps[i];
//                if (record->feature == GC_FT_EXON && record->gene == gene)
//                {
//                    record->containsStopCodon = true;
//                }
//            }
//        }
//
//        if (feature=="start_codon")
//        {
//            gencode_feature = GC_FT_STOP_CODON;
//
//            CHROM[chrom]->search(start1, end1, overlaps);
//
//            for (uint32_t i=0; i<overlaps.size(); ++i)
//            {
//                GENCODERecord* record = (GENCODERecord*)overlaps[i];
//                if (record->feature == GC_FT_EXON && record->gene == gene)
//                {
//                    record->containsStartCodon = true;
//                }
//            }
//        }
//        
//        GENCODERecord* record = new GENCODERecord(chrom, start1, end1, strand,
//                                             gene, gencode_feature, frame, exon_no,
//                                             fivePrimeConservedEssentialSpliceSite, threePrimeConservedEssentialSpliceSite,
//                                             containsStartCodon, containsStopCodon,
//                                             level);
        return NULL;
    }
    
    /**
     * Checks if a record is overlapping.  This uses the orderedness of the gtf file to search through the records.
     */
    bool overlaps_with(std::string& chrom, int32_t start1, int32_t end1, int32_t feature_type)
    {
        bool search_flag = false;
        bool fill_buffer = true;
        
        //is this the right chromosome?
        if (current_chrom!=chrom)
        {
            GenomeInterval interval(chrom);
            todr->jump_to_interval(interval);
        }

        std::list<GENCODERecord*>::iterator i;
        while ( i!=buffer.end())
        {
            if ((*i)->end<start1)
            {
                delete *i;
                i = buffer.erase(i);
            }    
            else if ((*i)->overlaps_with(start1, end1))   
            {
                if ((*i)->feature==feature_type)
                {
                    search_flag=true;
                }
            }
            else //after end1
            {
                fill_buffer = false;
                break;
            }
        }
        
        //read and store records till record is after target
        if (fill_buffer)
        {
            kstring_t s = {0,0,0};
            
            while (todr->read(&s))
            {
                GENCODERecord *rec = parse_gencode(&s);
                buffer.push_back(rec);
                
                if (rec->start>end1)
                {
                    break;
                }    
            }
            
            if (s.m) free(s.s);
        }    
        
        return search_flag;
    };

    private:
    std::string current_chrom;
    std::list<GENCODERecord*> buffer;
    TBXOrderedReader *todr;
};

#endif