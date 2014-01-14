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

#include "gencode.h"

GENCODERecord::GENCODERecord(std::string& _chrom, uint32_t _start, uint32_t _end, char _strand,
          std::string& _gene, int32_t _feature, int32_t _frame, int32_t _exonNo,
          bool _fivePrimeConservedEssentialSpliceSite, bool _threePrimeConservedEssentialSpliceSite,
          bool _containsStartCodon, bool _containsStopCodon,
          uint32_t _level)
{
    chrom = _chrom;
    start = _start;
    end = _end;
    strand = _strand;
    gene = _gene;
    feature = _feature;
    frame = _frame;
    exonNo = _exonNo;
    fivePrimeConservedEssentialSpliceSite = _fivePrimeConservedEssentialSpliceSite;
    threePrimeConservedEssentialSpliceSite = _threePrimeConservedEssentialSpliceSite;
    containsStartCodon = _containsStartCodon;
    containsStopCodon = _containsStopCodon;
    level = _level;
};

/**
 * Prints this GENCODE record to STDERR.
 */
void GENCODERecord::print()
{
    std::cerr << "chrom   : " << chrom << ":" << start << "-" << end << "\n";
    std::cerr << "strand                    : " << strand << "\n";
    std::cerr << "gene                      : " << gene << "\n";
    kstring_t s = {0,0,0};
    feature2string(feature, &s);
    std::cerr << "feature                   : " << s.s << "\n";
    std::cerr << "frame                     : " << frame << "\n";
    std::cerr << "exon number               : " << exonNo << "\n";
    std::cerr << "5' conserved splice site  : " << fivePrimeConservedEssentialSpliceSite << "\n";
    std::cerr << "3' conserved splice site  : " << threePrimeConservedEssentialSpliceSite << "\n";
    std::cerr << "contains start codon      : " << containsStartCodon << "\n";
    std::cerr << "contains stop codon       : " << containsStopCodon << "\n";
    std::cerr << "level                     : " << level << "\n";
    if (s.m) free(s.s);
};

/**
 * Converts feature to string.
 */
void GENCODERecord::feature2string(int32_t feature, kstring_t *s)
{
    s->l = 0;

    if (feature==GC_FT_EXON)
    {
        kputs("exon", s);
    }
    else if (feature==GC_FT_CDS)
    {
        kputs("CDS", s);
    }
    else if (feature==GC_FT_START_CODON)
    {
        kputs("start_codon", s);
    }
    else if (feature==GC_FT_STOP_CODON)
    {
        kputs("stop_codon", s);
    }
}

/**
 * Constructs and initialized a GENCODE object.
 */
GENCODE::GENCODE(std::string& gencode_gtf_file, std::string& ref_fasta_file, std::vector<GenomeInterval>& intervals)
{
    fai = fai_load(ref_fasta_file.c_str());
    this->gencode_gtf_file = gencode_gtf_file;
    initialize(intervals);
}
    
/**
 * Constructs a GENCODE object.
 */
GENCODE::GENCODE(std::string& gencode_gtf_file, std::string& ref_fasta_file)
{
    fai = fai_load(ref_fasta_file.c_str());
    this->gencode_gtf_file = gencode_gtf_file;
}    
    
/**
 * Initialize a vector of intervals.
 */
void GENCODE::initialize(std::vector<GenomeInterval>& intervals)
{
    std::vector<GenomeInterval> chromosomes;
    for (int32_t i=0; i<intervals.size(); ++i)
    {
        intervals[i].chromosomify();
        std::string chrom = intervals[i].to_string();
        if (CHROM.find(chrom)==CHROM.end())
        {
            std::clog << "Initializing GENCODE tree for chromosome " << chrom << "\n";
            CHROM[chrom] = new IntervalTree();
            chromosomes.push_back(intervals[i]);
        }
        
    }    
    
    TBXOrderedReader *todr = new TBXOrderedReader(gencode_gtf_file, chromosomes);
    std::vector<std::string> fields;
    std::map<std::string, std::string> attrib_map;
    kstring_t s = {0,0,0};

    //for storing returned overlapping intervals
    std::vector<Interval*> overlaps;
    
    while (todr->read(&s))
    {
        //populate interval trees with reference sets
        //gene transfer format
        //1   chr1
        //2   HAVANA - source
        //3   exon   - feature
        //4   13221  - start
        //5   14409  - end
        //6   .      - score
        //7   +      - strand
        //8   .      - frame
        //9   gene_id "ENSG - attribute
            //gene_id "ENSG00000149656.4";
            //transcript_id "ENST00000425473.1";
            //gene_type "processed_transcript";
            //gene_status "KNOWN";
            //gene_name "LINC00266-1";
            //transcript_type "processed_transcript";
            //transcript_status "KNOWN";
            //transcript_name "LINC00266-1-002"; exon_number 3;
            //exon_id "ENSE00002450675.1";
            //level 2; havana_gene "OTTHUMG00000033036.2";
            //havana_transcript "OTTHUMT00000080305.1";

        split(fields, "\t", s.s);

        std::string chrom = fields[0]=="M" ? std::string("MT") : fields[0];
        std::string& feature = fields[2];
        int32_t gencode_feature = GC_FT_CDS;    

        if (feature!="exon" && feature!="CDS" && feature!="start_codon" && feature!="stop_codon")
        {
            continue;
        }
        
        //process fields
        std::string& attrib = fields[8];
        split_gtf_attribute_field(attrib_map, attrib);

        int32_t start1; str2int32(fields[3], start1);
        int32_t end1; str2int32(fields[4], end1);
        char strand = fields[6].at(0);
        std::string& gene = attrib_map["gene_name"];

        int32_t frame;
        if (!str2int32(fields[7], frame)) frame = -1;

        int32_t level;
        if (!str2int32(attrib_map["level"], level))
        {
            level = -1;
        }

        bool fivePrimeConservedEssentialSpliceSite = false;
        bool threePrimeConservedEssentialSpliceSite = false;
        bool containsStartCodon = false;
        bool containsStopCodon = false;


        int32_t exon_no = -1;
        
        if (feature=="exon")
        {
            gencode_feature = GC_FT_EXON;

            if(attrib_map.find("exon_number")!=attrib_map.end() && !str2int32(attrib_map["exon_number"], exon_no))
            {
                exon_no = -1;
            }

            int32_t ref_len1 = 0;
            int32_t ref_len2 = 0;
            
            char *dnc1 = faidx_fetch_seq(fai, chrom.c_str(), start1-3, start1-2, &ref_len1);
            char *dnc2 = faidx_fetch_seq(fai, chrom.c_str(), end1, end1+1, &ref_len2);
            
            if(strand=='+')
            {
                if (!strcmp(dnc1,"AG"))
                {
                    fivePrimeConservedEssentialSpliceSite = true;
                }

                if (!strcmp(dnc2,"GT"))
                {
                    threePrimeConservedEssentialSpliceSite = true;
                }
            }

            if(strand=='-')
            {
                if (!strcmp(dnc2,"CT"))
                {
                    fivePrimeConservedEssentialSpliceSite = true;
                }

                if (!strcmp(dnc1,"AC"))
                {
                    threePrimeConservedEssentialSpliceSite = true;
                }
            }

            if (ref_len1) free(dnc1);
            if (ref_len2) free(dnc2);

        }

        if (feature=="stop_codon")
        {
            gencode_feature = GC_FT_START_CODON;
            
            CHROM[chrom]->search(start1, end1, overlaps);

            for (uint32_t i=0; i<overlaps.size(); ++i)
            {
                GENCODERecord* record = (GENCODERecord*)overlaps[i];
                if (record->feature == GC_FT_EXON && record->gene == gene)
                {
                    record->containsStopCodon = true;
                }
            }
        }

        if (feature=="start_codon")
        {
            gencode_feature = GC_FT_STOP_CODON;
            
            CHROM[chrom]->search(start1, end1, overlaps);

            for (uint32_t i=0; i<overlaps.size(); ++i)
            {
                GENCODERecord* record = (GENCODERecord*)overlaps[i];
                if (record->feature == GC_FT_EXON && record->gene == gene)
                {
                    record->containsStartCodon = true;
                }
            }
        }

        GENCODERecord* record = new GENCODERecord(chrom, start1, end1, strand,
                                             gene, gencode_feature, frame, exon_no,
                                             fivePrimeConservedEssentialSpliceSite, threePrimeConservedEssentialSpliceSite,
                                             containsStartCodon, containsStopCodon,
                                             level);

        CHROM[chrom]->insert(record);
    }
    
    todr->close();
}

/**
 * Initialize a chromosome in the GENCODE tree.
 */
void GENCODE::initialize(std::string& chrom)
{
    std::vector<GenomeInterval> intervals;
    intervals.push_back(GenomeInterval(chrom));
    initialize(intervals);
}

/**
 * Gets overlapping intervals with chrom:start1-end1.
 */
void GENCODE::search(std::string& chrom, int32_t start1, int32_t end1, std::vector<Interval*>& intervals)
{
    if (CHROM.find(chrom)==CHROM.end())
    {
        initialize(chrom);
    }
    
    CHROM[chrom]->search(start1, end1, intervals);
}

/**
 * Splits a line into a map - PERL style.
 */
void GENCODE::split_gtf_attribute_field(std::map<std::string, std::string>& map, std::string& str)
{
    map.clear();
    const char* tempStr = str.c_str();
    int32_t i=0, lastIndex = str.size()-1;
    std::string key;
    std::string val;
    
    if (lastIndex<0) return;

    token.str("");
    while (i<=lastIndex)
    {
        //read next character
        if(tempStr[i]!=';' && tempStr[i]!=' ')
        {
            token << tempStr[i];
        }

        //store key-value pair
        if (i==lastIndex || (tempStr[i]==';'))
        {
            val = token.str();
            if (val.at(0)=='"')
            {
                val = token.str();
                val = val.substr(1,val.size()-2);
            }
            map[key] = val;

            token.str("");
            key.clear();
        }

        //store key
        if (tempStr[i]==' ' && token.str().size()!=0)
        {
            key = token.str();
            token.str("");
        }

        ++i;
    }
};