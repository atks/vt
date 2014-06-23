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

GENCODERecord::GENCODERecord(std::string& chrom, int32_t start, int32_t end, char strand,
          std::string& gene, int32_t feature, int32_t frame, int32_t exonNo,
          bool fivePrimeConservedEssentialSpliceSite, bool threePrimeConservedEssentialSpliceSite,
          bool containsStartCodon, bool containsStopCodon,
          int32_t level)
{
    this->chrom = chrom;
    this->start = start;
    this->end = end;
    this->strand = strand;
    this->gene = gene;
    this->feature = feature;
    this->frame = frame;
    this->exonNo = exonNo;
    this->fivePrimeConservedEssentialSpliceSite = fivePrimeConservedEssentialSpliceSite;
    this->threePrimeConservedEssentialSpliceSite = threePrimeConservedEssentialSpliceSite;
    this->containsStartCodon = containsStartCodon;
    this->containsStopCodon = containsStopCodon;
    this->level = level;
};

/**
 * Checks if base at position position is synonymous.
 */
bool is_synonymous(int32_t pos1, char base)
{   
    return true;
}

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
    if (fai==NULL) 
    {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }
    this->gencode_gtf_file = gencode_gtf_file;
    initialize(intervals);
    
    khiter_t k;
    int32_t ret;
    codon2syn = kh_init(aadict);    
    
    //constructs the mapping for amino acids
    //ALA
    k = kh_put(aadict, codon2syn, "GCA", &ret); 
    kh_value(codon2syn, k) = (NT_G<<8) & (NT_C<<4) & (NT_A|NT_C|NT_G|NT_T);
    k = kh_put(aadict, codon2syn, "GCC", &ret); 
    kh_value(codon2syn, k) = (NT_G<<8) & (NT_C<<4) & (NT_A|NT_C|NT_G|NT_T);
    k = kh_put(aadict, codon2syn, "GCG", &ret); 
    kh_value(codon2syn, k) = (NT_G<<8) & (NT_C<<4) & (NT_A|NT_C|NT_G|NT_T);
    k = kh_put(aadict, codon2syn, "GCT", &ret); 
    kh_value(codon2syn, k) = (NT_G<<8) & (NT_C<<4) & (NT_A|NT_C|NT_G|NT_T);
    //ALA
    k = kh_put(aadict, codon2syn, "CGA", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_C)<<8) & (NT_G<<4) & (NT_A|NT_C|NT_G|NT_T);
    k = kh_put(aadict, codon2syn, "CGC", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_C)<<8) & (NT_G<<4) & (NT_A|NT_C|NT_G|NT_T);    
    k = kh_put(aadict, codon2syn, "CGG", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_C)<<8) & (NT_G<<4) & (NT_A|NT_C|NT_G|NT_T);    
    k = kh_put(aadict, codon2syn, "CGT", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_C)<<8) & (NT_G<<4) & (NT_A|NT_C|NT_G|NT_T);
    k = kh_put(aadict, codon2syn, "AGA", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_C)<<8) & (NT_G<<4) & (NT_A|NT_C|NT_G|NT_T);    
    k = kh_put(aadict, codon2syn, "AGG", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_C)<<8) & (NT_G<<4) & (NT_A|NT_C|NT_G|NT_T);
    //ASN
    k = kh_put(aadict, codon2syn, "AAC", &ret); 
    kh_value(codon2syn, k) = ((NT_A)<<8) & ((NT_A)<<4) & (NT_C|NT_T);
    k = kh_put(aadict, codon2syn, "AAT", &ret); 
    kh_value(codon2syn, k) = ((NT_A)<<8) & ((NT_A)<<4) & (NT_C|NT_T);
    //ASP
    k = kh_put(aadict, codon2syn, "GAC", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_A)<<4) & (NT_C|NT_T);
    k = kh_put(aadict, codon2syn, "GAT", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_A)<<4) & (NT_C|NT_T);
    //CYS
    k = kh_put(aadict, codon2syn, "TGT", &ret); 
    kh_value(codon2syn, k) = ((NT_T)<<8) & ((NT_G)<<4) & (NT_C|NT_T);
    k = kh_put(aadict, codon2syn, "TGC", &ret); 
    kh_value(codon2syn, k) = ((NT_T)<<8) & ((NT_G)<<4) & (NT_C|NT_T);       
    //GLN
    k = kh_put(aadict, codon2syn, "CAA", &ret); 
    kh_value(codon2syn, k) = ((NT_C)<<8) & ((NT_A)<<4) & (NT_A|NT_G);
    k = kh_put(aadict, codon2syn, "CAG", &ret); 
    kh_value(codon2syn, k) = ((NT_C)<<8) & ((NT_A)<<4) & (NT_A|NT_G);
    //GLU
    k = kh_put(aadict, codon2syn, "GAA", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_A)<<4) & (NT_A|NT_G);
    k = kh_put(aadict, codon2syn, "GAG", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_A)<<4) & (NT_A|NT_G);
    //GLY
    k = kh_put(aadict, codon2syn, "GGA", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_G)<<4) & (NT_A|NT_C|NT_G|NT_T);
    k = kh_put(aadict, codon2syn, "GGC", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_G)<<4) & (NT_A|NT_C|NT_G|NT_T);
    k = kh_put(aadict, codon2syn, "GGG", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_G)<<4) & (NT_A|NT_C|NT_G|NT_T);
    k = kh_put(aadict, codon2syn, "GGT", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_G)<<4) & (NT_A|NT_C|NT_G|NT_T);        
    //HIS
    k = kh_put(aadict, codon2syn, "CAT", &ret); 
    kh_value(codon2syn, k) = ((NT_C)<<8) & ((NT_A)<<4) & (NT_C|NT_T);
    k = kh_put(aadict, codon2syn, "CAC", &ret); 
    kh_value(codon2syn, k) = ((NT_C)<<8) & ((NT_A)<<4) & (NT_C|NT_T);
    //ILE
    k = kh_put(aadict, codon2syn, "ATT", &ret); 
    kh_value(codon2syn, k) = ((NT_A)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_T);    
    k = kh_put(aadict, codon2syn, "ATC", &ret); 
    kh_value(codon2syn, k) = ((NT_A)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_T);   
    k = kh_put(aadict, codon2syn, "ATA", &ret); 
    kh_value(codon2syn, k) = ((NT_A)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_T); 
    //LEU
    k = kh_put(aadict, codon2syn, "TTA", &ret); 
    kh_value(codon2syn, k) = ((NT_C|NT_T)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_G|NT_T); 
    k = kh_put(aadict, codon2syn, "TTG", &ret); 
    kh_value(codon2syn, k) = ((NT_C|NT_T)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_G|NT_T); 
    k = kh_put(aadict, codon2syn, "CTT", &ret); 
    kh_value(codon2syn, k) = ((NT_C|NT_T)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_G|NT_T); 
    k = kh_put(aadict, codon2syn, "CTC", &ret); 
    kh_value(codon2syn, k) = ((NT_C|NT_T)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_G|NT_T); 
    k = kh_put(aadict, codon2syn, "CTA", &ret); 
    kh_value(codon2syn, k) = ((NT_C|NT_T)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_G|NT_T); 
    k = kh_put(aadict, codon2syn, "CTG", &ret); 
    kh_value(codon2syn, k) = ((NT_C|NT_T)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_G|NT_T); 
    //LYS
    k = kh_put(aadict, codon2syn, "AAA", &ret); 
    kh_value(codon2syn, k) = ((NT_C|NT_A)<<8) & ((NT_A)<<4) & (NT_A|NT_G);     
    k = kh_put(aadict, codon2syn, "AAG", &ret); 
    kh_value(codon2syn, k) = ((NT_C|NT_A)<<8) & ((NT_A)<<4) & (NT_A|NT_G);      
    //MET
    k = kh_put(aadict, codon2syn, "AUG", &ret); 
    kh_value(codon2syn, k) = ((NT_C|NT_A)<<8) & ((NT_T)<<4) & (NT_G);  
    //PHE
    k = kh_put(aadict, codon2syn, "TTT", &ret); 
    kh_value(codon2syn, k) = ((NT_T)<<8) & ((NT_T)<<4) & (NT_T|NT_C);  
    k = kh_put(aadict, codon2syn, "TTC", &ret); 
    kh_value(codon2syn, k) = ((NT_T)<<8) & ((NT_T)<<4) & (NT_T|NT_C);  
    //PRO
    k = kh_put(aadict, codon2syn, "CCA", &ret); 
    kh_value(codon2syn, k) = ((NT_C)<<8) & ((NT_C)<<4) & (NT_A|NT_C|NT_G|NT_T);  
    k = kh_put(aadict, codon2syn, "CCC", &ret); 
    kh_value(codon2syn, k) = ((NT_C)<<8) & ((NT_C)<<4) & (NT_A|NT_C|NT_G|NT_T);  
    k = kh_put(aadict, codon2syn, "CCT", &ret); 
    kh_value(codon2syn, k) = ((NT_C)<<8) & ((NT_C)<<4) & (NT_A|NT_C|NT_G|NT_T);  
    k = kh_put(aadict, codon2syn, "CCG", &ret); 
    kh_value(codon2syn, k) = ((NT_C)<<8) & ((NT_C)<<4) & (NT_A|NT_C|NT_G|NT_T);          
    //SER
    k = kh_put(aadict, codon2syn, "TCT", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_T)<<8) & ((NT_C|NT_G)<<4) & (NT_A|NT_C|NT_G|NT_T);  
    k = kh_put(aadict, codon2syn, "TCC", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_T)<<8) & ((NT_C|NT_G)<<4) & (NT_A|NT_C|NT_G|NT_T);  
    k = kh_put(aadict, codon2syn, "TCA", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_T)<<8) & ((NT_C|NT_G)<<4) & (NT_A|NT_C|NT_G|NT_T);      
    k = kh_put(aadict, codon2syn, "TCG", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_T)<<8) & ((NT_C|NT_G)<<4) & (NT_A|NT_C|NT_G|NT_T);  
    k = kh_put(aadict, codon2syn, "AGT", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_T)<<8) & ((NT_C|NT_G)<<4) & (NT_A|NT_C|NT_G|NT_T);  
    k = kh_put(aadict, codon2syn, "AGC", &ret); 
    kh_value(codon2syn, k) = ((NT_A|NT_T)<<8) & ((NT_C|NT_G)<<4) & (NT_A|NT_C|NT_G|NT_T);   
    //THR
    k = kh_put(aadict, codon2syn, "ACA", &ret); 
    kh_value(codon2syn, k) = ((NT_A)<<8) & ((NT_C)<<4) & (NT_A|NT_C|NT_G|NT_T);  
    k = kh_put(aadict, codon2syn, "ACC", &ret); 
    kh_value(codon2syn, k) = ((NT_A)<<8) & ((NT_C)<<4) & (NT_A|NT_C|NT_G|NT_T);  
    k = kh_put(aadict, codon2syn, "ACG", &ret); 
    kh_value(codon2syn, k) = ((NT_A)<<8) & ((NT_C)<<4) & (NT_A|NT_C|NT_G|NT_T);             
    k = kh_put(aadict, codon2syn, "ACT", &ret); 
    kh_value(codon2syn, k) = ((NT_A)<<8) & ((NT_C)<<4) & (NT_A|NT_C|NT_G|NT_T);  
    //TRP
    k = kh_put(aadict, codon2syn, "TGG", &ret); 
    kh_value(codon2syn, k) = ((NT_T)<<8) & ((NT_G)<<4) & (NT_G); 
    //TYR
    k = kh_put(aadict, codon2syn, "TAT", &ret); 
    kh_value(codon2syn, k) = ((NT_T)<<8) & ((NT_A)<<4) & (NT_C|NT_T); 
    k = kh_put(aadict, codon2syn, "TAC", &ret); 
    kh_value(codon2syn, k) = ((NT_T)<<8) & ((NT_A)<<4) & (NT_C|NT_T); 
    //VAL
    k = kh_put(aadict, codon2syn, "GTA", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_G|NT_T); 
    k = kh_put(aadict, codon2syn, "GTC", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_G|NT_T); 
    k = kh_put(aadict, codon2syn, "GTG", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_G|NT_T); 
    k = kh_put(aadict, codon2syn, "GTT", &ret); 
    kh_value(codon2syn, k) = ((NT_G)<<8) & ((NT_T)<<4) & (NT_A|NT_C|NT_G|NT_T);             
    //STOP
    k = kh_put(aadict, codon2syn, "TAA", &ret); 
    kh_value(codon2syn, k) = ((NT_T)<<8) & ((NT_A|NT_G)<<4) & (NT_A|NT_G);     
    k = kh_put(aadict, codon2syn, "TGA", &ret); 
    kh_value(codon2syn, k) = ((NT_T)<<8) & ((NT_A|NT_G)<<4) & (NT_A|NT_G);  
    k = kh_put(aadict, codon2syn, "TAG", &ret); 
    kh_value(codon2syn, k) = ((NT_T)<<8) & ((NT_A|NT_G)<<4) & (NT_A|NT_G);                
}

/**
 * Constructs a GENCODE object.
 */
GENCODE::GENCODE(std::string& gencode_gtf_file, std::string& ref_fasta_file)
{
    fai = fai_load(ref_fasta_file.c_str());
    if (fai==NULL) 
    {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }
    this->gencode_gtf_file = gencode_gtf_file;
}

/**
 * Generate array for ease of checking synonymous, non synonymous SNPs.
 */
void GENCODE::fill_synonymous(GENCODERecord *g)
{
    if (g->feature==GC_FT_CDS)
    {
        //extract sequence
        int32_t ref_len;
//        char* seq = faidx_fetch_seq(fai, g->chrom.c_str(), g->start, g->end, &ref_len);
        
//        g->syn = new int32_t[ref_len];
//        kstring_t s = {0,0,0}; 
//        //check for each codon, check the nucleotides for each 
//        //frame that will not induce a non synonymous amino acid
//        for (int32_t i=g->frame; i<ref_len; i+=3)
//        {
//            //get the 3 bases
//            s.l=0;
//            kputc(seq[i], &s);
//            kputc(seq[i+1], &s);
//            kputc(seq[i+2], &s);
//            
//            //access the syn value 
//           // int32_t val = kh_get(hm;
//            int32_t val = 0;
//            
//            //populate syn
//            g->syn[i] = (val >> 8) & 15;
//            g->syn[i+1] = (val >> 4) & 15;
//            g->syn[i+2] = val & 15;
//        }
        
//        free(seq);
        //if (s.m) free(s.s);
    }    
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
            std::clog << "Initializing GENCODE tree for chromosome " << chrom << " ... ";
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
        int32_t gencode_feature = -1;

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

        if (feature=="CDS")
        {
            gencode_feature = GC_FT_CDS;
        }

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

        if (gencode_feature == GC_FT_CDS)
        {
            fill_synonymous(record);
        }

        CHROM[chrom]->insert(record);
    }

    std::clog << " done.\n";
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