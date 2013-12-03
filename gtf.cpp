/* The MIT License

   Copyright (c) 2013 Adrian Tan <atks@umich.edu>

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

#include "gtf.h"

class GTFRecord: public Interval
{
    public:
    std::string gene;
    std::string feature;
    std::string chrom;
    char strand;
    int32_t frame;
    int32_t exonNo;
    bool fivePrimeConservedEssentialSpliceSite;
    bool threePrimeConservedEssentialSpliceSite;
    bool containsStartCodon;
    bool containsStopCodon;
    uint32_t level;
    std::string attrib;

    GTFRecord(std::string& _chrom, uint32_t _start, uint32_t _end, char _strand,
              std::string& _gene, std::string& _feature, int32_t _frame, int32_t _exonNo,
              bool _fivePrimeConservedEssentialSpliceSite, bool _threePrimeConservedEssentialSpliceSite,
              bool _containsStartCodon, bool _containsStopCodon,
              uint32_t _level, std::string& _attrib)
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
        attrib = _attrib;
    };

    void print()
    {
        std::cerr << "chrom   : " << chrom << "\n";
        std::cerr << "[" << start << "," << end << "]\n";
        std::cerr << "strand                    : " << strand << "\n";
        std::cerr << "address                   : " << this << "\n";
        std::cerr << "gene                      : " << gene << "\n";
        std::cerr << "feature                   : " << feature << "\n";
        std::cerr << "frame                     : " << frame << "\n";
        std::cerr << "exon number               : " << exonNo << "\n";
        std::cerr << "5' conserved splice site  : " << fivePrimeConservedEssentialSpliceSite << "\n";
        std::cerr << "3' conserved splice site  : " << threePrimeConservedEssentialSpliceSite << "\n";
        std::cerr << "contains start codon      : " << containsStartCodon << "\n";
        std::cerr << "contains stop codon       : " << containsStopCodon << "\n";
        std::cerr << "level                     : " << level << "\n";
        std::cerr << "attrib                    : " << attrib << "\n";
    };

    private:
};

/**
 * Process
 */
void populate_gencode_tree(Igor& igor)
{
//    std::string line;
//  std::vector<std::string> fields;
//  std::vector<std::string> alleles;
//
//    std::clog << "Populating GENCODE tree ... \n";
//
//    std::vector<GTFRecord*> exons;
//    std::map<std::string, IntervalTree*> GENCODE;
//    std::vector<Interval*> intervals;
//    IntervalTree* G;
//    uint32_t noRecords = 0;
//
//  //for gencode records
//    boost::regex regex("chr(.+)");
//  boost::smatch regexResult;
//
//  bool allChromosomes = igor.chromosome == "all";
//
//  uint32_t noConservedSpliced = 0;
//  uint32_t noUnconservedSpliced = 0;
//
//  //construct tree
//    while (!ifeof(IN_GTF))
//  {
//      IN_GTF >> line;
//
//      //std::cerr << line ;
//
//      if(line.at(0)!='#')
//      {
//          split(fields, '\t', line, 9);
//            std::string chrom = fields[0];
//
//          if (boost::regex_match(chrom, regexResult, regex))
//          {
//              chrom = regexResult[1];
//          }
//
//            //std::cerr << "chrom: " << chrom << "\n";
//            //std::cerr << "chromosome: " << chromosome << "\n";
//
//            //if (chrom=="chr20" && feature=="CDS")
//            if (!allChromosomes && chrom!=chromosome)
//            {
//                continue;
//            }
//
//            std::string feature = fields[2];
//
//            if (feature!="exon" && feature!="CDS" && feature!="start_codon" && feature!="stop_codon")
//            {
//                continue;
//            }
//
//            //create tree for chromosome
//            if(!exists(GENCODE, chrom))
//            {
//                std::clog << "Creating tree for chromosome " << chrom << "\n";
//                GENCODE[chrom] = new IntervalTree();
//            }
//
//            //process fields
//            std::map<std::string, std::string> attrib_map;
//            splitGTFAttributeFields(attrib_map, fields[8]);
//
////        std::cerr << line << "\n";
////            for (std::map<std::string, std::string>::iterator i=attrib_map.begin(); i!=attrib_map.end() ; ++i)
////            {
////                std::cerr << i->first << " => " << i->second << "\n";
////            }
//
//          uint32_t start = boost::lexical_cast<int32_t>(fields[3]);
//          uint32_t end = boost::lexical_cast<int32_t>(fields[4]);
//          char strand = fields[6].at(0);
//            std::string gene = attrib_map["gene_id"];
//          int32_t frame = fields[7]=="." ? -1 : boost::lexical_cast<int32_t>(fields[7]);
//            int32_t exonNo = !exists(attrib_map, "exon_number") ? -1 : boost::lexical_cast<int32_t>(attrib_map["exon_number"]);
//          uint32_t level = boost::lexical_cast<uint32_t>(attrib_map["level"]);
//          bool fivePrimeConservedEssentialSpliceSite = false;
//          bool threePrimeConservedEssentialSpliceSite = false;
//          bool containsStartCodon = false;
//          bool containsStopCodon = false;
//          std::string attrib = fields[8];
//
//          if (feature=="exon")
//          {
//              if (attrib_map["gene_type"]!="protein_coding")
//                {
//                  continue;
//              }
//
//                std::string dnc1;
//                std::string dnc2;
//                uint32_t chromNo = myRefSeq->getChromosome(chrom=="M"?"MT":chrom.c_str());
//
//
//                //std::cerr << chromNo << " " <<  chrom << "\n";
//
//                if(strand=='+')
//                {
//                    myRefSeq->getString(dnc1, chromNo, (genomeIndex_t) start-2, 2);
//                    myRefSeq->getString(dnc2, chromNo, (genomeIndex_t) end+1, 2);
//
//                    if (dnc1=="AG")
//                    {
//                        fivePrimeConservedEssentialSpliceSite = true;
//                    }
//
//                    if (dnc2=="GT")
//                    {
//                        threePrimeConservedEssentialSpliceSite = true;
//                    }
//                }
//
//                if(strand=='-')
//                {
//                    myRefSeq->getString(dnc2, chromNo, (genomeIndex_t) start-2, 2);
//                    myRefSeq->getString(dnc1, chromNo, (genomeIndex_t) end+1, 2);
//
//                    if (dnc1=="CT")
//                    {
//                        fivePrimeConservedEssentialSpliceSite = true;
//                    }
//
//                    if (dnc2=="AC")
//                    {
//                        threePrimeConservedEssentialSpliceSite = true;
//                    }
//                }
//          }
//
//            if (feature=="stop_codon")
//            {
//                if (exists(GENCODE, chrom))
//                {
//                    GENCODE[chrom]->search(start, end, intervals);
//                }
//
//                for (uint32_t i=0; i<intervals.size(); ++i)
//                {
//                    GTFRecord* record = (GTFRecord*)intervals[i];
//                    if (record->feature == "exon" && record->gene == gene)
//                    {
//                        record->containsStopCodon = true;
//                    }
//                }
//            }
//
//            if (feature=="start_codon")
//            {
//                if (exists(GENCODE, chrom))
//                {
//                    GENCODE[chrom]->search(start, end, intervals);
//                }
//
//                for (uint32_t i=0; i<intervals.size(); ++i)
//                {
//                    GTFRecord* record = (GTFRecord*)intervals[i];
//                    if (record->feature == "exon" && record->gene == gene)
//                    {
//                        record->containsStartCodon = true;
//                    }
//                }
//            }
//
//            GTFRecord* gtfRecord = new GTFRecord(chrom, start, end, strand,
//                                                 gene, feature, frame, exonNo,
//                                                 fivePrimeConservedEssentialSpliceSite, threePrimeConservedEssentialSpliceSite,
//                                                 containsStartCodon, containsStopCodon,
//                                                 level, attrib);
//
//            GENCODE[chrom]->insert(gtfRecord);
//
//            if (feature=="exon")
//            {
//                exons.push_back(gtfRecord);
//            }
////            GENCODE.print();
////            if (GENCODE.size()==200)
////            {
////               std::cerr << "size: " << GENCODE.size() << "\n";
////               break;
////            }
////            if (GENCODE.size()%1000==0 && GENCODE.size()!=0)
////            if (GENCODE.size()==100)
////            {
////               std::cerr << "size: " << GENCODE.size() << "\n";
////               break;
////            }
////
//            ++noRecords;
//        }
//  }
//
//  std::clog << " ... completed\n";
//  if (GENCODE.size()==0)
//  {
//      std::cerr << "No reference GENCODE features!\n";
//        exit(1);
//  }
//
//  std::cerr << "start validation\n";
//    GENCODE["20"]->validate();
//    //G = GENCODE["chr20"];
//    //G->validate();
//    std::cerr << "end validation\n";
//    std::cerr << "height : " << GENCODE["20"]->height << "\n";
//    std::cerr << "size : " << GENCODE["20"]->size() << "\n";
//    std::cerr << "height : " << G->height << "\n";
//    std::cerr << "size : " << G->size() << "\n";


//    for (uint32_t i=0; i<exons.size(); ++i)
//    {
//        if (((exons[i]->exonNo==1 || exons[i]->fivePrimeConservedEssentialSpliceSite) && (exons[i]->containsStopCodon ||  exons[i]->threePrimeConservedEssentialSpliceSite)))
//        {
//            ++noConservedSpliced;
//        }
//        else
//        {
//            ++noUnconservedSpliced;
//        }

//        char strand = exons[i]->strand;
//        uint32_t start = exons[i]->start;
//        uint32_t end = exons[i]->end;
//        uint32_t chromNo = myRefSeq->getChromosome(exons[i]->chrom.c_str());
//        std::string dnc1;
//        std::string dnc2;
//
//        if(strand=='+')
//        {
//            myRefSeq->getString(dnc1, chromNo, (genomeIndex_t) start-2, 2);
//            myRefSeq->getString(dnc2, chromNo, (genomeIndex_t) end+1, 2);
//
//            //if (dnc1=="AG" && dnc2=="GT")
//            //if (((exons[i]->containsStartCodon || dnc1=="AG") && (exons[i]->containsStopCodon || dnc2=="GT")) ||
//            //     (!exons[i]->containsStopCodon && !exons[i]->containsStartCodon && (dnc1=="AG" || dnc2=="GT")))
//            if (((exons[i]->exonNo==1 || dnc1=="AG") && (exons[i]->containsStopCodon || dnc2=="GT")))
//            {
//                ++noConservedSpliced;
//                //std::cerr << "+ conserved " << dnc1 << " " << dnc2 << "\n";
//            }
//            else
//            {
//                ++noUnconservedSpliced;
//                //std::cerr << "+ unconserved " << dnc1 << " " << dnc2 << "\n";
//                //exons[i]->print();
//            }
//        }
//
//        if(strand=='-')
//        {
//            myRefSeq->getString(dnc2, chromNo, (genomeIndex_t) start-2, 2);
//            myRefSeq->getString(dnc1, chromNo, (genomeIndex_t) end+1, 2);
//
//            //if (dnc1=="CT" && dnc2=="AC")
//            //if ((exons[i]->containsStartCodon || dnc1=="CT") &&
//            //    (exons[i]->containsStopCodon || dnc2=="AC"))
//            //if (((exons[i]->containsStartCodon || dnc1=="CT") && (exons[i]->containsStopCodon || dnc2=="AC")) ||
//            //     (!exons[i]->containsStopCodon && !exons[i]->containsStartCodon && (dnc1=="CT" || dnc2=="AC")))
//            if (((exons[i]->exonNo==1 || dnc1=="CT") && (exons[i]->containsStopCodon || dnc2=="AC")))
//            {
//                ++noConservedSpliced;
//                //std::cerr << "- conserved " << dnc1 << " " << dnc2 << "\n";
//            }
//            else
//            {
//                ++noUnconservedSpliced;
//                //std::cerr << "- unconserved " << dnc1 << " " << dnc2 << "\n";
//                //exons[i]->print();
//            }
//        }
//    }


}

bool exists(std::map<std::string, IntervalTree*>& map, const std::string& key)
{
    return map.end()!=map.find(key);
}

/**
Splits a line into a map - PERL style
*/
void splitGTFAttributeFields(std::map<std::string, std::string>& map, std::string& str)
{
    map.clear();
    const char* tempStr = str.c_str();
    int32_t i=0, lastIndex = str.size()-1;
    std::string key;
    std::string val;
    std::stringstream token;

    if (lastIndex<0) return;

    while (i<=lastIndex)
    {
        //read next character
        if(tempStr[i]!=';' && tempStr[i]!=' ')
        {
            token << tempStr[i];
        }

        //store key-value pair
        if (i==lastIndex || (tempStr[i]==';' && tempStr[i+1]==' '))
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

            if (i!=lastIndex)
                ++i;
        }

        //store key
        if (tempStr[i]==' ')
        {
            key = token.str();
            token.str("");
        }

        ++i;
    }
};
