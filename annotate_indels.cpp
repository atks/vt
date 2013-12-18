#include "annotate_indels.h"

namespace
{

class BEDRecord: public Interval
{
    public:
    std::string chrom;

    BEDRecord(std::string& _chrom, uint32_t _start, uint32_t _end)
    {
        chrom = _chrom;
        start = _start;
        end = _end;
    };

    void print()
    {
        std::cerr << "chrom   : " << chrom << "\n";
        std::cerr << "[" << start << "," << end << "]\n";
    };

    private:
};

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

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string ref_fasta_file;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
//  std::string gencode_gtf_file;
//  std::string indel_hotspot_bed_file;
//  std::string str_ref_bed_file;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    /////////////////////
    //Reference Regions//
    /////////////////////
//    std::map<std::string, IntervalTree*> GENCODE;
//    std::map<std::string, IntervalTree*> INDEL_HOTSPOTS;
//    std::map<std::string, IntervalTree*> STR_REF;

    /////////
    //stats//
    /////////
    uint32_t no_indels_annotated;

    ////////////////
    //common tools//
    ////////////////
    VariantManip *vm;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        try
        {
            std::string desc = "Annotate indels.\n\n";

            std::cerr << desc;

            std::string version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file("i", "input_vcf", "Input VCF file", true, "", "string", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "output_vcf", "Output VCF file [default: stdout]", false, "-", "string", cmd);
//          TCLAP::ValueArg<std::string> arg_gencode_gtf_file("d", "gencode", "GENCODE GTF file", false, "/net/fantasia/home/atks/ref/encode/gencode.v15.annotation.gtf.gz", "string", cmd);
//          TCLAP::ValueArg<std::string> arg_indel_hotspot_bed_file("t", "hotspot", "Indel Hotspot (Gerton Lunter) BED file", false, "/net/fantasia/home/atks/ref/indel/human_g1k_v37_indelhotspots_20121212.bed.gz", "string", cmd);
//          TCLAP::ValueArg<std::string> arg_str_ref_bed_file("s", "str_ref", "STR Reference (lobSTR) BED file", false, "/net/fantasia/home/atks/ref/str/lobSTR_index_hg19_v2.0.3.bed.gz", "string", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("g", "ref_genome", "Reference Genome Fasta file", false, "/net/fantasia/home/atks/ref/genome/hs37d5.fa", "string", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
//          gencode_gtf_file = arg_gencode_gtf_file.getValue();
//          indel_hotspot_bed_file = arg_indel_hotspot_bed_file.getValue();
//          str_ref_bed_file = arg_str_ref_bed_file.getValue();
//          ref_fasta_file   = arg_ref_fasta_file.getValue();
//          regions = argRegions.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
            abort();
        }
    };

    ~Igor() {};

    void print_options()
    {
        std::clog << "annotate_indels v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file(s)     " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";

//        std::clog << "         GENCODE Annotation file        " << gencode_gtf_file << "\n";
//      std::clog << "         Indel Hotspot Annotation file  " << indel_hotspot_bed_file << "\n";
//      std::clog << "         STR Reference Annotation file  " << str_ref_bed_file << "\n";
//      std::clog << "         Reference Genome file          " << ref_fasta_file    << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: Total Number of Indels annotated     " << no_indels_annotated << "\n";
        std::clog << "\n";
    }

    void initialize()
    {
        //******************
        //i/o initialization
        //******************
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file);
        odw->link_hdr(odr->hdr);
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant Type - SNP, MNP, INDEL, HP, STR\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat unit in a STR or Homopolymer\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Repeat Length\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=NS,Number=0,Type=Flag,Description=\"Near to STR\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=FS,Number=0,Type=Flag,Description=\"Frameshift INDEL\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=NFS,Number=0,Type=Flag,Description=\"Non Frameshift INDEL\">");

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip(ref_fasta_file);

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_indels_annotated = 0;
    }

    void annotate_indels()
    {
//        populate_gencode_tree();
//        populate_str_ref_tree();
//        populate_indel_hotspots_tree();

        bcf1_t *v = odw->get_bcf1_from_pool();
        Variant variant;
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
            vm->detect_str(odr->hdr, v, variant);


            ++no_indels_annotated;
        

            v = odw->get_bcf1_from_pool();
        }
    };

    /**
     * Populate GENCODE tree
     */
//    void populate_gencode_tree()
//    {
//        std::string line;
//      std::vector<std::string> fields;
//      std::vector<std::string> alleles;
//
//        std::clog << "Populating GENCODE tree ... \n";
//
//        std::vector<GTFRecord*> exons;
//        std::vector<Interval*> intervals;
//        uint32_t noRecords = 0;
//
//      //for gencode records
//        boost::regex regex("chr(.+)");
//      boost::smatch regexResult;
//      uint32_t noConservedSpliced = 0;
//      uint32_t noUnconservedSpliced = 0;
//
//        tbx_t *tbx;
//      hts_itr_t *iter;
//      BGZF* gffgz = xbgzf_open(gencode_gtf_file.c_str(), "r");
//        if (!(tbx=tbx_index_load(gencode_gtf_file.c_str())))
//        {
//          std::cerr << "Fail to load the VCF GZ index\n";
//          abort();
//      }
//
//          std::stringstream region;
//      kstring_t s;
//        s.s = 0; s.l = s.m = 0;
//        for (int32_t i = 0; i < nseqs; ++i)
//      {
//            //choose records from chromosome i
//          region.str("");
//          region << "chr" << seqnames[i] ;
//          if (!(iter = tbx_itr_querys(tbx, region.str().c_str())))
//          {
//              continue;
//          }
//
//          region << "processing\n" ;
//
//          while (tbx_itr_next(gffgz, tbx, iter, &s) >= 0)
//            {
//
//                //populate interval trees with reference sets
//                //gene transfer format
//                //1   chr1
//                //2   HAVANA - source
//                //3   exon   - feature
//                //4   13221  - start
//                //5   14409  - end
//                //6   .      - score
//                //7   +      - strand
//                //8   .      - frame
//                //9   gene_id "ENSG - attribute
//
//                //std::cerr << s.s << "\n";
//
//                std::string line = std::string(s.s);
//
//              if(line.at(0)!='#')
//              {
//                  split(fields, '\t', line, 9);
//                    std::string chrom = fields[0];
//
//                  if (boost::regex_match(chrom, regexResult, regex))
//                  {
//                      chrom = regexResult[1];
//                  }
//
//                    std::string feature = fields[2];
//
//                    if (feature!="exon" && feature!="CDS" && feature!="start_codon" && feature!="stop_codon")
//                    {
//                        continue;
//                    }
//
//                    //create tree for chromosome
//                    if(!exists(GENCODE, chrom))
//                    {
//                        std::clog << "Creating tree for chromosome " << chrom << "\n";
//                        GENCODE[chrom] = new IntervalTree();
//                    }
//
//                    //process fields
//                    std::map<std::string, std::string> attrib_map;
//                    splitGTFAttributeFields(attrib_map, fields[8]);
//
////                    std::cerr << line << "\n";
////                    for (std::map<std::string, std::string>::iterator i=attrib_map.begin(); i!=attrib_map.end() ; ++i)
////                    {
////                        std::cerr << i->first << " => " << i->second << "\n";
////                    }
//                    //exit(1);
//                  uint32_t start = boost::lexical_cast<int32_t>(fields[3]);
//                  uint32_t end = boost::lexical_cast<int32_t>(fields[4]);
//                  char strand = fields[6].at(0);
//                    std::string gene = attrib_map["gene_id"];
//                  int32_t frame = fields[7]=="." ? -1 : boost::lexical_cast<int32_t>(fields[7]);
//                    int32_t exonNo = !existsInMap(attrib_map, "exon_number") ? -1 : boost::lexical_cast<int32_t>(attrib_map["exon_number"]);
//                  uint32_t level = boost::lexical_cast<uint32_t>(attrib_map["level"]);
//                  bool fivePrimeConservedEssentialSpliceSite = false;
//                  bool threePrimeConservedEssentialSpliceSite = false;
//                  bool containsStartCodon = false;
//                  bool containsStopCodon = false;
//                  std::string attrib = fields[8];
//
//                  if (feature=="exon")
//                  {
//                      if (attrib_map["gene_type"]!="protein_coding")
//                        {
//                          continue;
//                      }
//
//                        std::string dnc1;
//                        std::string dnc2;
//                        uint32_t chromNo = ref_seq->getChromosome(chrom=="M"?"MT":chrom.c_str());
//
//                        if(strand=='+')
//                        {
//                            ref_seq->getString(dnc1, chromNo, (genomeIndex_t) start-2, 2);
//                            ref_seq->getString(dnc2, chromNo, (genomeIndex_t) end+1, 2);
//
//                            if (dnc1=="AG")
//                            {
//                                fivePrimeConservedEssentialSpliceSite = true;
//                            }
//
//                            if (dnc2=="GT")
//                            {
//                                threePrimeConservedEssentialSpliceSite = true;
//                            }
//                        }
//
//                        if(strand=='-')
//                        {
//                            ref_seq->getString(dnc2, chromNo, (genomeIndex_t) start-2, 2);
//                            ref_seq->getString(dnc1, chromNo, (genomeIndex_t) end+1, 2);
//
//                            if (dnc1=="CT")
//                            {
//                                fivePrimeConservedEssentialSpliceSite = true;
//                            }
//
//                            if (dnc2=="AC")
//                            {
//                                threePrimeConservedEssentialSpliceSite = true;
//                            }
//                        }
//                  }
//
//                    if (feature=="stop_codon")
//                    {
//                        if (exists(GENCODE, chrom))
//                        {
//                            GENCODE[chrom]->search(start, end, intervals);
//                        }
//
//                        for (uint32_t i=0; i<intervals.size(); ++i)
//                        {
//                            GTFRecord* record = (GTFRecord*)intervals[i];
//                            if (record->feature == "exon" && record->gene == gene)
//                            {
//                                record->containsStopCodon = true;
//                            }
//                        }
//                    }
//
//                    if (feature=="start_codon")
//                    {
//                        if (exists(GENCODE, chrom))
//                        {
//                            GENCODE[chrom]->search(start, end, intervals);
//                        }
//
//                        for (uint32_t i=0; i<intervals.size(); ++i)
//                        {
//                            GTFRecord* record = (GTFRecord*)intervals[i];
//                            if (record->feature == "exon" && record->gene == gene)
//                            {
//                                record->containsStartCodon = true;
//                            }
//                        }
//                    }
//
//                    GTFRecord* gtfRecord = new GTFRecord(chrom, start, end, strand,
//                                                         gene, feature, frame, exonNo,
//                                                         fivePrimeConservedEssentialSpliceSite, threePrimeConservedEssentialSpliceSite,
//                                                         containsStartCodon, containsStopCodon,
//                                                         level, attrib);
//
//                    //std::cerr << "chromosome " << chrom << "\n";
//                    GENCODE[chrom]->insert(gtfRecord);
//        //
//                    ++noRecords;
//                }
//
//
//            }
//
//        }
//
//
//
//      std::clog << " ... completed\n";
//      if (GENCODE.size()==0)
//      {
//          std::cerr << "No reference GENCODE features!\n";
//            exit(1);
//      }
//
//      std::cerr << "start validation\n";
//        GENCODE["X"]->validate();
//        //G = GENCODE["chr20"];
//        //G->validate();
////        std::cerr << "end validation\n";
////        std::cerr << "height : " << GENCODE["X"]->height << "\n";
////        std::cerr << "size : " << GENCODE["X"]->size() << "\n";
////        std::cerr << "height : " << G->height << "\n";
////        std::cerr << "size : " << G->size() << "\n";
//    }

    /**
     * Populate STR_REF tree
     */
    void populate_str_ref_tree()
    {
//        std::string line;
//      std::vector<std::string> fields;
//      std::vector<std::string> alleles;
//
//        std::clog << "Populating STR_REF tree ... \n";
//
//        std::vector<Interval*> intervals;
//        uint32_t noRecords = 0;
//
//      //for gencode records
//        boost::regex regex("chr(.+)");
//      boost::smatch regexResult;

//        tbx_t *tbx;
//      hts_itr_t *iter;
//      BGZF* gffgz = xbgzf_open(str_ref_bed_file.c_str(), "r");
//        if (!(tbx=tbx_index_load(str_ref_bed_file.c_str())))
//        {
//          std::cerr << "Fail to load the VCF GZ index\n";
//          abort();
//      }
//
//          std::stringstream region;
//      kstring_t s;
//        s.s = 0; s.l = s.m = 0;
//        for (int32_t i = 0; i < nseqs; ++i)
//      {
//          //choose records from chromosome i
//          region.str("");
//          region << "chr" << seqnames[i] ;
//          if (!(iter = tbx_itr_querys(tbx, region.str().c_str())))
//          {
//              continue;
//          }
//
//          while (tbx_itr_next(gffgz, tbx, iter, &s) >= 0)
//            {
//                std::string line = std::string(s.s);
//
//              if(line.at(0)!='#')
//              {
//                  split(fields, '\t', line, 9);
//                    std::string chrom = fields[0];
//
//                  if (boost::regex_match(chrom, regexResult, regex))
//                  {
//                      chrom = regexResult[1];
//                  }
//
//                    //create tree for chromosome
//                    if(!exists(STR_REF, chrom))
//                    {
//                        std::clog << "Creating tree for chromosome " << chrom << "\n";
//                        STR_REF[chrom] = new IntervalTree();
//                    }
//
//                    //process fields
//                    std::map<std::string, std::string> attrib_map;
//
//
//                  uint32_t start = boost::lexical_cast<int32_t>(fields[1]);
//                  uint32_t end = boost::lexical_cast<int32_t>(fields[2]);
//
//                    BEDRecord* bedRecord = new BEDRecord(chrom, start, end);
//
//                    STR_REF[chrom]->insert(bedRecord);
//
//                    ++noRecords;
//                }
//            }
//        }

//      std::clog << " ... completed\n";
//      if (STR_REF.size()==0)
//      {
//          std::cerr << "No reference STR_REF features!\n";
//            exit(1);
//      }
//
//      std::cerr << "start validation\n";
//        STR_REF["X"]->validate();
//        //G = GENCODE["chr20"];
//        //G->validate();
//        std::cerr << "end validation\n";
//        std::cerr << "height : " << STR_REF["X"]->height << "\n";
//        std::cerr << "size : " << STR_REF["X"]->size() << "\n";
    }

    /**
     * Populate INDEL_HOTSPOTS tree
     */
//    void populate_indel_hotspots_tree()
//    {
//        std::string line;
//      std::vector<std::string> fields;
//      std::vector<std::string> alleles;
//
//        std::clog << "Populating INDEL_HOTSPOTS tree ... \n";
//
//        std::vector<Interval*> intervals;
//        uint32_t noRecords = 0;
//
//      //for gencode records
//        boost::regex regex("chr(.+)");
//      boost::smatch regexResult;
//
//        tbx_t *tbx;
//      hts_itr_t *iter;
//      BGZF* gffgz = xbgzf_open(indel_hotspot_bed_file.c_str(), "r");
//        if (!(tbx=tbx_index_load(indel_hotspot_bed_file.c_str())))
//        {
//          std::cerr << "Fail to load the VCF GZ index\n";
//          abort();
//      }
//
//          std::stringstream region;
//      kstring_t s;
//        s.s = 0; s.l = s.m = 0;
//        for (int32_t i = 0; i < nseqs; ++i)
//      {
//          //choose records from chromosome i
//          if (!(iter = tbx_itr_querys(tbx,  seqnames[i])))
//          {
//              continue;
//          }
//
//          while (tbx_itr_next(gffgz, tbx, iter, &s) >= 0)
//            {
//                std::string line = std::string(s.s);
//
//              if(line.at(0)!='#')
//              {
//                  split(fields, '\t', line, 9);
//                    std::string chrom = fields[0];
//
//                  if (boost::regex_match(chrom, regexResult, regex))
//                  {
//                      chrom = regexResult[1];
//                  }
//
//                    //create tree for chromosome
//                    if(!exists(INDEL_HOTSPOTS, chrom))
//                    {
//                        std::clog << "Creating tree for chromosome " << chrom << "\n";
//                        INDEL_HOTSPOTS[chrom] = new IntervalTree();
//                    }
//
//                    //process fields
//                    std::map<std::string, std::string> attrib_map;
//
//                  uint32_t start = boost::lexical_cast<int32_t>(fields[1]);
//                  uint32_t end = boost::lexical_cast<int32_t>(fields[2]);
//
//                    BEDRecord* bedRecord = new BEDRecord(chrom, start, end);
//
//                    INDEL_HOTSPOTS[chrom]->insert(bedRecord);
//
//                    ++noRecords;
//                }
//            }
//        }
//
//      std::clog << " ... completed\n";
//      if (INDEL_HOTSPOTS.size()==0)
//      {
//          std::cerr << "No reference INDEL_HOTSPOTS features!\n";
//            exit(1);
//      }
//
//      std::cerr << "start validation\n";
//        INDEL_HOTSPOTS["X"]->validate();
//        //G = GENCODE["chr20"];
//        //G->validate();
//        std::cerr << "end validation\n";
//        std::cerr << "height : " << INDEL_HOTSPOTS["X"]->height << "\n";
//        std::cerr << "size : " << INDEL_HOTSPOTS["X"]->size() << "\n";
//    }

    private:
//    bool exists(std::map<std::string, IntervalTree*>& map, const std::string& key)
//    {
//        return map.end()!=map.find(key);
//    }
//
//    /**
//     *Splits a line into a map - PERL style
//     */
//    void splitGTFAttributeFields(std::map<std::string, std::string>& map, std::string& str)
//    {
//      map.clear();
//      const char* tempStr = str.c_str();
//      int32_t i=0, lastIndex = str.size()-1;
//      std::string key;
//        std::string val;
//      std::stringstream token;
//
//      if (lastIndex<0) return;
//
//      while (i<=lastIndex)
//      {
//          //read next character
//          if(tempStr[i]!=';' && tempStr[i]!=' ')
//          {
//              token << tempStr[i];
//          }
//
//            //store key-value pair
//          if (i==lastIndex || (tempStr[i]==';' && tempStr[i+1]==' '))
//          {
//              val = token.str();
//              if (val.at(0)=='"')
//              {
//                  val = token.str();
//                  val = val.substr(1,val.size()-2);
//              }
//              map[key] = val;
//
//              token.str("");
//              key.clear();
//
//              if (i!=lastIndex)
//                  ++i;
//          }
//
//          //store key
//          if (tempStr[i]==' ')
//          {
//              key = token.str();
//              token.str("");
//          }
//
//          ++i;
//      }
//    };
};

//    bool exists(std::map<std::string, IntervalTree*>& map, const std::string& key)
//    {
//        return map.end()!=map.find(key);
//    }

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

}

void annotate_indels(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.annotate_indels();
    igor.print_stats();
};

//int annotate_indels1(int argc, char ** argv)
//{
//  //options
//  std::string inputVCFFileName;
//    std::string otxt_file;
//    std::string gencode_gtf_file;
//      std::string ref_fasta_file  ;
//  std::string chromosome;
//
//  try
//  {
//      std::string desc = "Example:\n\
//./vleftalign  -i pscalare.vcf -o pscalare.vcf -r hg37.fa \n\
//Left aligns variants in a VCF file.  Sorting and duplicate removal is required after left alignment. Multi-allelic variants are handled too\n";
//
//          std::string version = "0.5";
//      TCLAP::CmdLine cmd(desc, ' ', version);
//      TCLAP::ValueArg<std::string> arg_ivcf_file("i", "input_vcf", "Input Probe VCF file", true, "", "string");
//      TCLAP::ValueArg<std::string> arg_otxt_file("o", "output_txt", "Output text file", false, "-", "string");
//      TCLAP::ValueArg<std::string> arg_gencode_gtf_file("g", "gencode", "GENCODE GTF file", false, "/net/fantasia/home/atks/ref/encode/gencode.v15.annotation.gtf.gz", "string");
//      TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "ref_genome", "Reference Genome Fasta file", false, "/net/fantasia/home/atks/ref/genome/hs37d5.fa", "string");
//      TCLAP::ValueArg<std::string> argChrom("c", "chromosome", "Chromosome", false, "all", "string");
//
//      cmd.add(arg_ivcf_file);
//      cmd.add(arg_otxt_file);
//      cmd.add(arg_gencode_gtf_file);
//      cmd.add(arg_ref_fasta_file);
//      cmd.add(argChrom);
//      cmd.parse(argc, argv);
//
//      inputVCFFileName = arg_ivcf_file.getValue();
//      otxt_file = arg_otxt_file.getValue();
//      gencode_gtf_file = arg_gencode_gtf_file.getValue();
//      ref_fasta_file   = arg_ref_fasta_file.getValue();
//      chromosome = argChrom.getValue();
//
//      std::clog << "Profile Indels\n";
//      std::clog << "==============\n";
//      std::clog << "Options\n";
//      std::clog << "=======\n";
//      std::clog << "Input VCF file        : " << inputVCFFileName << "\n";
//      std::clog << "Output summary file   : " << otxt_file << "\n";
//      std::clog << "Gencode file          : " << gencode_gtf_file << "\n";
//      std::clog << "Reference Genome file : " << ref_fasta_file    << "\n";
//      std::clog << "Chromosome            : " << chromosome << "\n";
//      std::clog << "==============\n";
//  }
//  catch (TCLAP::ArgException &e)
//  {
//      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
//      abort();
//  }
//
//    //////////////////////
//    //i/o initialization//
//    //////////////////////
//
//  int32_t ftype = fileType(inputVCFFileName, 'r');
//  IFILE IN_VCF = ifopen(inputVCFFileName.c_str(), "r", (ftype==0?InputFile::UNCOMPRESSED:(ftype==1?InputFile::GZIP:InputFile::BGZF)));
//
//  ftype = fileType(otxt_file, 'w');
//  InputFile& OUT_VCF = *(ifopen(otxt_file.c_str(), "w", (ftype==0?InputFile::UNCOMPRESSED:(ftype==1?InputFile::GZIP:InputFile::BGZF))));
//
//
//
//    ftype = fileType(inputVCFFileName, 'r');
//  IFILE IN_GTF = ifopen(gencode_gtf_file.c_str(), "r", (ftype==0?InputFile::UNCOMPRESSED:(ftype==1?InputFile::GZIP:InputFile::BGZF)));
//
//    GenomeSequence* myRefSeq = new GenomeSequence(ref_fasta_file  );
//
//  std::string line;
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
//  bool allChromosomes = chromosome == "all";
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
////    std::cerr << "height : " << G->height << "\n";
////    std::cerr << "size : " << G->size() << "\n";
//
//    // GENCODE.print();
//    uint32_t CHROM = 0;
//  uint32_t POS = 1;
//  //uint32_t ID = 2;
//  uint32_t REF = 3;
//  uint32_t ALT = 4;
//  //uint32_t QUAL = 5;
//  //uint32_t FILTER = 6;
//  uint32_t INFO = 7;
//  uint32_t REST = 8;
//
////    exit(0);
////    std::cerr << "searching\n";
////    GENCODE.RBSearch(68000, 77000, intervals);
////
////    for (uint32_t i=0; i<intervals.size(); ++i)
////    {
////        std::cerr << ((GTFRecord*)intervals[i])->feature << "," << intervals[i]->start << "," << intervals[i]->end << "\n";
////    }
////    std::cerr << "search results printed\n";
//
//
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
//
////        char strand = exons[i]->strand;
////        uint32_t start = exons[i]->start;
////        uint32_t end = exons[i]->end;
////        uint32_t chromNo = myRefSeq->getChromosome(exons[i]->chrom.c_str());
////        std::string dnc1;
////        std::string dnc2;
////
////        if(strand=='+')
////        {
////            myRefSeq->getString(dnc1, chromNo, (genomeIndex_t) start-2, 2);
////            myRefSeq->getString(dnc2, chromNo, (genomeIndex_t) end+1, 2);
////
////            //if (dnc1=="AG" && dnc2=="GT")
////            //if (((exons[i]->containsStartCodon || dnc1=="AG") && (exons[i]->containsStopCodon || dnc2=="GT")) ||
////            //     (!exons[i]->containsStopCodon && !exons[i]->containsStartCodon && (dnc1=="AG" || dnc2=="GT")))
////            if (((exons[i]->exonNo==1 || dnc1=="AG") && (exons[i]->containsStopCodon || dnc2=="GT")))
////            {
////                ++noConservedSpliced;
////                //std::cerr << "+ conserved " << dnc1 << " " << dnc2 << "\n";
////            }
////            else
////            {
////                ++noUnconservedSpliced;
////                //std::cerr << "+ unconserved " << dnc1 << " " << dnc2 << "\n";
////                //exons[i]->print();
////            }
////        }
////
////        if(strand=='-')
////        {
////            myRefSeq->getString(dnc2, chromNo, (genomeIndex_t) start-2, 2);
////            myRefSeq->getString(dnc1, chromNo, (genomeIndex_t) end+1, 2);
////
////            //if (dnc1=="CT" && dnc2=="AC")
////            //if ((exons[i]->containsStartCodon || dnc1=="CT") &&
////            //    (exons[i]->containsStopCodon || dnc2=="AC"))
////            //if (((exons[i]->containsStartCodon || dnc1=="CT") && (exons[i]->containsStopCodon || dnc2=="AC")) ||
////            //     (!exons[i]->containsStopCodon && !exons[i]->containsStartCodon && (dnc1=="CT" || dnc2=="AC")))
////            if (((exons[i]->exonNo==1 || dnc1=="CT") && (exons[i]->containsStopCodon || dnc2=="AC")))
////            {
////                ++noConservedSpliced;
////                //std::cerr << "- conserved " << dnc1 << " " << dnc2 << "\n";
////            }
////            else
////            {
////                ++noUnconservedSpliced;
////                //std::cerr << "- unconserved " << dnc1 << " " << dnc2 << "\n";
////                //exons[i]->print();
////            }
////        }
//    }
//
//  uint32_t noVariants = 0;
//  uint32_t noFSIndels = 0;
//  uint32_t noNFSIndels = 0;
//  uint32_t noSpliceOverlaps = 0;
//    uint32_t noNonSpliceOverlaps = 0;
//
//  //std::vector<Interval*> intervals;
//    //std::vector<std::string> alleles;
//
//
//
//    //process VCF
//  while (!ifeof(IN_VCF))
//  {
//      IN_VCF >> line;
//
//      //std::cerr << line ;
//
//      if(line.at(0)=='#')
//      {
//          //OUT_VCF << line << "\n";
//      }
//      else
//      {
//          split(fields, '\t', line, 9);
//
//          std::string chrom = fields[CHROM];
//          int32_t pos = boost::lexical_cast<int32_t>(fields[POS]);
//
//          std::string ref = fields[REF];
//          std::string alts = fields[ALT];
//          std::map<std::string,std::string> info;
//          int32_t posEnd = pos + ref.size()-1;
//
//
//          splitInfoFields(info, fields[INFO]);
//
//            alleles.clear();
//          alleles.push_back(ref);
//
//
//          split(alleles, ',', alts, UINT_MAX, false);
//
////          std::cerr << "alleles " << alleles.size() << "\n";
////          std::cerr << "search results printed\n";
////            if ($F[$FEATURE] eq "CDS")
////            {
////                if ($start<=$F[$END]-1 && $end>=$F[$START])
////                {
////                    my @alts = split(",", $alt);
////                    my $refLength = length($ref);
////                    my %diff = ();
////                    $diff{($refLength-length($_))%3} = 1 foreach @alts;
////
////                    if (exists($diff{1})||exists($diff{2}))
////                    {
////                        $GENCODE{"CDS_FS"} = 1;
////                    }
////                    else
////                    {
////                        $GENCODE{"CDS_NFS"} = 1;
////                    }
////                }
////            }
//
//            if (alleles.size()==2)
//            {
//                std::string alt = alleles[1];
//                //std::cerr << "caught " << ref << " " << alt << "\n";
//
//                //SNP
//              if (ref.size()==1 && alt.size()==1)
//              {
//                 // std::cerr << "SNP\n";
//                  //OUT_VCF << line << "\n";
//              }
//              //MNP
//              else if (ref.size()==alt.size())
//              {
//                 // std::cerr << "MNP\n";
//                  //OUT_VCF << line << "\n";
//              }
//              //Indel + Cplxsub
//              else
//              {
//                  //std::cerr << "INDEL ";
//                  //std::cerr << "searching [" << (pos+1) << "," <<  posEnd << "]\n";
//                    //std::cerr << "RU " << info["RU"] << "\n";
//
//
//
//                    intervals.clear();
//                    if (exists(GENCODE, chrom))
//                    {
//                        GENCODE[chrom]->search(pos+1-2, posEnd+2, intervals);
//                  }
//
//                  //G = GENCODE[chrom];
//                  //G->search(pos+1, pos+ref.size()-1, intervals);
//
//                  //1-nfs
//                  //2-fs
//                  uint32_t frameShiftState = 0;
//                  //1-splice overlap
//                  //2-non overlap
//                  uint32_t spliceOverlapState = 0;
//
//                  for (uint32_t i=0; i<intervals.size(); ++i)
//                    {
//                        GTFRecord* record = (GTFRecord*)intervals[i];
//                      //std::cerr << ((GTFRecord*)intervals[i])->feature << "," << intervals[i]->start << "," << intervals[i]->end << "\n";
//
//                        if (record->feature == "CDS" &&
//                            record->end>=pos+1 && record->start<=posEnd)
//                        {
//                            if (abs((int32_t)ref.size()-(int32_t)alt.size())%3==0 || info["RU"].size()==3)
//                            {
//                                frameShiftState = 2;
//                            }
//                            else
//                            {
//                                frameShiftState = 1;
//                            }
//
//                            //std::cerr << "CDS: " << ref.size() << " " << alt.size() << ":" << ref << " " << alt << " " << state << "\n";
//                        }
//
//                        if (record->feature == "exon")
//                        {
//                            if ((1||(record->strand=='+' && record->exonNo!=1 && record->fivePrimeConservedEssentialSpliceSite) ||
//                                 (record->strand=='-' && !record->containsStopCodon && record->threePrimeConservedEssentialSpliceSite)) &&
//                                 record->start-1>=pos+1 && record->start-2<=posEnd)
//                            {
//                                spliceOverlapState = 1;
//                            }
//                            else if ((1||(record->strand=='+' && !record->containsStopCodon && record->threePrimeConservedEssentialSpliceSite) ||
//                                      (record->strand=='-' && record->exonNo!=1 && record->fivePrimeConservedEssentialSpliceSite)) &&
//                                       record->end+2>=pos+1 && record->end+1<=posEnd)
//                            {
//                                spliceOverlapState = 1;
//                            }
//                            else
//                            {
//                                spliceOverlapState = 2;
//                            }
//
//                            //std::cerr << "CDS: " << ref.size() << " " << alt.size() << ":" << ref << " " << alt << " " << state << "\n";
//                        }
//
//                    }
//
//                  if (frameShiftState ==1)
//                    {
//                        ++noFSIndels;
//
// //                     std::cerr << "CDS: " << ref.size() << " " << alt.size() << ":" << ref << " " << alt << " " << state  << " " << fields[INFO]<< "\n";
//
//                    }
//                    else if (frameShiftState ==2)
//                    {
//                        ++noNFSIndels;
//                    }
//
//                    if (spliceOverlapState ==1)
//                    {
//                        ++noSpliceOverlaps;
//
// //                     std::cerr << "CDS: " << ref.size() << " " << alt.size() << ":" << ref << " " << alt << " " << state  << " " << fields[INFO]<< "\n";
//
//                    }
//                    else if (spliceOverlapState ==2)
//                    {
//                        ++noNonSpliceOverlaps;
//                    }
//
////                    if (!changed)
////                    {
////                        OUT_VCF << line << "\n";
////                    }
////                    else
////                    {
////                        if (alleles.size()==2)
////                        {
////                            ++noLeftAligned;
////                        }
////                        else
////                        {
////                            ++noMultiLeftAligned;
////                        }
////                        OUT_VCF << fields[CHROM] << "\t"
////                                << pos << "\t"
////                                << fields[ID] << "\t"
////                                << alleles[alleles.size()-1] << "\t"
////                                << alleles[0];
////
////                        for (uint32_t i=1; i<(alleles.size()-1); ++i)
////                        {
////                            OUT_VCF << "," << alleles[i];
////                        }
////
////                        OUT_VCF << "\t" << fields[QUAL] << "\t"
////                                        << fields[FILTER] << "\t"
////                                        << fields[INFO] << ";"
////                                        << "OLD_VARIANT=" << fields[CHROM] << ":" << fields[POS] << ":" << fields[REF] << ":" << fields[ALT]
////                                        << (fields.size()>8 ? "\t" : "" )
////                                        << (fields.size()>8 ? fields[REST] : "") << "\n";
////                    }
//
//
//                  ++noVariants;
//              }
//
//
//            }
//            //MULTIALLELICS
//            else
//            {
//            }
//      }
//  }
//
//    ifclose(IN_VCF);
//    ifclose(&OUT_VCF);
//
//    for (std::map<std::string, IntervalTree*>::iterator i = GENCODE.begin();  i!=GENCODE.end(); ++i)
//    {
//        delete i->second;
//    }
//
//    std::clog << "No. biallelic variants : " << noVariants << "\n";
//    std::clog << "No. biallelic FS Indels : " << noFSIndels << "\n";
//    std::clog << "No. biallelic NFS Indels : " << noNFSIndels << "\n";
//    std::clog << "No. biallelic Splice Indels : " << noSpliceOverlaps << "\n";
//    std::clog << "No. biallelic Non Splice Indels : " << noNonSpliceOverlaps << "\n";
//
//    std::clog << "No. Conserved Spliced Sites   : " << noConservedSpliced << "\n";
//    std::clog << "No. Unconserved Spliced Sites : " << noUnconservedSpliced << "\n";
//
//
//    //generate overlap analysis
//    //type 1-3 matchings
//
//
//    return 0;
//};
