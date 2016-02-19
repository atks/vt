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

#include "compute_rl_dist.h"

namespace
{
class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string ref_fasta_file;
    std::vector<GenomeInterval> intervals;
    ReferenceSequence *rs;
    std::map<std::string, std::vector<int32_t> > stats;

    bool debug;
    uint32_t no;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "vt test  -m detect_motif -s ACTGACT \n";

            std::string version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::ValueArg<std::string> arg_chromosome("c", "c", "chromosome fasta file []", false, "", "str", cmd);

            cmd.parse(argc, argv);
            
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            ref_fasta_file = arg_ref_fasta_file.getValue();
            debug = arg_debug.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
        //////////////////////
        //i/o initialization//
        //////////////////////


        ////////////////////////
        //tools initialization//
        ////////////////////////
        rs = new ReferenceSequence(ref_fasta_file);
    }


    /**
     * Return the canonical representation of a motif.
     */
    std::string canonicalize(std::string& motif)
    {
        std::string cmotif = motif;

        for (uint32_t i=0; i<motif.size(); ++i)
        {
            std::string shifted_motif = motif.substr(i) + motif.substr(0,i);

            if (shifted_motif < cmotif)
            {
                cmotif = shifted_motif;
            }
        }

        return cmotif;
    }

    /**
     * Updates joint allele distribution.
     */
    void update_stats(std::string& del, int32_t rl)
    {
        std::string motif = canonicalize(del);
        
        if (stats.find(motif)==stats.end())
        {
            std::vector<int32_t> v;
            stats[motif] = v;
            stats[motif].resize(rl, 0);

            ++stats[motif][rl-1];
        }
        else
        {
            std::vector<int32_t>& v = stats[motif];

            if (rl-1>v.size()-1)
            {
                int32_t to_add = (rl - 1) - (v.size()-1);
                for (int32_t i=0; i<to_add; ++i)
                {
                    stats[motif].push_back(0);
                }
            }

            ++v[rl-1];
        }
    }

   /**
     * Print stats.
     */
    void print_dist()
    {
        std::map<std::string, std::vector<int32_t> >::iterator i = stats.begin();

        while (i!=stats.end())
        {
            std::string motif = i->first;

            std::vector<int32_t>& v = i->second;
            for (int32_t j=0; j<v.size(); ++j)
            {
                if (v[j]!=0)
                {
                    std::cout << motif << "\t" << (j+1) << "\t" << v[j] << "\n";
                }
            }
            
            ++i;
        }
    }


   /**
     * Print stats.
     */
    void print_dist_err()
    {
        std::map<std::string, std::vector<int32_t> >::iterator i = stats.begin();

        while (i!=stats.end())
        {
            std::string motif = i->first;

            std::vector<int32_t>& v = i->second;
            for (int32_t j=0; j<v.size(); ++j)
            {
                if (v[j]!=0)
                {
                    std::cerr << motif << "\t" << (j+1) << "\t" << v[j] << "\n";
                }
            }
            
            ++i;
        }
    }
    
    /**
     * Gets number of genotypes from number of alleles and ploidy.
     */
    void compute_rl_dist()
    {
        int32_t L = 3;

        if (intervals.size()==0)
        {
            parse_intervals(intervals, "", "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y");
            print_dist_err();
        }    


        for (uint32_t i = 0; i<intervals.size(); ++i)
        {
            std::string chrom = intervals[i].seq;
            if (intervals[i].end1 == ((1<<29) - 1))
            {
                intervals[i].end1 = rs->fetch_seq_len(chrom);
            }    
            
            int32_t beg1 = intervals[i].start1;
            int32_t end1 = intervals[i].end1;
            std::string del;

            std::cerr << "#" << chrom << ":" << beg1 << "-" << end1 << "\n";
//            print_dist_err();
        
            //for each position
            for (int32_t pos1=beg1; pos1<=end1; ++pos1)
            {
                char b = rs->fetch_base(chrom, pos1);
                    
                    

                if (b=='N') continue;

                if (debug) std::cerr << chrom << ":" << pos1 << " " << b << "\n";
                if (pos1%1000000==0)
                {
                    std::cerr << "#" << chrom << ":" << pos1 << " " << b << "\n";
//                    print_dist_err();
                }
                for (int32_t offset = 0; offset<L; ++offset)
                {
                    rs->fetch_seq(chrom.c_str(), pos1, pos1+offset, del);
                    
                    if (del.find('N')!=std::string::npos || del.size()!=L)
                    {
                        continue;
                    }

                    //left alignment
                    int32_t ref_pos1 = pos1 + offset;
                    int32_t alt_pos1 = pos1 - 1;
                    char ref = del.at(offset);
                    char alt = rs->fetch_base(chrom, alt_pos1);
                    while (ref==alt)
                    {
                        --ref_pos1;
                        --alt_pos1;

                        ref = rs->fetch_base(chrom, ref_pos1);
                        alt = rs->fetch_base(chrom, alt_pos1);
                    }
                    int32_t lflank_pos1 = alt_pos1;
                    
                    //right alignment
                    ref_pos1 = pos1;
                    alt_pos1 = pos1 + offset + 1;
                    ref = del.at(0);
                    alt = rs->fetch_base(chrom, alt_pos1);
                    while (ref==alt)
                    {
                        ++ref_pos1;
                        ++alt_pos1;

                        ref = rs->fetch_base(chrom, ref_pos1);
                        alt = rs->fetch_base(chrom, alt_pos1);
                    }
                    int32_t rflank_pos1 = alt_pos1;

                    int32_t rl = rflank_pos1 - lflank_pos1 - 1;

                    update_stats(del, rl);

                    if (debug) std::cerr << "\t" << del << " [" << lflank_pos1 << "," << rflank_pos1 << "] (" << rl << ")\n";
                }
            }
        }
        
        print_dist();
    }

    void print_options()
    {
//        if (!print) return;
//
//        std::clog << "normalize v" << version << "\n";
//        std::clog << "\n";
//        std::clog << "options:     input VCF file                                  " << input_vcf_file << "\n";
//        std::clog << "         [o] output VCF file                                 " << output_vcf_file << "\n";
//        std::clog << "         [w] sorting window size                             " << window_size << "\n";
//        print_str_op("         [f] filter                                          ", fexp);
//        std::clog << "         [n] no fail on reference inconsistency for non SNPs " << (!strict ? "true" : "false") << "\n";
//        std::clog << "         [q] quiet                                           " << (!print ? "true" : "false") << "\n";
//        std::clog << "         [d] debug                                           " << (debug ? "true" : "false")  << "\n";
//        std::clog << "         [r] reference FASTA file                            " << ref_fasta_file << "\n";
//        print_int_op("         [i] intervals                                       ", intervals);
//        std::clog << "\n";
    }

    void print_stats()
    {
//        if (!print) return;
//
//        int32_t no_biallelic_normalized = no_lt+no_rt+no_lt_rt+no_rt_la+no_la;
//        int32_t no_multiallelic_normalized = no_multi_lt+no_multi_rt+no_multi_lt_rt+no_multi_rt_la+no_multi_la;
//        int32_t no_normalized = no_biallelic_normalized + no_multiallelic_normalized;
//
//        std::clog << "\n";
//        std::clog << "stats: biallelic\n";
//        std::clog << "          no. left trimmed                      : " << no_lt << "\n";
//        std::clog << "          no. right trimmed                     : " << no_rt << "\n";
//        std::clog << "          no. left and right trimmed            : " << no_lt_rt << "\n";
//        std::clog << "          no. right trimmed and left aligned    : " << no_rt_la << "\n";
//        std::clog << "          no. left aligned                      : " << no_la << "\n";
//        std::clog << "\n";
//        std::clog << "       total no. biallelic normalized           : " << no_biallelic_normalized << "\n";
//        std::clog << "\n";
//        std::clog << "       multiallelic\n";
//        std::clog << "          no. left trimmed                      : " << no_multi_lt << "\n";
//        std::clog << "          no. right trimmed                     : " << no_multi_rt << "\n";
//        std::clog << "          no. left and right trimmed            : " << no_multi_lt_rt << "\n";
//        std::clog << "          no. right trimmed and left aligned    : " << no_multi_rt_la << "\n";
//        std::clog << "          no. left aligned                      : " << no_multi_la << "\n";
//        std::clog << "\n";
//        std::clog << "       total no. multiallelic normalized        : " << no_multiallelic_normalized << "\n";
//        std::clog << "\n";
//        std::clog << "       total no. variants normalized            : " << no_normalized << "\n";
//        std::clog << "       total no. variants observed              : " << no_variants << "\n";
//        std::clog << "       total no. reference observed             : " << no_refs << "\n";
//        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void compute_rl_dist(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.compute_rl_dist();
    igor.print_stats();
};
