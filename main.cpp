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

#include "normalize.h"
#include "merge_duplicate_variants.h"
#include "peek.h"
#include "construct_probes.h"
#include "discover.h"
#include "annotate_variants.h"
#include "genotype.h"
#include "merge_candidate_variants.h"
#include "merge.h"
#include "concat.h"
#include "partition.h"
#include "view.h"
#include "index.h"
#include "profile_indels.h"
#include "profile_mendel_errors.h"

void print_time(double t)
{
    if (t<60)
    {
        fprintf(stderr, "Time elapsed: %.2fs\n\n", t);
    }
    else if (t<60*60) //less than an hour
    {
        fprintf(stderr, "Time elapsed: %dm %ds\n\n", ((int32_t)(t/60)), ((int32_t)fmod(t, 60)));
    }
    else if (t<60*60*24) //less than a day
    {
        double m = fmod(t, 60*60); //remaining minutes
        fprintf(stderr, "Time elapsed: %dh %dm %ds\n\n", ((int32_t)(t/(60*60))), ((int32_t)(m/60)), ((int32_t)fmod(m, 60)));
    }
    else if (t<60*60*24*365) //less than a year
    {
        double h = fmod(t, 60*60*24); //remaining hours
        double m = fmod(h, 60*60); //remaining minutes
        fprintf(stderr, "Time elapsed: %dd %dh %dm %ds\n\n", ((int32_t)(t/(60*60*24))), ((int32_t)(h/(60*60))), ((int32_t)(m/60)), ((int32_t)fmod(m, 60)));
    }
};

void help()
{
    std::clog << "Help page on http://statgen.sph.umich.edu/wiki/Vt\n";
    std::clog << "\n";
    std::clog << "view                      view vcf/vcf.gz/bcf files\n";
    std::clog << "index                     index vcf.gz/bcf files\n";
    std::clog << "normalize                 normalize variants\n";
    std::clog << "mergedups                 merge duplicate variants\n";
    std::clog << "merge                     merge VCF files\n";
    std::clog << "concat                    concatenate VCF files\n";
    std::clog << "annotate_variants         annotate variants\n";
    std::clog << "peek                      summary of variants in the vcf file\n";
    std::clog << "partition                 partition variants\n";
    std::clog << "profile_indels            profile indels\n";
    std::clog << "profile_mendel_errors     profile indels\n";
    std::clog << "discover                  discover variants\n";
    std::clog << "merge_candidate_variants  merge candidate variants\n";
    std::clog << "construct_probes          construct probes for each variant\n";
    std::clog << "genotype                  genotype variants\n";
    std::clog << "\n";
}

int main(int argc, char ** argv)
{
    clock_t t0;
    t0 = clock();
    bool print = true;
        
    if (argc==1)
    {
        help();
        exit(0);
    }
    
    std::string cmd(argv[1]);
    
    //primitive programs that do not require help pages and summary statistics by default
    if (argc>1 && cmd=="view")
    {
        print = view(argc-1, ++argv);
    } 
    else if (argc>1 && cmd=="index")
    {
        print = index(argc-1, ++argv);
    } 
    else if (argc>1 && cmd=="merge")
    {
        print = merge(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="concat")
    {
        print = concat(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="normalize")
    {
        normalize(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="mergedups")
    {
        merge_duplicate_variants(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="peek")
    {
        peek(argc-1, ++argv);
    }  
    else if (argc>1 && cmd=="partition")
    {
        partition(argc-1, ++argv);
    } 
    else if (argc>1 && cmd=="annotate_variants")
    {
        annotate_variants(argc-1, ++argv);
    }     
    else if (argc>1 && cmd=="discover")
    {
        discover(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="merge_candidate_variants")
    {
        merge_candidate_variants(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="genotype")
    {
        genotype(argc-1, ++argv);
    } 
    else if (argc>1 && cmd=="construct_probes")
    {
        construct_probes(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_indels")
    {
        profile_indels(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_mendel_errors")
    {
        profile_mendel_errors(argc-1, ++argv);
    }    
    else
    {
        std::clog << "Command not found: " << argv[1] << "\n\n";
        help();
        exit(1);
    }

    if (print)
    {
        clock_t t1;
        t1 = clock();
        print_time((float)(t1-t0)/CLOCKS_PER_SEC);
    }

    return 0;
}
