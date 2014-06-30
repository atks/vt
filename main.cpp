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
#include "annotate_str.h"
#include "annotate_dbsnp_rsid.h"
#include "genotype.h"
#include "genotype2.h"
#include "merge_candidate_variants.h"
#include "merge.h"
#include "concat.h"
#include "paste.h"
#include "subset.h"
#include "partition.h"
#include "view.h"
#include "index.h"
#include "profile_indels.h"
#include "profile_mendelian.h"
#include "decompose.h"
#include "remove_overlap.h"
#include "profile_na12878.h"
#include "profile_snps.h"
#include "align.h"
#include "compute_features.h"
#include "profile_afs.h"
#include "profile_hwe.h"
#include "profile_len.h"
#include "profile_chrom.h"
#include "annotate_regions.h"
#include "annotate_str.h"
#include "consolidate_variants.h"
#include "test.h"
#include "config.h"

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
    std::clog << "Useful tools:\n";
    std::clog << "view                      view vcf/vcf.gz/bcf files\n";
    std::clog << "index                     index vcf.gz/bcf files\n";
    std::clog << "normalize                 normalize variants\n";
    std::clog << "mergedups                 merge duplicate variants\n";
    std::clog << "concat                    concatenate VCF files\n";
    std::clog << "paste                     paste VCF files\n";
    std::clog << "subset                    subset VCF file to variants polymorphic in a sample\n";
    std::clog << "\n";
    std::clog << "peek                      summary of variants in the vcf file\n";
    std::clog << "partition                 partition variants\n";
    std::clog << "annotate_variants         annotate variants\n";
    std::clog << "annotate_regions          annotate regions\n";
    std::clog << "compute_concordance       compute genotype concordance between 2 call sets\n";
    std::clog << "compute_features          compute genotype likelihood based statistics\n";
    std::clog << "\n";
    std::clog << "discover                  discover variants\n";
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
    else if (argc>1 && cmd=="paste")
    {
        print = paste(argc-1, ++argv);
    }    
    else if (argc>1 && cmd=="concat")
    {
        print = concat(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="subset")
    {
        subset(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="decompose")
    {
        decompose(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="normalize")
    {
        print = normalize(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="config")
    {
        config(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="mergedups")
    {
        merge_duplicate_variants(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="remove_overlap")
    {
        remove_overlap(argc-1, ++argv);
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
    else if (argc>1 && cmd=="annotate_regions")
    {
        annotate_regions(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="annotate_dbsnp_rsid")
    {
        annotate_dbsnp_rsid(argc-1, ++argv);
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
        genotype2(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="characterize")
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
    else if (argc>1 && cmd=="profile_snps")
    {
        profile_snps(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_mendelian")
    {
        profile_mendelian(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_na12878")
    {
        profile_na12878(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_chrom")
    {
        profile_chrom(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="align")
    {
        align(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="compute_features")
    {
        compute_features(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_afs")
    {
        profile_afs(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_hwe")
    {
        profile_hwe(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_len")
    {
        profile_len(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="annotate_str")
    {
        annotate_str(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="consolidate_variants")
    {
        consolidate_variants(argc-1, ++argv);
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
