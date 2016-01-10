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

#include "align.h"
#include "annotate_1000g.h"
#include "annotate_dbsnp_rsid.h"
#include "annotate_indels.h"
#include "annotate_regions.h"
#include "annotate_variants.h"
#include "cat.h"
#include "compute_features.h"
#include "compute_concordance.h"
#include "config.h"
#include "consolidate.h"
#include "consolidate_vntrs.h"
#include "construct_probes.h"
#include "decompose_blocksub.h"
#include "decompose.h"
#include "decompose2.h"
#include "discover.h"
#include "estimate.h"
#include "genotype2.h"
#include "genotype.h"
#include "hfilter.h"
#include "index.h"
#include "merge_candidate_variants.h"
#include "merge_candidate_variants2.h"
#include "merge.h"
#include "multi_partition.h"
#include "normalize.h"
#include "partition.h"
#include "paste.h"
#include "paste_and_compute_features_sequential.h"
#include "peek.h"
#include "profile_afs.h"
#include "profile_chm1.h"
#include "profile_chrom.h"
#include "profile_fic_hwe.h"
#include "profile_hwe.h"
#include "profile_indels.h"
#include "profile_len.h"
#include "profile_mendelian.h"
#include "profile_na12878.h"
#include "profile_snps.h"
#include "profile_vntrs.h"
#include "remove_overlap.h"
#include "rminfo.h"
#include "rpartition.h"
#include "seq.h"
#include "sort.h"
#include "subset.h"
#include "svm_predict.h"
#include "svm_train.h"
#include "test.h"
#include "union_variants.h"
#include "uniq.h"
#include "validate.h"
#include "view.h"
#include "vntrize.h"

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
    std::clog << "decompose                 decompose variants\n";
    std::clog << "uniq                      drop duplicate variants\n";
    std::clog << "cat                       concatenate VCF files\n";
    std::clog << "paste                     paste VCF files\n";
    std::clog << "sort                      sort VCF files\n";
    std::clog << "subset                    subset VCF file to variants polymorphic in a sample\n";
    std::clog << "\n";
    std::clog << "peek                      summary of variants in the vcf file\n";
    std::clog << "partition                 partition variants\n";
    std::clog << "multi_partition           partition variants from multiple VCF files\n";
    std::clog << "annotate_variants         annotate variants\n";
    std::clog << "annotate_db_rsid          annotate variants with dbSNP rsid\n";
    std::clog << "annotate_1000g            annotate variants with 1000 Genomes variants\n";
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
    else if (argc>1 && cmd=="paste_and_compute_features_sequential")
    {
        paste_and_compute_features_sequential(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="cat")
    {
        print = cat(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="subset")
    {
        subset(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="decompose")
    {
        decompose(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="decompose2")
    {
        decompose2(argc-1, ++argv);
    }    
    else if (argc>1 && cmd=="decompose_blocksub")
    {
        decompose_blocksub(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="normalize")
    {
        print = normalize(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="vntrize")
    {
        print = vntrize(argc-1, ++argv);
    }    
    else if (argc>1 && cmd=="validate")
    {
        print = validate(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="rminfo")
    {
        print = rminfo(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="sort")
    {
        print = sort(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="config")
    {
        config(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="uniq")
    {
        uniq(argc-1, ++argv);
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
    else if (argc>1 && cmd=="rpartition")
    {
        rpartition(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="svm_train")
    {
        svm_train(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="svm_predict")
    {
        svm_predict(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="multi_partition")
    {
        multi_partition(argc-1, ++argv);
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
    else if (argc>1 && cmd=="merge_candidate_variants2")
    {
        merge_candidate_variants2(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="union_variants")
    {
        union_variants(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="genotype")
    {
        genotype(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="genotype2")
    {
        genotype2(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="construct_probes")
    {
        construct_probes(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_vntrs")
    {
        profile_vntrs(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_indels")
    {
        profile_indels(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_snps")
    {
        profile_snps(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="seq")
    {
        print = seq(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_mendelian")
    {
        profile_mendelian(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_na12878")
    {
        profile_na12878(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="profile_chm1")
    {
        profile_chm1(argc-1, ++argv);
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
    else if (argc>1 && cmd=="compute_concordance")
    {
        compute_concordance(argc-1, ++argv);
    }    
    else if (argc>1 && cmd=="estimate")
    {
        estimate(argc-1, ++argv);
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
    else if (argc>1 && cmd=="annotate_indels")
    {
        annotate_indels(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="annotate_1000g")
    {
        annotate_1000g(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="consolidate")
    {
        consolidate(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="consolidate_vntrs")
    {
        consolidate_vntrs(argc-1, ++argv);
    }    
    else if (argc>1 && cmd=="filter")
    {
        hfilter(argc-1, ++argv);
    }
    else if (argc>1 && cmd=="test")
    {
        test(argc-1, ++argv);
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
