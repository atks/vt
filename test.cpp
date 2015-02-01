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

#include "align.h"

namespace
{
void print_time(float t)
{
    if (t<60)
    {
        fprintf(stderr, "Alignment time elapsed: %fs\n\n", t);
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
        fprintf(stderr, "Alignment time elapsed: %dd %dh %dm %.6fs\n\n", ((int32_t)(t/(60*60*24))), ((int32_t)(h/(60*60))), ((int32_t)(m/60)), (fmod(m, 60)));
    }
};

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string method;
    std::vector<uint32_t> x;


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
            TCLAP::UnlabeledMultiArg<uint32_t> arg_x("ap", "#ploidy #alleles", true, "", cmd);

            cmd.parse(argc, argv);

            x = arg_x.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }

        no = 1;
    };

    /**
     * n choose r.
     */
    uint32_t choose(uint32_t n, uint32_t r)
    {
        if (r>n)
        {
            return 0;
        }
        else if (r==n)
        {
            return 1;
        }
        else if (r==0)
        {
            return 1;
        }
        else
        {
            if (r>(n>>1))
            {
                r = n-r;
            }

            uint32_t num = n;
            uint32_t denum = 1;

            for (uint32_t i=1; i<r; ++i)
            {
                num *= n-i;
                denum *= i+1;
            }

            return num/denum;
        }
    }

    /**
     * Gets number of genotypes from number of alleles and ploidy.
     */
    uint32_t bcf_ap2g(uint32_t no_allele, uint32_t no_ploidy)
    {
        if (no_ploidy<=no_allele)
        {
            uint32_t no_genotypes = 0;

            for (uint32_t i=1; i<=no_ploidy; ++i)
            {
                uint32_t n = 0;

                if (i<no_ploidy)
                {
                    n += choose(no_allele, i)*choose(no_ploidy-1, i-1);
                }
                else if (i==no_ploidy)
                {
                    n += choose(no_allele, no_ploidy);
                }

                no_genotypes += n;
            }

            return no_genotypes;
        }
        else // alleles less than ploidy
        {
            return choose(no_ploidy+no_allele-1, no_allele-1);
        }

        return 0;
    }

    /**
     * Gets number of genotypes from number of alleles and ploidy.
     */
    uint32_t bcf_g2i(std::string genotype)
    {
        uint32_t allele = genotype.at(genotype.size()-1)-65;
        
        if (genotype.size()==1)
        {
            return allele;
        }    
        else
        {
            return bcf_ap2g(allele, genotype.size()) +  bcf_g2i(genotype.substr(0,genotype.size()-1));
        }
    }

    void print_genotypes(uint32_t A, uint32_t P, std::string genotype)
    {
        if (genotype.size()==P)
        {
            std::cerr << no << ") " << genotype << " " << bcf_g2i(genotype) << "\n";
            ++no;
        }
        else
        {
            for (uint32_t a=0; a<A; ++a)
            {
                std::string s(1,(char)(a+65));
                s.append(genotype);
                print_genotypes(a+1, P, s);
            }
        }
    }

    void test()
    {
        print_genotypes(x[0], x[1], "");
        uint32_t g = bcf_ap2g(x[0], x[1]);

        std::cerr << "A: " << x[0] << " P: " << x[1] << " G: " << g << "\n";
    };

    void print_stats()
    {
        std::clog << "\n";

    };

    ~Igor() {};

    private:
};

}

void test(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.test();

};
