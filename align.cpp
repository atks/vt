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
    std::string x;
    std::string y;
    std::string lflank;
    std::string ru;
    std::string rflank;
    float delta;
    float epsilon;
    float tau;
    float eta;
    float mismatch_penalty;

    bool debug;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "vt align  -m nw -x ACTGACT -y CCCCCCACTGACTGGGGGG\n"
                     "          vt align  -m sw -x ACTGACT -y CCCCCCACTGACTGGGGGG\n"
                     "          vt align  -m ghmm -x ACTGACT -y CCCCCCACTGACTGGGGGG\n"
                     "\n";

            std::string version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_method("m", "method", "alignment Method", true, "", "string");
            TCLAP::ValueArg<std::string> arg_x("x", "x_seq", "X Sequence", false, "", "string");
            TCLAP::ValueArg<std::string> arg_y("y", "y_seq", "Y Sequence", true, "", "string");
            TCLAP::ValueArg<std::string> arg_lflank("l", "l", "left flanks", false, "", "string");
            TCLAP::ValueArg<std::string> arg_ru("u", "u", "repeat unit", false, "", "string");
            TCLAP::ValueArg<std::string> arg_rflank("r", "r", "right flanks", false, "", "string");
            TCLAP::ValueArg<float> arg_delta("d", "d", "delta", false, 0.01, "float");
            TCLAP::ValueArg<float> arg_epsilon("e", "e", "epsilon", false, 0.05, "float");
            TCLAP::ValueArg<float> arg_tau("t", "t", "tau", false, 0.01, "float");
            TCLAP::ValueArg<float> arg_eta("n", "n", "eta", false, 0.01, "float");
            TCLAP::ValueArg<float> arg_mismatch_penalty("p", "p", "mismatch penalty", false, 1, "float");
            TCLAP::SwitchArg arg_debug("v", "v", "debug mode", false);

            cmd.add(arg_method);
            cmd.add(arg_x);
            cmd.add(arg_y);
            cmd.add(arg_lflank);
            cmd.add(arg_ru);
            cmd.add(arg_rflank);
            cmd.add(arg_epsilon);
            cmd.add(arg_delta);
            cmd.add(arg_tau);
            cmd.add(arg_eta);
            cmd.add(arg_mismatch_penalty);
            cmd.add(arg_debug);
            cmd.parse(argc, argv);

            method = arg_method.getValue();
            x = arg_x.getValue();
            y = arg_y.getValue();
            lflank = arg_lflank.getValue();
            ru = arg_ru.getValue();
            rflank = arg_rflank.getValue();
            delta = arg_delta.getValue();
            epsilon = arg_epsilon.getValue();
            tau = arg_tau.getValue();
            eta = arg_eta.getValue();
            mismatch_penalty = arg_mismatch_penalty.getValue();
            debug = arg_debug.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void align()
    {
        if (method=="lhmm")
        {
            LHMM lhmm;
            double llk;
            std::string qual;
            for (int32_t i=0; i<y.size(); ++i)
            {
                qual += 'K';
            }
            clock_t t0 = clock();
            lhmm.align(llk, x.c_str(), y.c_str(), qual.c_str());
            lhmm.print_alignment();
            clock_t t1 = clock();
            print_time((float)(t1-t0)/CLOCKS_PER_SEC);
        }
        else if (method=="chmm")
        {
            CHMM hmm(debug);
            double llk;
            std::string qual;
            for (int32_t i=0; i<y.size(); ++i)
            {
                qual += 'K';
            }
            hmm.set_model(lflank.c_str(), ru.c_str(), rflank.c_str());
            hmm.set_delta(delta);
            hmm.set_epsilon(epsilon);
            hmm.set_tau(tau);
            hmm.set_eta(eta);
            hmm.set_mismatch_penalty(mismatch_penalty);
            hmm.initialize_T();
            clock_t t0 = clock();
            hmm.align(y.c_str(), qual.c_str());
            hmm.print_alignment();
            clock_t t1 = clock();
            std::cerr << "\n";
            print_time((float)(t1-t0)/CLOCKS_PER_SEC);
        }
        else if (method=="lfhmm")
        {
            LFHMM hmm(debug);
            double llk;
            std::string qual;
            for (int32_t i=0; i<y.size(); ++i)
            {
                qual += 'K';
            }
            clock_t t0 = clock();
            hmm.set_model(lflank.c_str(), ru.c_str());
            hmm.set_delta(delta);
            hmm.set_epsilon(epsilon);
            hmm.set_tau(tau);
            hmm.set_eta(eta);
            hmm.initialize_T();
            hmm.set_mismatch_penalty(mismatch_penalty);
            hmm.align(y.c_str(), qual.c_str());
            hmm.print_alignment();
            clock_t t1 = clock();
            print_time((float)(t1-t0)/CLOCKS_PER_SEC);
        }
        else if (method=="rfhmm")
        {
            RFHMM hmm(debug);
            std::string qual;
            for (int32_t i=0; i<y.size(); ++i)
            {
                qual += 'K';
            }
            clock_t t0 = clock();
            hmm.set_model(ru.c_str(), rflank.c_str());
            hmm.set_delta(delta);
            hmm.set_epsilon(epsilon);
            hmm.set_tau(tau);
            hmm.set_eta(eta);
            hmm.set_mismatch_penalty(mismatch_penalty);
            hmm.initialize_T();
            hmm.align(y.c_str(), qual.c_str());
            hmm.print_alignment();
            clock_t t1 = clock();
            print_time((float)(t1-t0)/CLOCKS_PER_SEC);
        }
        else if (method=="ahmm")
        {
            AHMM hmm(debug);
            std::string qual;
            for (int32_t i=0; i<y.size(); ++i)
            {
                qual += 'K';
            }
            clock_t t0 = clock();
            hmm.set_model(ru.c_str());
            hmm.set_delta(delta);
            hmm.set_epsilon(epsilon);
            hmm.set_tau(tau);
            hmm.set_eta(eta);
            hmm.set_mismatch_penalty(mismatch_penalty);
            hmm.initialize_T();
            hmm.align(y.c_str(), qual.c_str());
            hmm.print_alignment();
            clock_t t1 = clock();
            print_time((float)(t1-t0)/CLOCKS_PER_SEC);
        }
        else if (method=="lhmm1")
        {
            LHMM1 lhmm1;
            double llk;
            std::string qual;
            for (int32_t i=0; i<y.size(); ++i)
            {
                qual += 'K';
            }
            clock_t t0 = clock();
            lhmm1.align(llk, x.c_str(), y.c_str(), qual.c_str());
            lhmm1.printAlignment();
            clock_t t1 = clock();
            print_time((float)(t1-t0)/CLOCKS_PER_SEC);
        }
    };

    void print_options()
    {
        std::clog << "align v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     method      " << method << "\n";
        std::clog << "         [x] x        " << x << "\n";
        std::clog << "         [l] lflank   " << lflank << "\n";
        std::clog << "         [u] repeat   " << ru << "\n";
        std::clog << "         [r] rflank   " << rflank << "\n";
        std::clog << "         [y] y        " << y << "\n";
        std::clog << "         [d] delta    " << delta << "\n";
        std::clog << "         [e] epsilon  " << epsilon << "\n";
        std::clog << "         [t] tau      " << tau << "\n";
        std::clog << "         [n] eta      " << eta << "\n";
        std::clog << "         [p] p        " << mismatch_penalty << "\n";
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";

    };

    ~Igor() {};

    private:
};

}

void align(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.align();

};
