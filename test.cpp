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

            TCLAP::SwitchArg arg_debug("d", "debug", "Debug", false);

            cmd.add(arg_method);
            cmd.add(arg_x);
            cmd.add(arg_y);
            cmd.add(arg_lflank);
            cmd.add(arg_ru);
            cmd.add(arg_rflank);
            cmd.add(arg_debug);
            cmd.parse(argc, argv);

            method = arg_method.getValue();
            x = arg_x.getValue();
            y = arg_y.getValue();
            lflank = arg_lflank.getValue();
            ru = arg_ru.getValue();
            rflank = arg_rflank.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void test()
    {
        
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
