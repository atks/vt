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

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string method;
    std::string x;
    std::string y;
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
            TCLAP::ValueArg<std::string> arg_x("x", "x_seq", "X Sequence", true, "", "string");
            TCLAP::ValueArg<std::string> arg_y("y", "y_seq", "Y Sequence", true, "", "string");
            TCLAP::SwitchArg arg_debug("d", "debug", "Debug", false);

            cmd.add(arg_method);
            cmd.add(arg_x);
            cmd.add(arg_y);
            cmd.add(arg_debug);
            cmd.parse(argc, argv);

            method = arg_method.getValue();
            x = arg_x.getValue();
            y = arg_y.getValue();
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
            lhmm.align(llk, x.c_str(), y.c_str(), qual.c_str());
            lhmm.printAlignment();
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
            lhmm1.align(llk, x.c_str(), y.c_str(), qual.c_str());
            lhmm1.printAlignment();
        }
    };

    void print_options()
    {
        std::clog << "align v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     method   " << method << "\n";
        std::clog << "         [x] x        " << x << "\n";
        std::clog << "         [y] y        " << y << "\n";
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";

        std::clog << "what's wrong\n";
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
