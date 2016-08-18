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

#include "program.h"

void VTOutput::failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e)
{
    std::clog << "\n";
    std::clog << "  " << e.what() << "\n\n";
    usage(c);
    exit(1);
}

void VTOutput::usage(TCLAP::CmdLineInterface& c)
{
    std::string s = "";
    std::list<TCLAP::Arg*> args = c.getArgList();
    //prints unlabeled arument list first
    for (TCLAP::ArgListIterator it = args.begin(); it != args.end(); it++)
    {
        TCLAP::Arg& arg = **it;
        if (typeid(arg)==typeid(TCLAP::UnlabeledValueArg<std::string>))
        {
            TCLAP::UnlabeledValueArg<std::string> *i = (TCLAP::UnlabeledValueArg<std::string> *) (*it);
            s = i->getName();
        }
        else if (typeid(arg)==typeid(TCLAP::UnlabeledMultiArg<std::string>))
        {
            TCLAP::UnlabeledMultiArg<std::string> *i = (TCLAP::UnlabeledMultiArg<std::string> *) (*it);
            s = i->getName();
        }
    }

    std::clog << c.getProgramName() << " v" << c.getVersion() << "\n\n";
    std::clog << "description : " << c.getMessage() << "\n\n";
    std::clog << "usage : vt "  << c.getProgramName() << " [options] " << s << "\n\n";

    //prints rest of arguments
    for (TCLAP::ArgListIterator it = args.begin(); it != args.end(); it++)
    {
        if (it==args.begin())
        {
            std::clog << "options : ";
        }
        else
        {
            std::clog << "          ";
        }

        TCLAP::Arg& arg = **it;
        if (typeid(arg)==typeid(TCLAP::ValueArg<std::string>) ||
            typeid(arg)==typeid(TCLAP::ValueArg<uint32_t>) ||
            typeid(arg)==typeid(TCLAP::ValueArg<int32_t>) ||
            typeid(arg)==typeid(TCLAP::ValueArg<double>) ||
            typeid(arg)==typeid(TCLAP::ValueArg<float>))
        {
            TCLAP::ValueArg<std::string> *i = (TCLAP::ValueArg<std::string> *) (*it);

            std::clog  << "-" << (i->getFlag()=="" ? i->getName() : i->getFlag())
                       << "  " << i->getDescription() << "\n";
        }
        else if (typeid(arg)==typeid(TCLAP::SwitchArg))
        {
            TCLAP::SwitchArg *i = (TCLAP::SwitchArg *) (*it);

            std::clog  << "-" << i->getFlag()
                       << "  " << i->getDescription() << "\n";
        }
        else if (typeid(arg)==typeid(TCLAP::UnlabeledValueArg<std::string>))
        {
            //ignored
        }
        else if (typeid(arg)==typeid(TCLAP::UnlabeledMultiArg<std::string>))
        {
            //ignored
        }
        else
        {
            std::clog << "oops, argument type not handled\n";
        }
    }

    std::clog  <<  "\n";
}

/**
 * Parse multiple files from command line unlabeled arguments or -L denoted file list.  If both are defined, the files are merged.
 *
 * @files          - file names are stored in this vector
 * @argument_files - vector of input files
 * @file_list      - file names stored in a file
 *
 */
void Program::parse_files(std::vector<std::string>& files, const std::vector<std::string>& arg_files, std::string file_list)
{
    files.clear();

    if (arg_files.size()!=0)
    {
        files = arg_files;
    }

    if (file_list != "")
    {
        htsFile *file = hts_open(file_list.c_str(), "r");
        if (file==NULL)
        {
            std::cerr << "cannot open " << file_list << "\n";
            exit(1);
        }
        kstring_t *s = &file->line;
        while (hts_getline(file, '\n', s) >= 0)
        {
            if (s->s[0]!='#')
            {
                files.push_back(std::string(s->s));
            }
        }
        hts_close(file);
    }
}

/**
 * Parse intervals. Processes the interval list first followed by the interval string. Duplicates are dropped.
 *
 * @intervals       - intervals stored in this vector
 * @interval_list   - file containing intervals
 * @interval_string - comma delimited intervals in a string
 *
 * todo: merge overlapping sites?
 */
void Program::parse_intervals(std::vector<GenomeInterval>& intervals, std::string interval_list, std::string interval_string)
{
    intervals.clear();
    std::map<std::string, uint32_t> m;

    if (interval_list!="")
    {
        htsFile *file = hts_open(interval_list.c_str(), "r");
        if (file)
        {
            kstring_t *s = &file->line;
            while (hts_getline(file, '\n', s)>=0)
            {
                std::string ss = std::string(s->s);
                if (m.find(ss)==m.end())
                {
                    m[ss] = 1;
                    GenomeInterval interval(ss);
                    intervals.push_back(interval);
                }
            }
            hts_close(file);
        }
    }

    std::vector<std::string> v;
    if (interval_string!="")
        split(v, ",", interval_string);

    for (size_t i=0; i<v.size(); ++i)
    {
        if (m.find(v[i])==m.end())
        {
            m[v[i]] = 1;
            GenomeInterval interval(v[i]);
            intervals.push_back(interval);
        }
    }
}

/**
 * Parse filters. Processes the filter list.
 *
 * @filters       - filters stored in this vector
 * @filter_string - comma delimited filters in a string
 * @n             - ensure that filters vector had n filters.
 *                  if there are less, just pad with empty strings
 *                  if there are more, thrown an error. 
 *                  if n is 0, ignore the previous contraints.
 * @pad           - if there are less than expected variant expressions
 *                      when true, the remaining filter expressions are padded with the empty string.
 *                      when false and only one expression is observed, the remaining filter expressions
 *                      duplicated with that filter expression.
 */
void Program::parse_filters(std::vector<std::string>& filters, std::string filter_string, int32_t n, bool pad)
{
    filters.clear();
    if (filter_string!="")
        split(filters, ",", filter_string);
    
    if (n && filters.size()!=0)
    {
        if (filters.size()<n)
        {
            if (pad)
            {    
                while(filters.size()!=n) filters.push_back("");
                fprintf(stderr, "[%s:%d %s] Number of filters less than expected, padding remaining filters with empty string\n", __FILE__, __LINE__, __FUNCTION__);
            }
            else
            {
                if (filters.size()==1)
                {
                    filters.resize(n, filters[0]);
                }    
                else
                {
                    fprintf(stderr, "[%s:%d %s] %d filter expressions are expected : %s\n", __FILE__, __LINE__, __FUNCTION__, n, filter_string.c_str());
                    exit(1);
                }
            }
        }
        else if (filters.size()>n)
        {
            fprintf(stderr, "[%s:%d %s] %d filter expressions are expected : %s\n", __FILE__, __LINE__, __FUNCTION__, n, filter_string.c_str());
            exit(1);
        }
        else
        {
            //all is good
        }
    }  
    
    if (filters.size()==0)
    {
        filters.push_back("");
    }     
}

/**
 * Parse a list of strings delimited by commas.
 *
 * @strings        - list of strings
 * @string_list    - comma delimited strings
 */
void Program::parse_string_list(std::vector<std::string>& strings, std::string string_list)
{
    strings.clear();
    if (string_list!="")
        split(strings, ",", string_list);
}

/**
 * Print reference FASTA file option.
 */
void Program::print_ref_op(const char* option_line, std::string ref_fasta_file)
{
    if (ref_fasta_file!="")
    {
        std::clog << option_line << ref_fasta_file << "\n";
    }
}

/**
 * Print string option, hide if not present.
 */
void Program::print_str_op(const char* option_line, std::string str_value)
{
    if (str_value!="")
    {
        std::clog << option_line << str_value << "\n";
    }
}

/**
 * Print number option, hide if 0.
 */
void Program::print_num_op(const char* option_line, uint32_t num_value)
{
    if (num_value)
    {
        std::clog << option_line << num_value << "\n";
    }
}

/**
 * Print switch option, hide if not switched on.
 */
void Program::print_boo_op(const char* option_line, bool value)
{
    if (value)
    {
        std::clog << option_line << "true" << "\n";
    }
    else
    {
        std::clog << option_line << "false" << "\n";
    }
}

/**
 * Print intervals option.
 */
void Program::print_int_op(const char* option_line, std::vector<GenomeInterval>& intervals)
{
    if (intervals.size()!=0)
    {
        std::clog << option_line;
        for (size_t i=0; i<std::min((uint32_t)intervals.size(),(uint32_t)5); ++i)
        {
            if (i) std::clog << ",";
            std::clog << intervals[i].to_string();
        }
        if (intervals.size()>5)
        {
            std::clog << " and " << (intervals.size()-5) <<  " other intervals\n";
        }
        else
        {
            std::clog << "\n";
        }
    }
}

/**
 * Print string vector.
 */
void Program::print_strvec(const char* option_line, std::vector<std::string>& vec)
{
    if (vec.size()!=0)
    {
        std::clog << option_line;
        for (size_t i=0; i<std::min((uint32_t)vec.size(),(uint32_t)4); ++i)
        {
            if (i) std::clog << ",";
            std::clog << vec[i];
        }
        
        if (vec.size()>4)
        {
            std::clog << " and " << (vec.size()-4) <<  " other values\n";
        }
        else
        {
            std::clog << "\n";
        }
    }
}

/**
 * Print input files.
 */
void Program::print_ifiles(const char* option_line, std::vector<std::string>& files)
{
    if (files.size()!=0)
    {
        std::clog << option_line;
        for (size_t i=0; i<std::min((uint32_t)files.size(),(uint32_t)2); ++i)
        {
            if (i) std::clog << ",";
            std::clog << files[i];
        }
        if (files.size()>2)
        {
            std::clog << " and " << (files.size()-2) <<  " other files\n";
        }
        else
        {
            std::clog << "\n";
        }
    }
}

/**
 * Parse samples. Processes the sample list. Duplicates are dropped.
 *
 * @nsamples     - number of unique samples found in list
 * @sample_list  - file containing sample names
 */
char** Program::read_sample_list(int32_t& nsamples, std::string sample_list)
{
    std::vector<std::string> vsamples;
    std::map<std::string, int32_t> map;

    if (sample_list!="")
    {
        htsFile *file = hts_open(sample_list.c_str(), "r");
        if (file)
        {
            kstring_t *s = &file->line;
            while (hts_getline(file, '\n', s)>=0)
            {
                std::string ss = std::string(s->s);
                if (map.find(ss)==map.end())
                {
                    map[ss] = 1;
                    vsamples.push_back(ss);
                }
            }
            hts_close(file);
        }

        nsamples = vsamples.size();
        char** samples = (char**) malloc(sizeof(char*)*nsamples);

        for (int32_t i=0; i<vsamples.size(); ++i)
        {
            samples[i] = strdup(vsamples[i].c_str());
        }

        return samples;
    }

    return NULL;
}