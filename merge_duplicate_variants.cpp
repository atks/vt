/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>

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

#include "merge_duplicate_variants.h"

int merge_duplicate_variants(int argc, char ** argv)
{
	//options
	std::string ivcf_file;
    std::string ovcf_file;
    bool merge_by_pos;

	try
	{
		std::string desc =
"Merges duplicate variants by position with the option of considering alleles.  (This just discards the duplicate variant that appears later in the VCF file)\n\
$path = /net/fantasia/home/atks/programs/vt\n\
e.g. $path/vt merge_duplicate_variants -i $path/test/8904indels.dups.genotypes.vcf -o out.vcf\n\
e.g. $path/vt merge_duplicate_variants -p -i $path/test/8904indels.dups.genotypes.vcf -o out.vcf\n";

   		std::string version = "0.5";
		TCLAP::CmdLine cmd(desc, ' ', version);
		TCLAP::ValueArg<std::string> arg_ovcf_file("o", "output-vcf", "Output VCF file [-]", false, "-", "string", cmd);
		TCLAP::SwitchArg arg_merge_by_position("p", "merge-by-position", "Merge by position [false]", cmd, false);
        TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("input-vcf-files", "Input VCF Files", true, "string", cmd);
    
		cmd.parse(argc, argv);

        if (arg_input_vcf_files.getValue().size()==1)
		{
		    ivcf_file = arg_input_vcf_files.getValue()[0];
		}
		else
	    {
	        //error
	    }
		ovcf_file = arg_ovcf_file.getValue();
		merge_by_pos = arg_merge_by_position.getValue();

		std::clog << "merge_duplicate_variants v0.57\n\n";

		std::clog << "Options: [_] Input VCF File    " << ivcf_file << "\n";
		std::clog << "         [o] Output VCF File   " << ovcf_file << "\n";
        std::clog << "         [p] Merge by          " << (merge_by_pos?"position":"position and alleles") << "\n\n";

	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		abort();
	}

    bcf_srs_t *sr;
    sr =  bcf_sr_init();
    if (!bcf_sr_add_reader(sr, ivcf_file.c_str()))
            std::cerr << "Failed to open or the file not indexed: " << ivcf_file << "\n";

    bcf_hdr_t *hdr = sr->readers[0].header;
    bcf_hdr_fmt_text(hdr);

    htsFile *out = hts_open(ovcf_file.c_str(), "w", 0);
    bcf_add_hs37d5_contig_headers(hdr);
    vcf_hdr_write(out, hdr);
    std::map<std::string, uint32_t> variants;
    std::stringstream ss;
    uint32_t no_unique_variants = 0;
    uint32_t no_total_variants = 0;
    int32_t last_rid = -1;
    int32_t last_pos0 = -1;

    while (bcf_sr_next_line(sr))
    {
        bcf1_t *line = sr->readers[0].buffer[0];
        bcf_unpack(line, BCF_UN_STR);

        if (line->rid==last_rid && line->pos==last_pos0)
        {
            if (!merge_by_pos)
            {            
                ss.str("");
                for (uint32_t i=0; i<line->n_allele; ++i)
                {
                    ss << (i?":":"")  << (line->d).allele[i];
                }
            
                if (variants.find(ss.str()) == variants.end())
                {
                    variants[ss.str()] = 1;
                    vcf_write1(out, hdr, line);
                    ++no_unique_variants;
                }
            }
        }
        else
        {
            variants.clear();
            
            last_rid = line->rid;
            last_pos0 = line->pos;
                       
            if (!merge_by_pos)
            {    
                ss.str("");
                for (uint32_t i=0; i<line->n_allele; ++i)
                {
                    ss << (i?":":"")  << (line->d).allele[i];
                }
                
                variants[ss.str()] = 1;
            }
            
            vcf_write1(out, hdr, line);
            ++no_unique_variants;
        }

        ++no_total_variants;
    }

    hts_close(out);
    bcf_sr_destroy(sr);

    std::cerr << "Stats: Total Number of Observed Variants   " << no_total_variants << "\n";
    std::cerr << "       Total Number of Unique Variants     " << no_unique_variants << "\n\n";

    return 0;
};
