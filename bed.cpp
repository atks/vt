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

#include "bed.h"

/**
 * Constructor.
 */
BEDRecord::BEDRecord(kstring_t *s)
{
    std::vector<std::string> fields;
    split(fields, "\t", s->s);

    chrom = fields[0];
    str2int32(fields[1], start1);
    str2int32(fields[2], end1);
};

/**
 * Constructor.
 */
BEDRecord::BEDRecord(char *s)
{
    std::vector<std::string> fields;
    split(fields, "\t", s);

    chrom = fields[0];
    str2int32(fields[1], start1);
    str2int32(fields[2], end1);
};

/**
 * Constructor.
 */
BEDRecord::BEDRecord(std::string& s)
{
    std::vector<std::string> fields;
    split(fields, "\t", s.c_str());

    chrom = fields[0];
    str2int32(fields[1], start1);
    str2int32(fields[2], end1);
};

/**
 * Constructor.
 */
BEDRecord::BEDRecord(std::string& chrom, int32_t start1, int32_t end1)
{
    this->chrom = chrom;
    this->start1 = start1;
    this->end1 = end1;
};

/**
 * Prints this BED record to STDERR.
 */
void BEDRecord::print()
{
    std::cerr << this->chrom << ":" << this->start1 << "-" <<this->end1 << "\n";
};

/**
 * String version of BED record.
 */
std::string BEDRecord::to_string()
{
    kstring_t s = {0,0,0};

    kputs(this->chrom.c_str(), &s);
    kputc(':', &s);
    kputw(this->start1, &s);
    kputc('-', &s);
    kputw(this->end1, &s);

    std::string str(s.s);
    if (s.m) free(s.s);
    return str;
};