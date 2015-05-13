# ---------------------------------------------------------------------------
# Function strip_stderr()
# ---------------------------------------------------------------------------
# Utility that removes the time from vt stderr output.

strip_stderr()
{
    sed 's/time.*?s/time <stripped>/g' \
        | sed 's/input VCF file.*/input VCF file <stripped>/g'
}
export -f strip_stderr
