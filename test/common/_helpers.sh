# ---------------------------------------------------------------------------
# Function strip_stderr()
# ---------------------------------------------------------------------------
# Utility that removes the time from vt stderr output.

strip_stderr()
{
    sed 's/Time elapsed.*/Time elapsed <stripped>/g' | sed 's/file.*/file <stripped>/g' 
}
export -f strip_stderr
