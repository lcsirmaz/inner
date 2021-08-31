## Source files

Special thanks to [Elod P. Csirmaz](https://github.com/csirmaz/inner) for the thread routines.

* [data.c](data.c), [data.h](data.h) &ndash; parsing and reading character input
* [glp_oracle.c](glp_oracle.c), [glp_oracle.h](glp_oracle.h) &ndash; implementing the vertex separation oracle based on glpk library
* [inner.c](inner.c), [inner.h](inner.h) &ndash; the main loop executing the inner approximation algorithm
* [main.c](main.c), [main.h](main.h) &ndash; the main program
* [params.c](params.c), [params.h](params.h) &ndash; reading and processing configuration parameters and comand line options
* [patch-4.57.txt](patch-4.57.txt), [patch-4.60.txt](patch-4.60.txt), [patch-4.65.txt](patch-4.65.txt) &ndash; patches to the glpk library
* [poly.c](poly.c), [poly.h](poly.h) &ndash; polytope algorithms implementing the double description vertex enumeration
* [report.c](report.c), [report.h](report.h) &ndash; all output, reporting, saving results
* [version.h](version.h) &ndash; version information

