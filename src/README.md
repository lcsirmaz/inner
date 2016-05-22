## Source files

Special thanks to [https://github.com/csirmaz/inner](Elod Csirmaz) for supplying the thread
routines.

* glp_oracle.c, glp_oracle.h &ndash; implementing the vertex separation oracle based on glpk library
* inner.c, inner.h &ndash; the main loop executing the inner approximation algorithm
* main.c, main.h &ndash; the main program
* params.c, params.h &ndash; reading and processing configuration parameters and comand line options
* patch-4.57.txt, patch-4.60.txt &ndash; patches to the glpk library
* poly.c, poly.h &ndash; polytope algorithms implementing the double description vertex enumeration
* report.c, report.h &ndash; all output, reporting, saving results
* version.h &ndash; version information

