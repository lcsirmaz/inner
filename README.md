## inner - a Multi Objective Linear Program (MOLP) solver

A MOLP is a linear program with more than one objective function. The linear
'constraints' are typically written in a matrix form where the columns
correspond to (constraint) variables. 


#### USAGE

The program is invoked as

    inner [options] <vlp-file>

The only obligatory argument is the file name which contains the descritpion
of the problem is vlp format. Accepted options are

| Option | meaning |
|:-------|:--------|
| -h          | print a short help and quit |
| --help      | print all options |
| --help=vlp  | describe the vlp file format |
| --help=out  | describe output format |
| --version   | version and copyright information |
| --dump     | dump the default config file and quit |
| --config=<config-file> | |
| -c <config-file>  | read configuration from the given file.
                      Use \'--dump\' to show the default config file |
| -o <file>  | save result (both vertices and facets) to <file> |
| -ov <file> | save vertices to <file> |
| -of <file> | save facets to <file> |
| --name=NAME | |
| -n NAME    | specify the problem name |
| -m[0..3]   | set message level: 0: none, 1: errors, 2: all, 3: verbose |
| -q         | quiet, same as -m0 |
| -p T       | progress report in every T seconds (default: T=5) |
| -p0        | no progress report |
| -y+        | report vertices immediately when generated (default) |
| -y-        | do not report vertices when generated |
| -r N       | recalculate facet equations after N rounds (default: N=100) |
| -k N       | check numerical consistency after N rounds (default: N=0) |
| --KEYWORD=value | change value of a config keyword |



#### COMPILATION

The program uses a patched version of 'glpk', the GNU Linear Program Kit. 
First, glpk should be compiled after the patch has been applied (Version
4.57 only). `<GLPK>` is the main directory of the glpk distribution, and 
`<NAME>` is the complied program's name. To compile it using static version
of the patched glpk, use the command

    gcc -O3 -W -I <GLPK>/src -o <NAME> -DPROG=<NAME> *.c <GLPK>/src/.libs/libglpk.a -lm



#### AUTHOR

Laszlo Csirmaz, <csirmaz@ceu.edu>

#### DATE

April 7, 2016

