## inner - a MultiObjective Linear Program (MOLP) solver

A MOLP is a linear program with more than one objective function. The linear
*constraints* are arranged into a matrix form with c columns and r rows.
Columns correspond to variables x<sub>1</sub>, x<sub>2</sub>, . . . ,
x<sub>c</sub>, which are subject to the restrictions

<table><tbody><tr>
<td>L<sub>i</sub> &le; x<sub>i</sub> &le; U<sub>i</sub></td>
<td>where 1 &le; i &le; c</td>
</tr></tbody></table>

The lower bound L<sub>i</sub> can be -&#x221e; and the upper bound can be +&#x221e;.

For each row 1 &le; j &le; r the j-th constraint is

<table><tbody><tr>
<td>b<sub>j</sub> = a<sub>1,j</sub> x<sub>1</sub> + a<sub>2,j</sub> x<sub>2</sub> + . . . + a<sub>c,j</sub> x<sub>c</sub></td>
<td>l<sub>j</sub> &le; b<sub>j</sub> &le; u<sub>j</sub></td>
<td>where 1 &le; j &le; r</td>
</tr></tbody></table>

A *feasible solution* is a tuple **x** = &lt;x<sub>1</sub>, x<sub>2</sub>, . . . , x<sub>c</sub>&gt; 
which satisfies all constraints.

There are m &ge; 1 *objectives*, which are given as

<table><tbody><tr>
<td>y<sub>k</sub> = o<sub>1,k</sub> x<sub>1</sub> + o<sub>2,k</sub> x<sub>2</sub> + . . . + o<sub>c,k</sub> x<sub>c</sub></td>
<td>where 1 &le; k &le; m</td>
</tr><tbody></table>

The m-dimensional vector **y** = &lt;y<sub>1</sub>, . . ., y<sub>m</sub>&gt; is 
*achievable* if there is a feasible solution **x** which yields exactly these objective values.

The *solution of a minimizing* MOLP is a minimal list of achievable
objective vectors &ndash; called *extremal solutions* &ndash; such that any
other achievable objective vector is coordinatewise &ge; than some
non-negative linear combination of extremal solutions.

The solution set is always unique, and the MOLP solver's task is to find it. 
Extremal solutions are called *vertices* as they form the vertices of an
(unbounded) convex m-dimensional polytope.

This MOLP solver finds the extremal solutions by *iterations*. In each 
iteration step one more extremal solution is added to the final list. The
time required for an iteration varies widely from a couple of microseconds
to several days. After each iteration the solver checks if the process
was interrupted. If yes, it switches to a quick and dirty method which
might find further extremal solutions (but not necessarily all of them).


#### RESTRICTIONS

The linear constraints must have a feasible solution. Moreover, in a
minimizing MOLP all objectives must be bounded from below, and in a
maximizing MOLP all objectives must be bounded from above.


#### USAGE

The program is invoked as

    inner [options] <vlp-file>

The only obligatory argument is the file name which contains the description
of the problem in vlp format. Accepted options are

| Option | meaning |
|:-------|:--------|
| `-h`          | display a short help and quit |
| `--help`      | display all options |
| `--help=vlp`  | describe the vlp file format |
| `--help=out`  | describe the output format |
| `--version`   | version and copyright information |
| `--dump`      | dump the default config file and quit |
| `--config=<config-file>` or <br> `-c <config-file>`  | read configuration from the given file <br> use `--dump` to show the default config file |
| `-o <file>`  | save result (both vertices and facets) to \<file\> |
| `-ov <file>` | save vertices to \<file\> |
| `-of <file>` | save facets to \<file\> |
| `--name=NAME` or <br> `-n NAME`    | specify the problem name |
| `-m[0..3]`   | set message level: 0: none, 1: errors, 2: all, 3: verbose |
| `-q`         | quiet, same as `-m0`. Implies `--PrintStatistics=0` |
| `-p T`       | progress report in every T seconds (default: T=5) |
| `-p 0`       | no progress report |
| `-y+`        | report vertices immediately when generated (default) |
| `-y-`        | do not report vertices when generated |
| `--KEYWORD=value` | change value of a config keyword |


#### COMPILATION

The program uses a patched version of 'glpk', the GNU Linear Program Kit. 
First, glpk should be compiled after the patch has been applied. Unpack the
glpk source. Change to the '<GLPK>/src' directory, and execute the command

    patch -p1 < <INNER>/patch-X.Y.txt

where `<INNER>` is the inner directory. Then change to '<GLPK>`, and run 
the 'configure' script and then run 'make':

    ./configure
    ./make -CFLAGS='-O3 -DCSL'

You must define `CSL` as all patches to glpk are encapsulated in '#ifdef CSL`
blocks.

Going back to `<INNER>`, the following command complies a static version of
this program with name <NAME>:

    gcc -O3 -W -I <GLPK>/src -o <NAME> -DPROG=<NAME> *.c <GLPK>/src/.libs/libglpk.a -lm

#### EXAMPLES

The 'examples' directory contains vlp files describing different MOLPs. File
names contain three numbers separated by dashes: the number
of *objectives*, the number of *columns*, and the number of *rows*.  Each
file starts with comment lines describing the problem.

Solutions are in the 'solution' directory. The same file name with extension `.res`
contains the list of vertices and facets. `.out` files contain progress report,
statistics, and parameter settings.


#### AUTHOR

Laszlo Csirmaz, <csirmaz@ceu.edu>

#### DATE

April 10, 2016

