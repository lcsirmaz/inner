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

#### CONFIGURATION PARAMETERS

Fine tuning the algorithm, the underlying scalar LP solver and 
specifying the amount and type of saved information is done by giving
values of several keywords. Each keyword has a default value, which is
overwritten by the values in the config file (if specified), and those
values are overwritten by the `--KEYWORD=value` command line options.
Change tolerances with great care.

| Algorithm parameters | |
|:--------|:------------|
|`RandomFacet=0`<br>&nbsp; | 0 = no, 1 = yes <br>  pick a random facet which is then passed to the oracle. |
|`ExactFacetEq=0`<br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br>  when a facet is created, recompute its equation immediately from the set of adjacent <br> vertices. |
|`RecalculateFacets=100`<br>&nbsp;<br>&nbsp; | non-negative integers <br> after that many iterations recalculate all facet equations from the set of its adjacent vertices. The number should be zero (meaning never), or at least 5. |
|`CheckConsistency=0`<br>&nbsp;<br>&nbsp; | non-negative integer <br> after that many iterations check the consistency of the data structure against numerical errors. The number should be zero (meaning never), or at least 5. |
|`ExtractAfterBreak=1`<br>&nbsp;<br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> when the program is interrupted by Ctrl+C, continue extracting new vertices by asking <br> the oracle about every facet of the actual approximating polyhedron. Second Ctrl+C <br> aborts this post-processing. |
|**Oracle parameters**| |
|`OracleMessage=1`<br>&nbsp; | 0 = quiet, 1 = error, 2 = on, 3 = verbose <br> oracle (glpk) message level. |
|`OracleMethod=0`<br>&nbsp;  | 0 = primal, 1 = dual <br> the LP method used by the oracle. |
|`OraclePricing=1`<br>&nbsp; | 0 = standard, 1 = steepest edge <br> the LP pricing method. |
|`OracleRatioTest=1`<br>&nbsp; | 0 = standard, 1 = Harris' two pass <br> the LP ratio test. |
|`OracleTimeLimit=20` <br>&nbsp; | non-negative integer <br> time limit for each oracle call in seconds; 0 means unlimited. |
|`OracleItLimit=10000` <br>&nbsp; | non-negative integer <br> iteration limit for each oracle call; 0 means unlimited. |
|`OracleScale=1` <br>&nbsp; | 0 = no, 1 = yes <br> scale the constraint matrix; helps numerical stability. |
|`ShuffleMatrix=1` <br>&nbsp; | 0 = no, 1 = yes <br> shuffle rows and columns of the constraint matrix randomly. |
|`RoundVertices=1` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> when the oracle report a result vertex, round its coordinates to the nearest rational <br> number with small denominator. |
|**Reporting**| |
|`MessageLevel=3` <br>&nbsp;<br>&nbsp; | 0 = quiet, 1 = error, 2 = on, 3 = verbose <br> report level; quiet means no messages at all. Command line option `-m[0..3]` overrides <br> this value. |
|`Progressreport=5` <br>&nbsp;<br>&nbsp; | non-negative integer <br> minimum time between two progress reports in seconds. Should be zero for no progress <br> report, or at least 5. Use command line option `-p T` to override this value. |
|`VertexReport=1` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> report each vertex (extremal solution) immediately as it is found. Use command line <br> option `-y-` (no) or `-y+` (yes) to override the value defined here. |
|`MemoryReport=0` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> report the size and location of the 11 memory blocks storing the combinatorial data <br> structure whenever it changes. |
|`VertexAsFraction=1` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> if possible, print (and save) vertex coordinates as fraction with small denominator <br> rather than floating point numerals. |
|`PrintStistics=1` <br>&nbsp; | 0 = no, 1 = yes <br> print out resources used (number of iterations, ridge tests, etc.) when the program stops. |
|`PrintParams=0` <br>&nbsp; | 0 = no, 1 = yes <br> print out algorithm parameters which are not equal to their default values. |
|`PrintVertices=2` <br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> print out (again) all known vertices when the program terminates. |
|`PrintFacets=0` <br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> print out all known (relevant) facets when the program terminates. |
|`SaveVertices=2` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> when the program terminates, save known vertices to the file specified after command <br> line option `-o`. For file specified after `-ov` both 0 and 1 means &quot;save on normal exit only&quot;. |
|`SaveFacets=1` <br>&nbsp;<br>&nbsp;<br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> when the program terminates, save known (relevant) facets to the file specified after <br> command line option `-o`. For file specified after `-of` both 0 and 1 means &quot;save on <br> normal exit only&quot;. |
|**Tolerances**| |
|`PolytopeEps=1.3e-8` <br>&nbsp; | positive real number <br> a facet and a vertex are considered adjacent if their distance is smaller than this value. |
|`ScaleEps=3e-9` <br>&nbsp;<br>&nbsp; | positive real number <br> coefficients in the scaled facet equation are rounded to the nearest integer if they are <br> closer to it than this value. |
|`LineqEps=8e-8` <br>&nbsp;<br>&nbsp; | positive real number <br> when solving a system of linear equations for a facet equation, a coefficient smaller <br> than this is considered to be zero. |
|`RoundEps=1e-9` <br>&nbsp;<br>&nbsp; | positive real number <br> if the oracle reports vertices with rounded coordinates (`RoundVertices=1`), this is the <br> tolerance in the rounding algorithm. |
|`FacetRecalcEps=1e-6` <br>&nbsp;<br>&nbsp; | positive real number <br> when recalculating facet equation, report numerical instability if the new and old <br> coordinates differ at least that much. |


#### COMPILATION

The program uses a patched version of glpk, the GNU Linear Program Kit. 
First, glpk should be compiled after the patch has been applied. Unpack the
glpk source. Change to the `glpk-X.Y` directory, and execute the command

    patch -p1 < ../patch-X.Y.txt

assuming you have unpacked glpk in your `INNER` the directory. Then change to
`glpk-X.Y` and run 'configure' and 'make' as follows:

    ./configure
    ./make CFLAGS='-O3 -DCSL'

You must define `CSL` as all patches to glpk are encapsulated in `#ifdef CSL`
blocks.

Going back to `INNER`, the following command compiles a static version of
this program with name NAME:

    gcc -O3 -W -I glpk-X.Y/src -o NAME -DPROG=NAME *.c glpk-X.Y/src/.libs/libglpk.a -lm

#### EXAMPLES

The 'examples' directory contains vlp files describing different MOLPs. File
names contain three numbers separated by dashes: the number of *objectives*,
the number of *rows*, and the number of *columns*.  Each file starts with
some comment lines.

Solutions are in the 'solution' directory. The same file name with extension
`.res`contains the list of vertices and facets.  `.out` files contain
progress report, statistics, and parameter settings.


#### AUTHOR

Laszlo Csirmaz, <csirmaz@ceu.edu>

#### DATE

April 10, 2016

