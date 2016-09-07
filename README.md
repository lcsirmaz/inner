# inner - a MultiObjective Linear Program (MOLP) solver

A MOLP is a linear program with more than one objective functions. The linear
*constraints* are arranged in a matrix form with c columns and r rows.
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
<td>where 1 &le; k &le; m.</td>
</tr><tbody></table>

The m-dimensional objective vector **y** = &lt;y<sub>1</sub>, . . ., y<sub>m</sub>&gt; is 
*achievable* if there is a feasible solution **x** which yields exactly these objective values.

The *solution* of a MOLP is the list of all extremal objective vectors. For a 
minimizing (maximizing) MOLP the objective vector **y** is *extremal* if there 
is no other achievable objective vector **y**' which would be coordinatewise &le; 
**y** (or &ge; when maximizing).

The task of MOLP solver is to find the collection of all extremal objective vectors.
These extremal vectors are called *vertices* as they form the vertices of an
(unbounded) convex m-dimensional polytope.

This MOLP solver finds the extremal vectors by *iteration*. In each 
step one more extremal vector is added to the final list. The
time required for an iteration varies widely from a couple of microseconds
to several days. After each iteration the solver checks if the process
has been interrupted by an USR1 signal.
If it has, it switches to a quick and dirty method which
might find further extremal vectors (but not necessarily all of them).


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
| `-o <file>`  | save results (both vertices and facets) to \<file\> |
| `-ov <file>` | save vertices to \<file\> |
| `-of <file>` | save facets to \<file\> |
| `-oc <stub>` | file stub for checkpoint files |
| `--name=NAME` or <br> `-n NAME`    | specify the problem name |
| `--boot=<vertex-list>` | start the algorithm with these vertices |
| `--resume=<chk-file>` | resume computation from a checkpoint file |
| `-m[0..3]`   | set message level: 0: none, 1: errors, 2: all, 3: verbose |
| `-q`         | quiet, same as `-m0`. Implies `--PrintStatistics=0` |
| `-p T`       | progress report in every T seconds (default: T=5) |
| `-p 0`       | no progress report |
| `-y+`        | report extremal solutions (vertices) immediately when generated (default) |
| `-y-`        | do not report extremal solutions when generated |
| `--KEYWORD=value` | change value of a config keyword |

#### CONFIGURATION PARAMETERS

Fine tuning the algorithm and the underlying scalar LP solver, and 
specifying the amount and type of information saved is done by giving
values of several keywords. Each keyword has a default value, which is
overwritten by the values in the config file (if specified), and those
values are overwritten by the `--KEYWORD=value` command line options.
Change tolerances with great care.

|**Algorithm parameters**| |
|:--------|:------------|
|`RandomFacet=1`<br>&nbsp; | 0 = no, 1 = yes <br>  pick a random facet which is then passed to the oracle. |
|`ExactFacetEq=0`<br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br>  when a facet is created, recompute its equation immediately from the set of adjacent <br> vertices. |
|`RecalculateFacets=100`<br>&nbsp;<br>&nbsp; | non-negative integer <br> after this many iterations recalculate all facet equations from the set of its adjacent vertices. The number should be zero (meaning never), or at least 5. |
|`CheckConsistency=0`<br>&nbsp;<br>&nbsp; | non-negative integer <br> after this many iterations check the consistency of the data structure against numerical errors. The number should be zero (meaning never), or at least 5. |
|`ExtractAfterBreak=1`<br>&nbsp;<br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> when the program receives a `SIGUSR1` signal, continue extracting new vertices by asking<br> the oracle about every facet of the actual approximating polyhedron. Second signal <br> aborts this post-processing. |
|`VertexPoolSize=0`<br>&nbsp;<br>&nbsp;<br>&nbsp; | non-negative integer <br> size of the vertex pool: add the vertex to the approximation which discards the largest <br> number of existing facets. Should be zero (don't use it), or at least 5. Using vertex pool <br> adds more oracle calls, but can simplify the approximating polytope. |
|`CheckPoint=10000`<br>&nbsp;<br>&nbsp;<br>&nbsp;<br>&nbsp; | positive integer <br> frequency (in seconds) for creating checkpoint files when the option <br> `-oc <filestub>` is given. The filename is got from the stub by appending `NNN.chk` where `NNN` starts <br> with 000 and increases. The computation can be resumed from a checkpoint file <br> by calling inner with `--resume=<checkpoint>` The value should be more than 500. |
|`OracleCallLimit=1`<br>&nbsp;<br>&nbsp; | non-negative integer <br> the maximal number of unsuccessful oracle calls during an iteration when filling the <br> vertex pool. Zero means no limit; otherwise should be less than 100. |
|`Threads=0`<br>&nbsp;<br>&nbsp; | non-negative integer <br> number of threads to use; should be less than 64. Zero means use as many threads as <br> are available; 1 means don't use threads. |
|**Oracle parameters**| |
|`OracleMessage=1`<br>&nbsp; | 0 = quiet, 1 = error, 2 = on, 3 = verbose <br> oracle (glpk) message level. |
|`OracleMethod=0`<br>&nbsp;  | 0 = primal, 1 = dual <br> the LP method used by the oracle. |
|`OraclePricing=1`<br>&nbsp; | 0 = standard, 1 = steepest edge <br> the LP pricing method. |
|`OracleRatioTest=1`<br>&nbsp; | 0 = standard, 1 = Harris' two pass <br> the LP ratio test. |
|`OracleTimeLimit=20` <br>&nbsp; | non-negative integer <br> time limit for each oracle call in seconds; 0 means unlimited. |
|`OracleItLimit=10000` <br>&nbsp; | non-negative integer <br> iteration limit for each oracle call; 0 means unlimited. |
|`OracleScale=1` <br>&nbsp; | 0 = no, 1 = yes <br> scale the constraint matrix; helps numerical stability. |
|`ShuffleMatrix=1` <br>&nbsp; | 0 = no, 1 = yes <br> shuffle rows and columns of the constraint matrix randomly. |
|`RoundVertices=1` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> when the oracle reports a result vertex, round its coordinates to the nearest rational <br> number with small denominator. |
|**Reporting**| |
|`MessageLevel=2` <br>&nbsp;<br>&nbsp; | 0 = quiet, 1 = error, 2 = on, 3 = verbose <br> report level; quiet means no messages at all. The command line option `-m[0..3]` overrides <br> this value. |
|`Progressreport=5` <br>&nbsp;<br>&nbsp; | non-negative integer <br> minimum time between two progress reports in seconds. Should be zero for no progress <br> reports, or at least 5. The command line option `-p T` overrides this value. |
|`VertexReport=1` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> report each vertex (extremal solution) immediately when it is found. The command line <br> option `-y-` (no) or `-y+` (yes) overrides the value defined here. |
|`MemoryReport=0` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> report the size and location, whenever they change, of the 11 memory blocks storing <br> the combinatorial data structure. |
|`VertexAsFraction=1` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> if possible, print (and save) vertex coordinates as fractions with small denominators <br> rather than floating point numerals. |
|`PrintStatistics=1` <br>&nbsp; | 0 = no, 1 = yes <br> print resources used (number of iterations, ridge tests, etc.) when the program stops. |
|`PrintParams=1` <br>&nbsp; | 0 = no, 1 = yes <br> print algorithm parameters which are not equal to their default values. |
|`PrintVertices=2` <br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> print (again) all known vertices when the program terminates. |
|`PrintFacets=0` <br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> print all known (relevant) facets when the program terminates. |
|`SaveVertices=2` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> when the program terminates, save known vertices to the file specified after command <br> line option `-o`. For file specified after `-ov` both 0 and 1 means &quot;save on normal exit only&quot;. |
|`SaveFacets=1` <br>&nbsp;<br>&nbsp;<br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> when the program terminates, save known (relevant) facets to the file specified after <br> the command line option `-o`. For the file specified after `-of` both 0 and 1 means &quot;save on <br> normal exit only&quot;. |
|**Tolerances**|**Change these values with care**|
|`PolytopeEps=1.3e-8` <br>&nbsp; | positive real number <br> a facet and a vertex are considered adjacent if their distance is smaller than this value. |
|`ScaleEps=3e-9` <br>&nbsp;<br>&nbsp; | positive real number <br> coefficients in the scaled facet equation are rounded to the nearest integer if they are <br> closer to it than this value. |
|`LineqEps=8e-8` <br>&nbsp;<br>&nbsp; | positive real number <br> when solving a system of linear equations for a facet equation, a coefficient smaller <br> than this is considered to be zero. |
|`RoundEps=1e-9` <br>&nbsp;<br>&nbsp; | positive real number <br> if the oracle reports vertices with rounded coordinates (`RoundVertices=1`), this is the <br> tolerance in the rounding algorithm. |
|`FacetRecalcEps=1e-6` <br>&nbsp;<br>&nbsp; | positive real number <br> when recalculating facet equation, report numerical instability if the new and old <br> coordinates differ at least this much. |

#### EXIT STATUS

When the program terminates, the exit status is zero when the problem has
been solved successfully; otherwise it indicates the failure condition:

|**Exit value**| | 
|:------|:-----|
|0    | problem solved |
|1    | error before starting the algorithm, such as missing argument, wrong or missing data file, out of memory, etc. |
|2    | there is no feasible solution for the MOLP problem |
|3    | the solution space is unbounded in some objective direction |
|4    | error while executing the algorithm (oracle failure, numerical error, etc) |
|5    | program execution was interrupted by using a `USR1` signal |


#### SIGNALS

When the program receives a `USR1` signal, it stops the main cycle of
iterations, and switches to a &quot;quick and dirt&quot; method to generate
additional extremal solutions. (Actually, the oracle is asked about all
facets of the actual approximation; it returns an extremal solution, if
exists, on or beyond that facet.) The method might miss extremal solutions,
so the result is not known (can need not be) complete. A second `USR1`
signal aborts this post-processing.

When receiving a `USR2` signal (since version 2.8), the program dumps all
vertices and facets of the actual approximation. Similarly to checkpoint
files, the dump file can be used to resume the computation from that point.

The name of the dump file is created as follows. If a checkpoint stub is provided
after the `-oc` argument, then `000.dmp` is appended to the stub. If no
`-oc` argument is provided but there is a `-o` output file, then the 
extension of that filename is replaced by `dmp`. If neither `-oc` nor `-o`
argument is present, then no dump file is created. When receiving a second
`USR2` signal, the same filename is used: the previous content is silently
overwritten.

Signals are checked only when an iteration step has been completed, that is a
new vertex has been added to the approximation. (These are the points when
the data is guaranteed to be consistent.) Depending on the complexity of the
problem and the number of intermediate vertices and facets, the time
between these points vary from milliseconds to several days.

#### USING MULTIPLE CORES

The threaded version of `inner` can use all available cores on the computer 
to speed up the computation. Multiple cores can help in the combinatorial part, 
when the work is distributed approximately uniformly among them. The LP
part (oracle calls) is inherently single-threaded and cannot be speed up
this way. Thus LP intensive problems (large number of columns and rows
and only a few objectives) will not gain too much from multiple threads.
Specifying more cores than available comes with a significant performance
penalty. 


#### COMPILATION

The program uses a patched version of glpk, the GNU Linear Program Kit. 
First, glpk should be compiled after the patch has been applied. Unpack the
glpk source. In the `glpk-X.Y` directory execute the command

    patch -p1 < ../src/patch-X.Y.txt

assuming you have unpacked glpk in the root of this repository. Still in the
`glpk-X.Y` directory run 'configure' and 'make' as follows:

    ./configure
    ./make CFLAGS='-O3 -DCSL'

You must define `CSL` as all patches to glpk are encapsulated in `#ifdef CSL`
blocks.

Changing to the directory `INNER/src`, the following command compiles
**inner** linking the patched glpk routines statically:

    gcc -O3 -W -I ../glpk-X.Y/src -o inner *.c ../glpk-X.Y/src/.libs/libglpk.a -lm

The threaded version requires defining `USETHREADS` and adding the flag `-pthread`:

    gcc -O3 -W -I ../glpk-X.Y/src -o innerth -DUSETHREADS -pthread *.c ../glpk-X.Y/src/.libs/libglpk.a -lm

#### EXAMPLES

The 'examples' directory contains vlp files describing different MOLPs. File
names contain three numbers separated by dashes: the number of *objectives*,
the number of *rows*, and the number of *columns*.  Each file starts with
some comment lines.

Solutions are in the 'solution' directory. The same file name with the extension
`.res` contains the list of vertices and facets.  `.out` files contain
progress reports, statistics, and parameter settings.


#### AUTHOR

Laszlo Csirmaz, <csirmaz@ceu.edu>

#### DATE

April 10, 2016

