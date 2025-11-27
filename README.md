<h2>Generating Smooth 3-Polytopes</h2>
<p>This repository contains the source code and generated dataset of my master's thesis "Generating Smooth 3-Polytopes", which will be available online soon. Implementing and running a noval algorithm for generating smooth 3-polytopes took about half a day to extend the pre-existing database of 103 smooth 3-polytopes with <=16 lattice points (by Lundman) to 12589 smooth 3-polytopes with <=50 lattice points (available in the file <code>prune_output</code>). Oda's Conjecture was verified on all generated polytopes. Moreover, all polytopes were found to have unimodular triangulations (see the repository <a href="https://github.com/krhuang/UniTriSat">UniTriSat</a>).</p>

<p>The theoretical background of the algorithm is explained in Section 3 of the thesis. The algorithm and its exact implementation is described in Sections 4 and 5. Detailed references are also contained therein.</p>

<h2>Where is the Dataset?</h2>
<p>The generated dataset of polytopes is available in the file <code>prune_output</code>.</p>

<p>Each smooth 3-polytope is given by a series of lines. The first line is the number of vertices, following by lines of triples of integers denoting the coordinates, such as follows. </p>

```
4 
0 0 0
1 0 0
0 1 0
0 0 1
```

<h3>Tables</h3>
The distribution of various properties (obtained via running <code>julia postprocess.jl</code> is shown below: 
<pre>
Distribution of Number of Vertices
-------------------------------------
               # Vertices | Frequency
-------------------------------------
                        4 | 4
                        6 | 3108
                        8 | 2514
                       10 | 1471
                       12 | 2788
                       14 | 1527
                       16 | 983
                       18 | 150
                       20 | 37
                       22 | 5
                       24 | 2
-------------------------------------
Distribution of Total Number of Lattice Points
-------------------------------------
   # Total Lattice Points | Frequency
-------------------------------------
                        4 | 1
                        6 | 1
                        7 | 1
                        8 | 3
                        9 | 4
                       10 | 6
                       11 | 5
                       12 | 12
                       13 | 10
                       14 | 17
                       15 | 14
                       16 | 29
                       17 | 21
                       18 | 39
                       19 | 30
                       20 | 54
                       21 | 42
                       22 | 63
                       23 | 56
                       24 | 94
                       25 | 75
                       26 | 113
                       27 | 91
                       28 | 154
                       29 | 103
                       30 | 186
                       31 | 140
                       32 | 247
                       33 | 180
                       34 | 292
                       35 | 227
                       36 | 353
                       37 | 256
                       38 | 411
                       39 | 335
                       40 | 541
                       41 | 431
                       42 | 651
                       43 | 556
                       44 | 817
                       45 | 652
                       46 | 903
                       47 | 772
                       48 | 1188
                       49 | 978
                       50 | 1435
-------------------------------------
Distribution of Interior Lattice Points
-------------------------------------
# Interior Lattice Points | Frequency
-------------------------------------
                        0 | 10519
                        1 | 61
                        2 | 227
                        3 | 607
                        4 | 564
                        5 | 321
                        6 | 221
                        7 | 55
                        8 | 13
                        9 | 1
-------------------------------------
</pre>

<h2>Running the Computation Yourself</h2>
<p>The algorithm utilizes other mathematical computation projects, namely plantri, Bohnert and Springer's large database of lattice polygons, and the many packages included in polymake. </p>

<p>For your convenience, plantri output and Bohnert-Springer's database have been included inside the plantri_output folder and the file SmoothPolytopesDB_test.txt, respectively.</p>

<p>All that remains is to run SmoothGeneration, and then prune the list for isomorphic copies and polytopes with too many lattice points, via the included Pruning Jupyter Notebook.</p>

<p>You can get started via the following commands:</p>

```

git clone https://github.com/krhuang/SmoothGeneration.git
cd SmoothGeneration

g++ -O3 -std=c++17 -Wall -fopenmp -I. SmoothGeneration.cpp helper_functions.cpp -o SmoothGeneration && ./SmoothGeneration

```

<p>If you wish to rerun the computation with smaller or larger numbers you will have to take care to increase the database of smooth polygons, have the needed plantri files, and modify various global variables in the programs.</p>

<p>If you plan to redo the computation with larger numbers please message me! I likely am either in the process of doing so myself or at least have some ideas for computational optimizations for you.</p>



