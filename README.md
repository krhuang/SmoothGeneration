<h2>Generating Smooth 3-Polytopes</h2>
<p>This repository contains the code and generated dataset of my master's thesis "Generating Smooth 3-Polytopes", which will be available online soon. Implementing and running a noval algorithm for generating smooth 3-polytopes took about half a day to extend the pre-existing database of 103 smooth 3-polytopes with <=16 lattice points (by Lundman) to 1364 smooth 3-polytopes with <=31 lattice points (available in the file prune_output). Oda's Conjecture was verified on all generated polytopes. There is ongoing computational work to extend this database to all smooth 3-polytopes with <=50 lattice points. </p>

<p>The theoretical background of the algorithm is explained in Section 3 of the thesis. The algorithm and its exact implementation is described in Sections 4 and 5. Detailed references are also contained therein.</p>

<h2>Where is the Dataset?</h2>
<p>The generated dataset of polytopes is available in the file `prune_output`.</p>

<p>Each smooth 3-polytope is given by a series of lines. The first line is the number of vertices, following by lines of triples of integers denoting the coordinates</p>

<h2>Running the Computation Yourself</h2>
<p>The algorithm utilizes other mathematical computation projects, namely plantri, Balletti's database of smooth polygons, and the many packages included in polymake. </p>

<p>For your convenience, plantri output and Balletti's database have been included inside the plantri_output folder and the file SmoothPolytopesDB.txt, respectively.</p>

<p>All that remains is to run SmoothGeneration, and then prune the list for isomorphic copies and polytopes with too many lattice points, via the included Pruning Jupyter Notebook.</p>

<p>You can get started via the following commands:</p>

```

git clone https://github.com/krhuang/SmoothGeneration.git

```

```

g++ -O3 -std=c++17 -Wall -fopenmp -I. SmoothGeneration.cpp helper_functions.cpp -o SmoothGeneration && ./SmoothGeneration

```

<p>If you wish to rerun the computation with smaller or larger numbers you will have to take care to increase the database of smooth polygons, have the needed plantri files, and modify various global variables in the programs.</p>

<p>If you plan to redo the computation with larger numbers please message me! I likely am either in the process of doing so myself or at least have many ideas for computational optimizations for you.</p>



