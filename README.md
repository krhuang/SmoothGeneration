<h2>Generating Smooth 3-Polytopes</h2>
<p>This repository contains the code and generated dataset of my master's thesis "Generating Smooth 3-Polytopes", which will be available online soon. Implementing and running a noval algorithm for generating smooth 3-polytopes took about half a day to extend the pre-existing database of 103 smooth 3-polytopes with <=16 lattice points (by Lundman) to 1364 smooth 3-polytopes with <=31 lattice points (available in the file prune_output). Oda's Conjecture was verified on all generated polytopes.</p>

<p>The theoretical background of the algorithm is explained in Section 3 of the thesis. The algorithm and its exact implementation is described in Sections 4 and 5. Detailed references are also contained therein.</p>

<h2>Where is the Dataset?</h2>
<p>The generated dataset of polytopes is available in the file `prune_output`.</p>

<p>Each smooth 3-polytope is given by a series of lines. The first line is the number of vertices, following by lines of triples of integers denoting the coordinates</p>

<h2>How to Run the Computation Yourself</h2>
The algorithm utilizes other mathematical computation projects, namely plantri, Balletti's database of smooth polygons, and the many packages included in polymake. 
