#include <iostream>
#include <vector>
using namespace std;

class Smooth_Polygon{
	public:
		int number_vertices;
		int number_interior_lattice_points;
		vector<int> edge_lengths; //edge lengths are given clockwise from the 0 0 vertex and in lattice-length format
		vector<pair<int, int>> coordinates; //vertex coordinates
};

//Triangulations, usually given by plantri in text form. 
class Triangulation{
	public:
		int number_vertices;
		vector<pair<int, int>> edges;  
};

//Smooth 3-Polytopes
class Smooth_3Polytope{

};

void unimodular3simplexexample(){
	Triangulation K_4; 
	K_4.number_vertices = 4;
	cout << K_4.number_vertices << endl;
	K_4.edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}; 
}

int main(){
	/*
		Algorithm Steps:
			Construct/Import a library of smooth polytopes
			Construct/Import a database of triangulations of the sphere (ie from plantri)
			For a specific triangulation
				Assign edge lengths to the triangulation
					[STOP CONDITION: Check if number of vertices + edge lenghts - #edges + minimum interior lattice points is >= N]
						minimum interior lattice points - for example, every pentagon has at least 1 interior lattice point
				Assign polygons to vertices of triangulation with corresponding lengths
					[STOP CONDITION: Check if number of vertices + edge lengths - #edges + #interior lattice opints is >= N]
				Try to close up a 3D Polytope using the triangulation
	*/
	/* 
		Example of the unimodular 3-simplex. It has vertices (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)
		Its triangulation is K_4
		Its smooth polygons are all the unimodular 2-simplex
	*/
	unimodular3simplexexample();
}