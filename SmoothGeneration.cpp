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

class Triangulation{

};

class Smooth_3Polytope{

};


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
	cout << "Hello World!";
	return 0;	
}