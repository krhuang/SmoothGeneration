#include <iostream>
#include <vector>
using namespace std;
int MAX_LATTICE_POINTS = 8;

class Smooth_Polygon{
	public:
		//constructors
		Smooth_Polygon(){ //default constructor
		}

		Smooth_Polygon(int init_number_vertices, int init_number_interior_lattice_points, vector<int> init_edge_lengths, vector<pair<int, int>> init_coordinates)
			: number_vertices(init_number_vertices), number_interior_lattice_points(init_number_interior_lattice_points), edge_lengths(init_edge_lengths), coordinates(init_coordinates)
		{}

		int number_vertices;
		int number_interior_lattice_points;
		vector<int> edge_lengths; //edge lengths are given clockwise from the 0 0 vertex and in lattice-length format. The first edge is the longest one. 
		vector<pair<int, int>> coordinates; //vertex coordinates
	void print(){
		cout << "A Smooth Polygon with " << number_vertices << " vertices and " << number_interior_lattice_points << " interior lattice points." << endl;
		cout << "Its vertices are " << endl;
		for(int i=0; i < number_vertices; i ++){
			cout << coordinates[i].first << ", " << coordinates[i].second << endl;
		}
		cout << "and it has clockwise edge lengths ";
		for(int i=0; i < number_vertices; i++){
			if(i < number_vertices - 1){
				cout << edge_lengths[i] << ", "; 
			}
			else{
				cout << edge_lengths[i] << "." << endl;
			}
		}
	}\

	vector<vector<int>> Embedding_Vertex_Coordinates(int first_vertex, vector<int> first_point, vector<int> second_point, vector<int> third_point){
		//returns the vertices of the smooth polygon as embedded according to assigning the vertex to the first_coordinate and the one after it to the second_coordinate
		//make sure to mod by polygon number_vertices to make cyclic behaviour
			int second_vertex = (first_vertex + 1)%number_vertices;
			int third_vertex = (second_vertex + 1)%number_vertices;
			int first_length = edge_lengths[first_vertex];
			int second_length = edge_lengths[second_vertex];
		return {{1, 1}};
	}
};
Smooth_Polygon Smooth_Polygon_Database[1];

//Triangulations, usually given by plantri in text form. 
//Each triangulation has [number_vertices] vertices and its edges are given clockwise
class Triangulation{
	public:
		int number_vertices;
		vector<vector<int>> adjacencies;
		vector<vector<int>> edge_weights;
	int vertex_degree(int vertex_number){
		return adjacencies[vertex_number].size();
	}
	void print(){
		cout << "A Triangulation with " << number_vertices << " vertices." << endl;
		cout << "Its Adjacency Graph is given by the 0-1s matrix" << endl;
		
	}
};

class SimpleGraph{
	public:
		int number_vertices;
		vector<vector<int>> adjacencies; //for each vertex, denotes the vertices to which it is adjacent? 
	void add_edge(int source, int target){
		//adds an edge to the Simple Graph
	}
};

//Smooth 3-Polytopes
class Smooth3Polytope{
	public:	
		int number_vertices;
		vector<vector<int>> coordinates;

	void construct_from_triangulation(Triangulation input_triangulation){

	}
};

void unimodular3simplexexample(){
	//Initialization
	Triangulation K_4; 
	K_4.number_vertices = 4;
	K_4.adjacencies = {{1, 3, 2}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}};
	/*for(int row = 0; row < K_4.number_vertices; row ++){
		for(int col = 0; col < K_4.number_vertices; col ++){
			if(row != col){
				K_4.adjacency_graph[row][col] = 1;
			}
			else{
				K_4.adjacency_graph[row][col] = 0;
			}
		}
	}*/
	//K_4.print();
	K_4.edge_weights = K_4.adjacencies;
	for(int i=0; i<K_4.number_vertices; i++){
		for(int j = 0; j < K_4.adjacencies[i].size(); j++){
			K_4.edge_weights[i][j] = 1;
		}
	}

	Smooth_Polygon Simplex_2 = Smooth_Polygon(3, 0, {1, 1, 1}, {{0, 0}, {0, 1}, {1, 0}});
	/*Simplex_2.number_vertices = 3;
	Simplex_2.number_interior_lattice_points = 0;
	Simplex_2.edge_lengths = {1, 1, 1}; 
	Simplex_2.coordinates = {{0, 0}, {0, 1}, {1, 0}}; */
	Smooth_Polygon_Database[0] = Simplex_2;
	Smooth_Polygon_Database[0].print();


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