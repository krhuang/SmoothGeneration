#include <iostream>
#include <set>
#include <map>
#include <bits/stdc++.h>
using namespace std;
int MAX_LATTICE_POINTS = 8;
void print_matrix(vector<vector<int>> input_matrix){ //prints a matrix for me, to test functions. 
	for(int row=0; row<input_matrix.size(); row++){
		for(int col=0; col<input_matrix[row].size(); col++){
			cout << input_matrix[row][col];
		}
		cout << endl;
	}
}

void print_dictionary(map<set<int>, vector<int>> dictionary){
	map <set<int>, vector<int>> :: iterator iter;
	cout << "A dictionary with" << endl;
	cout << "keys & entries" << endl;
	for(iter = dictionary.begin(); iter != dictionary.end(); iter++){
		for(auto j : (*iter).first){
			cout << j << " ";
		}
		//cout << 5 << endl;
		for(int i =0; i < (*iter).second.size(); i++){
			//cout << 5 << endl;
			cout << (*iter).second[i];
		}
		cout << endl;
	}
}

vector<int> divide_vector(vector<int> input_vector, int mod_factor){
	for(int i = 0; i < input_vector.size(); i++){
		if(input_vector[i] % mod_factor != 0){
			cout << "Unclean division warning!" << endl;
		}
		input_vector[i] = input_vector[i] / mod_factor;
	}
	return input_vector;
}

vector<int> subtract_vector(vector<int> minuend_vector, vector<int> subtrahend_vector){
	if(minuend_vector.size() != subtrahend_vector.size()){
		cout << "Subtracting vectors of different size!" << endl;
	}
	vector<int> difference_vector; 
	for(int i=0; i<minuend_vector.size(); i++){
		difference_vector[i] = minuend_vector[i] - subtrahend_vector[i];
	}
	return difference_vector; 
}

class Smooth_Polygon{
	public:
		//constructors
		Smooth_Polygon(){ //default constructor
		}

		Smooth_Polygon(int init_number_vertices, int init_number_interior_lattice_points, vector<int> init_edge_lengths, vector<vector<int>> init_coordinates)
			: number_vertices(init_number_vertices), number_interior_lattice_points(init_number_interior_lattice_points), edge_lengths(init_edge_lengths), vertex_coordinates(init_coordinates)
		{}

		int number_vertices{ 0 };
		int number_interior_lattice_points{ 0 };
		vector<int> edge_lengths{ {} }; //edge lengths are given clockwise from the 0 0 vertex and in lattice-length format. The first edge is the longest one. 
		vector<vector<int>> vertex_coordinates{ {0, 0} }; //vertex coordinates
	void print(){
		cout << "A Smooth Polygon with " << number_vertices << " vertices and " << number_interior_lattice_points << " interior lattice points." << endl;
		cout << "Its vertices are " << endl;
		print_matrix(vertex_coordinates);
		cout << "and it has clockwise edge lengths ";
		for(int i=0; i < number_vertices; i++){
			if(i < number_vertices - 1){
				cout << edge_lengths[i] << ", "; 
			}
			else{
				cout << edge_lengths[i] << "." << endl;
			}
		}
	}

	vector<vector<int>> Embedding_Vertex_Coordinates(int first_vertex, vector<int> first_point, vector<int> second_point, vector<int> third_point){
		//returns the vertices of the smooth polygon as embedded according to assigning the vertex to the first_coordinate and the one after it to the second_coordinate
		//make sure to mod by polygon number_vertices to ensure cyclic behaviour
			int second_vertex = (first_vertex + 1)%number_vertices;
			int third_vertex = (second_vertex + 1)%number_vertices;
			int first_length = edge_lengths[first_vertex];
			int second_length = edge_lengths[second_vertex];


		return {{1, 1}}; //placeholder
	}
};
Smooth_Polygon Smooth_Polygon_Database[1];//Database of Smooth polygons. This should have multiple entries for different embeddings of the same smooth polygon with various vertices being 0,0

//Triangulations, usually given by plantri in text form. 
//Each triangulation has [number_vertices] vertices and its edges are given clockwise
class Triangulation{
	public:
		Triangulation(int input_number_vertices, vector<vector<int>> input_adjacencies) //Triangulation constructor
			: number_vertices(input_number_vertices), adjacencies(input_adjacencies), number_edges(0)
		{
			adjacency_matrix = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}; //TODO: Fix this!! 
			for(int vertex = 0; vertex < number_vertices; vertex++){
				for(int adjacency = 0; adjacency < adjacencies[vertex].size(); adjacency++){
					adjacency_matrix[vertex][adjacencies[vertex][adjacency]] = 1;
					number_edges++;
				}
			}
			number_edges = number_edges / 2; 
			edge_weights_matrix = adjacency_matrix;
			edge_weights = adjacencies;
			for(int vertex = 0; vertex < number_vertices; vertex++){
				for(int adjacency = 0; adjacency< edge_weights[vertex].size(); adjacency++){
					edge_weights[vertex][adjacency] = 1;
				}
			}
		}
		int number_vertices{ 0 };
		int number_edges{ 0 };
		vector<vector<int>> adjacencies{}; //plantri-form
		vector<vector<int>> edge_weights{}; //same as plantri-form but now with weights?
		vector<vector<int>> adjacency_matrix{}; //0 1 symmetric matrix of incidences
		vector<vector<int>> edge_weights_matrix{}; //initialized to have weight one on edges
	int vertex_degree(int vertex_number){ //returns the degree of a particular vertex
		return adjacencies[vertex_number].size();
	}
	void print(){
		cout << "A Triangulation with " << number_vertices << " vertices and " << number_edges << " edges." << endl;
		cout << "Its Adjacency Graph is given by the 0-1s matrix" << endl;
		print_matrix(adjacency_matrix);
		cout << "Its Edge Weights are given by" << endl;
		print_matrix(edge_weights);
		
	}
};

//Smooth 3-Polytopes
class Smooth3Polytope{
	public:	
		int current_vertices; //counts the number of vertices as the polytope is being built up
		int total_vertices; //computed via euler characteristic
		map<set<int>, vector<int>> vertex_coordinates;

	Smooth3Polytope(Triangulation triangulation, vector<int> shelling_order){ //Constructs Smooth 3 Polytope(s) based on a weighted input-triangulation and a shelling_order
		current_vertices = 0; 
		total_vertices = 2 - triangulation.number_vertices + triangulation.number_edges; //#faces of the triangulation, via euler characteristic
		vertex_coordinates[{shelling_order[0], shelling_order[1], shelling_order[2]}] = {0, 0, 0};
		vector<vector<int>> new_vertices; //new vertices to be added to dictionary

		//0th vertex in shelling order. This corresponds to a face on the xy-plane
		if(triangulation.edge_weights[shelling_order[0]] == Smooth_Polygon_Database[0].edge_lengths){
			new_vertices = Smooth_Polygon_Database[0].vertex_coordinates;
			for(int i = 0; i < Smooth_Polygon_Database[0].number_vertices; i++){
				new_vertices[i].push_back(0);
			}
			for(int i = 0; i < Smooth_Polygon_Database[0].number_vertices; i++){
				vertex_coordinates[{shelling_order[0], triangulation.adjacencies[shelling_order[0]][i], triangulation.adjacencies[shelling_order[0]][(i-1 + Smooth_Polygon_Database[0].number_vertices) % Smooth_Polygon_Database[0].number_vertices]}] = new_vertices[i];
				//cout << triangulation.adjacencies[shelling_order[0]][i] << endl;
			}
			print_dictionary(vertex_coordinates);
		}

		//1st vertex
		if(triangulation.edge_weights[shelling_order[1]] == Smooth_Polygon_Database[0].edge_lengths){
			new_vertices = Smooth_Polygon_Database[0].vertex_coordinates;
			for(int i = 0; i < Smooth_Polygon_Database[0].number_vertices; i++){
				new_vertices[i].insert(new_vertices[i].begin(), 0);
			}
			for(int i = 0; i < Smooth_Polygon_Database[0].number_vertices; i++){
				vertex_coordinates[{shelling_order[1], triangulation.adjacencies[shelling_order[1]][i], triangulation.adjacencies[shelling_order[1]][(i-1 + Smooth_Polygon_Database[0].number_vertices) % Smooth_Polygon_Database[0].number_vertices]}] = new_vertices[i];
			}
			print_dictionary(vertex_coordinates);
		}

		//2nd vertex
		if(triangulation.edge_weights[shelling_order[2]] == Smooth_Polygon_Database[0].edge_lengths){
			new_vertices = Smooth_Polygon_Database[0].vertex_coordinates; 
			for(int i = 0; i < Smooth_Polygon_Database[0].number_vertices; i++){
				new_vertices[i].insert(new_vertices[i].begin()+1, 0);
			}
			for(int i = 0; i < Smooth_Polygon_Database[0].number_vertices; i++){
				vertex_coordinates[{shelling_order[2], triangulation.adjacencies[shelling_order[2]][i], triangulation.adjacencies[shelling_order[2]][(i-1 + Smooth_Polygon_Database[0].number_vertices) % Smooth_Polygon_Database[0].number_vertices]}] = new_vertices[i];
			}
			print_dictionary(vertex_coordinates);
		}

		//rest of the vertices
		int shelling_num = 3;
		while(shelling_num < triangulation.number_vertices){
			if(triangulation.edge_weights[shelling_order[shelling_num]] == Smooth_Polygon_Database[0].edge_lengths){
				vector<vector<int>> lin_transform_matrix = {{0, 0, 0}, {0, 0, 0}};
				vector<int> translation_matrix = {0, 0, 0};
				int y_length = Smooth_Polygon_Database[0].edge_lengths[0];
				int x_length = Smooth_Polygon_Database[0].edge_lengths[Smooth_Polygon_Database[0].edge_lengths.size()-1];
				//x and y-lengths of the Smooth polygon we want to insert
				//vertex
				translation_matrix = vertex_coordinates[{shelling_num, triangulation.adjacencies[shelling_order[shelling_num]][0], triangulation.adjacencies[shelling_order[shelling_num]][triangulation.adjacencies[shelling_num].size()-1]}];
				lin_transform_matrix = 
					{{vertex_coordinates[{shelling_num, triangulation.adjacencies[shelling_order[shelling_num]][0]}] }};
			}

			shelling_num++;
		}
	}
};


void unimodular3simplexexample(){
	//Initialization
	Triangulation K_4(4, {{1, 3, 2}, {2, 3, 0}, {0, 3, 1}, {1, 2, 0}}); 
	Smooth3Polytope(K_4, {0, 1, 2, 3}); //only input ""fixed"" triangulations!!
}

vector<vector<int>> unimodular_matrix_inverse(vector<vector<int>> input_matrix){ // inverts unimodular 2x2 matrices using the closed-form formula. The determinant is 1 already
	return {{input_matrix[1][1], -input_matrix[0][1]},{-input_matrix[1][0], input_matrix[0][0]}}; 
}

int main(){
	/*
		Algorithm Steps:
			Construct/Import a library of smooth polytopes
			Construct/Import a database of triangulations of the sphere (ie from plantri)
				Massage the Data!!
					-Duplicate rotated non-isomorphic embeddings
					-Fix Shelling order so that the first three are in the correct order
					-and rotate adjacencies so the first three are in the previous 
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
	Smooth_Polygon Simplex_2 = Smooth_Polygon(3, 0, {1, 1, 1}, {{0, 0}, {0, 1}, {1, 0}});
	Smooth_Polygon_Database[0] = Simplex_2; //initializing the smooth polytope database
	//Smooth_Polygon_Database[0].print();
	unimodular3simplexexample(); //running the example for constructing a unimodular simplex
}