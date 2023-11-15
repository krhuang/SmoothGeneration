#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <time.h>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <map>
using namespace std;
//The maximum # of lattice points in the 3-polytopes we generate. Previous work of Lundman has gone up to 16
int MAX_LATTICE_POINTS = 8;

//Checks if a number is an element of a vector
bool element_of_vector(int num, const vector<int>& vect) {
    return find(vect.begin(), vect.end(), num) != vect.end();
}

//Prints out a vector
void print_vector(const vector<int>& input_vect) {
    for (const auto& element : input_vect) {
        cout << element << " ";
    }
    cout << "\n";
}

//Prints out a matrix 
void print_matrix(const vector<vector<int>>& input_matrix) {
    for (const auto& row : input_matrix) {
        for (const auto& element : row) {
            cout << element << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

//Prints out a dictionary
void print_dictionary(const map<set<int>, vector<int>>& dictionary) {
    cout << "A dictionary with " << "\n";
    cout << "keys => entries" << "\n";

    for (const auto& [key, value] : dictionary) {
        for (const auto& iter : key) {
            cout << iter << " ";
        }
        cout << "=> ";
        for (const auto& iter : value) {
            cout << iter;
        }
        cout << "\n";
    }
}

//Multiplies a vector by a matrix. This only handles if the vector is 2-dimensional and the matrix is either 2x2 or 2x3.
vector<int> matrix_multiply(vector<vector<int>> matrix, vector<int> vect){
	if(matrix[0].size()==3){
		return {{matrix[0][0]*vect[0] + matrix[1][0]*vect[1], matrix[0][1]*vect[0] + matrix[1][1]*vect[1], matrix[0][2]*vect[0] + matrix[1][2]*vect[1]}};
	}
	else{
		return {{matrix[0][0]*vect[0] + matrix[1][0]*vect[1], matrix[0][1]*vect[0] + matrix[1][1]*vect[1]}};
	}

}

//Inverts a 2x2 matrix which has determinant 1
vector<vector<int>> matrix_inverse(vector<vector<int>> matrix){
	return {{matrix[1][1], -matrix[0][1]},{-matrix[1][0], matrix[0][0]}};
}

//Divides a entries of a vector by a factor
vector<int> divide_vector(vector<int> input_vector, int mod_factor){
	for(int i = 0; i < input_vector.size(); i++){
		if(input_vector[i] % mod_factor != 0){
			cout << "Unclean division warning!" << "\n";
		}
		input_vector[i] = input_vector[i] / mod_factor;
	}
	return input_vector;
}

//Adds together two vectors of equal size
vector<int> add_vector(vector<int> first_summand_vector, vector<int> second_summand_vector){
	if(first_summand_vector.size() != second_summand_vector.size()){
		cout << "Adding vectors of different size!" << "\n";
	}
	vector<int> sum_vector; 
	for(int i=0; i<first_summand_vector.size(); i++){
		sum_vector.push_back(first_summand_vector[i] + second_summand_vector[i]);
	}
	return sum_vector;
}

//Subtracts two vectors of equal size
vector<int> subtract_vector(vector<int> minuend_vector, vector<int> subtrahend_vector){
	if(minuend_vector.size() != subtrahend_vector.size()){
		cout << "Subtracting vectors of different size!" << "\n";
	}
	vector<int> difference_vector; 
	for(int i=0; i<minuend_vector.size(); i++){
		difference_vector.push_back(minuend_vector[i] - subtrahend_vector[i]);
	}
	return difference_vector; 
}

//Given the vertex coordinates of a Smooth Polygon, computes its edge lengths in clockwise order, starting from the origin
vector<int> compute_edge_lengths(vector<vector<int>> vertex_coordinates){
	//Computes the clockwise lattice edge-lengths of Smooth Polygon, for its initialization. 
	//Example:
	//>print_vector({{0, 0}, {0, 3}, {2,0}})
	//>3 1 2 
	vector<int> edge_lengths;
	vector<int> difference;
	for(int i=0; i<vertex_coordinates.size()-1; i++){
		difference = subtract_vector(vertex_coordinates[i+1], vertex_coordinates[i]);
		edge_lengths.push_back(gcd(difference[0], difference[1]));
	}
	difference = vertex_coordinates[vertex_coordinates.size()-1];
	edge_lengths.push_back(gcd(difference[0], difference[1]));
	return edge_lengths;
}

//Checks if two maps have overlapping keys with differing entries. If so, returns false. 
bool mergable(map<set<int>, vector<int>> map1, map<set<int>, vector<int>> map2){
	map1.insert(map2.begin(), map2.end());
	map2.insert(map1.begin(), map1.end());
	if(map1==map2){
		return true;
	}
	cout << "Dictionaries not combinable: vertex assignment does not line up (not an error)" << "\n";
	return false;
}

//The class of Smooth Polygons
class Smooth_Polygon{
public:
	int number_vertices{ 0 };
	int number_interior_lattice_points{ 0 };
	vector<int> edge_lengths{  }; //Edge lengths are given clockwise from the 0 0 vertex and in lattice-length format. The first edge is the longest one. 
	vector<vector<int>> vertex_coordinates{}; //Vertex coordinates

	//Default constructor
	Smooth_Polygon(int init_number_vertices, int init_number_interior_lattice_points, vector<int> init_edge_lengths, vector<vector<int>> init_coordinates)
		: number_vertices(init_number_vertices), number_interior_lattice_points(init_number_interior_lattice_points), edge_lengths(init_edge_lengths), vertex_coordinates(init_coordinates)
	{}

	bool operator<(const Smooth_Polygon& other) const {

		return vertex_coordinates < other.vertex_coordinates;
	}

	void print(){
		cout << "A Smooth Polygon with " << number_vertices << " vertices and " << number_interior_lattice_points << " interior lattice points." << "\n";
		cout << "Its vertices are " << "\n";
		print_matrix(vertex_coordinates);
		cout << "and it has clockwise edge lengths ";
		for(int i=0; i < number_vertices; i++){
			if(i < number_vertices - 1){
				cout << edge_lengths[i] << ", "; 
			}
			else{
				cout << edge_lengths[i] << "." << "\n";
			}
		}
	}

	vector<vector<int>> Affine_Transf(vector<int> origin_destination, vector<int> x_destination, vector<int> y_destination) const {
		//returns the vertices of the smooth polygon as embedded according to assigning the origin to origin_destination, the (a, 0) vertex to x_destination, and the (0, b) vertex to the y_destination
		vector <int> translation_vector = origin_destination;
		vector<vector<int>> new_vertices;
		int y_length = edge_lengths[0];
		int x_length = edge_lengths[edge_lengths.size() - 1];
		vector<int> first_col = divide_vector(subtract_vector(x_destination, translation_vector), x_length);
		vector<int> second_col = divide_vector(subtract_vector(y_destination, translation_vector), y_length);
		vector<vector<int>> lin_transform_matrix = {first_col, second_col}; 
		for(int i=0; i<edge_lengths.size(); i++){
			new_vertices.push_back(add_vector(matrix_multiply(lin_transform_matrix, vertex_coordinates[i]), translation_vector));
		}
		return new_vertices; 
	}
};
//Database of Smooth polygons. This should have multiple entries for different embeddings of the same smooth polygon with various vertices being 0,0
set<Smooth_Polygon> Smooth_Polygon_DB;

//Triangulations, usually given by plantri in text form. 
//Each triangulation has [number_vertices] vertices and its edges are given clockwise
class Triangulation{
	public:
		int number_vertices{ 0 };
		int number_edges{ 0 };
		vector<vector<int>> adjacencies{}; //plantri-form
		vector<vector<int>> edge_weights{}; //same as plantri-form but now with weights?
		int total_edge_weight{ 0 }; 
		//vector<vector<int>> adjacency_matrix{}; //0 1 symmetric matrix of incidences
		//vector<vector<int>> edge_weights_matrix{}; //initialized to have weight one on edges
		vector<int> shelling_order;
		Triangulation(int input_number_vertices, vector<vector<int>> input_adjacencies) //Triangulation constructor
			: number_vertices(input_number_vertices), number_edges(0), adjacencies(input_adjacencies)
		{
			/*adjacency_matrix = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}; //TODO: Fix this!! 
			for(int vertex = 0; vertex < number_vertices; vertex++){
				for(int adjacency = 0; adjacency < adjacencies[vertex].size(); adjacency++){
					adjacency_matrix[vertex][adjacencies[vertex][adjacency]] = 1;
					number_edges++;
				}
			}
			number_edges = number_edges / 2; 
			edge_weights_matrix = adjacency_matrix;*/
			edge_weights = adjacencies;
			for(int vertex = 0; vertex < number_vertices; vertex++){
				for(int adjacency = 0; adjacency< edge_weights[vertex].size(); adjacency++){
					edge_weights[vertex][adjacency] = 1;
					number_edges++;
				}
			}
			number_edges = number_edges / 2;
			total_edge_weight = number_edges; 
		}
	//Computes an arbitrary shelling order on the triangulation
	//Here a shelling requires that the first three vertices form a triangle, and that every new vertex thereafter must form a triangle with two of the previous vertices of the shelling
	void compute_a_shelling(){
		shelling_order.push_back(0); 
		for(int i = 0; shelling_order.size() < number_vertices; i++){
			int current_vertex = shelling_order[i];
			int j = adjacencies[current_vertex].size()-1;
			for(int j = adjacencies[current_vertex].size()-1; j >= 0; j--){ //going clockwise through the adjacencies of the current vertex
				if(element_of_vector(adjacencies[current_vertex][j], shelling_order) == false){
					shelling_order.push_back(adjacencies[current_vertex][j]);
				}
			}
		}
	}

	//Rotates the various adjacency vectors so that the first and last entry are previously inside the shelling_order
	void rotate_adjacencies(){ 
		vector<int> rebuilt_shelling = {};
		//The first three vertices. The origin of the smooth polygon should always be mapped to the origin in 3-space. 
		rebuilt_shelling.push_back(shelling_order[0]);
		rebuilt_shelling.push_back(shelling_order[1]);
		rebuilt_shelling.push_back(shelling_order[2]);
		while(adjacencies[shelling_order[0]][0] != shelling_order[1]){
			rotate(adjacencies[shelling_order[0]].begin(), adjacencies[shelling_order[0]].begin() + 1, adjacencies[shelling_order[0]].end());
		}

		while(adjacencies[shelling_order[1]][0] != shelling_order[2]){
			rotate(adjacencies[shelling_order[1]].begin(), adjacencies[shelling_order[1]].begin() + 1, adjacencies[shelling_order[1]].end());
		}
		while(adjacencies[shelling_order[2]][0] != shelling_order[0]){
			rotate(adjacencies[shelling_order[2]].begin(), adjacencies[shelling_order[2]].begin() + 1, adjacencies[shelling_order[2]].end());
		}
		while(rebuilt_shelling != shelling_order){
			int current_vertex = rebuilt_shelling.size();
			while((element_of_vector(adjacencies[current_vertex][0], rebuilt_shelling) && element_of_vector(adjacencies[current_vertex][adjacencies[current_vertex].size()-1], rebuilt_shelling)) == false){
				rotate(adjacencies[current_vertex].begin(), adjacencies[current_vertex].begin()+1, adjacencies[current_vertex].end()); 
			}
			rebuilt_shelling.push_back(shelling_order[current_vertex]);
			current_vertex++;
		}
	}
	//Returns the degree of a particular vertex
	int vertex_degree(int vertex_number){ 
		return adjacencies[vertex_number].size();
	}

	bool has_degree_3_or_degree_4_vertex(){
		//TODO

		return true;
	}

	void print(){
		cout << "A Triangulation with " << number_vertices << " vertices and " << number_edges << " edges." << "\n";
		cout << "Its Adjacencies are given by" << "\n";
		print_matrix(adjacencies);
		cout << "Its Edge Weights are given by" << "\n";
		print_matrix(edge_weights);
		cout << "and it has total edge weight " << total_edge_weight << "\n";
		cout << "Its Shelling Order is ";
		print_vector(shelling_order);
	}
	void build_polytopes(map<set<int>, vector<int>> vertex_coordinates, int shelling_num){ 
		//a recursive function that builds polytopes and appends them to the global variable
		vector<vector<int>> new_vertices = {};
		map<set<int>, vector<int>> new_vertex_coordinates = {};
		if(shelling_num == number_vertices){
			cout << "Finished iterating through the triangulation" << "\n";
		}
		else if(shelling_num == 0){
			for(auto& polygon:Smooth_Polygon_DB){
				if(edge_weights[shelling_order[0]] == polygon.edge_lengths){
					new_vertices = polygon.vertex_coordinates; 
					new_vertex_coordinates = vertex_coordinates; 
					for(auto& coordinate:new_vertices){
						coordinate.push_back(0); 
					}
					int end = edge_weights[shelling_order[0]].size();
					int neighbor;
					int prev;
					for(neighbor = 0, prev = end - 1; neighbor < end; prev = neighbor, neighbor++){
						new_vertex_coordinates[{shelling_order[0], adjacencies[shelling_order[0]][neighbor], adjacencies[shelling_order[0]][prev]}] = new_vertices[neighbor];
					}
					print_dictionary(new_vertex_coordinates);
					build_polytopes(new_vertex_coordinates, shelling_num+1); 
				}
			}
		}
		else if(shelling_num == 1){
			for(auto& polygon:Smooth_Polygon_DB){
				if(edge_weights[shelling_order[1]] == polygon.edge_lengths){
					new_vertices = polygon.vertex_coordinates;
					new_vertex_coordinates = vertex_coordinates;
					for(auto& coordinate:new_vertices){
						coordinate.insert(coordinate.begin(),0); 
					}
					int end = edge_weights[shelling_order[1]].size();
					int neighbor;
					int prev;
					for(neighbor = 0, prev = end - 1; neighbor < end; prev = neighbor, neighbor++){
						new_vertex_coordinates[{shelling_order[1], adjacencies[shelling_order[1]][neighbor], adjacencies[shelling_order[1]][prev]}] = new_vertices[neighbor];
					}
					print_dictionary(new_vertex_coordinates);
					build_polytopes(new_vertex_coordinates, shelling_num + 1);
				}
			}

		}
		else if(shelling_num == 2){
			for(auto& polygon:Smooth_Polygon_DB){
				if(edge_weights[shelling_order[2]] == polygon.edge_lengths){
					new_vertex_coordinates = vertex_coordinates;
					int y_length = polygon.edge_lengths[0];
					int x_length = polygon.edge_lengths[polygon.edge_lengths.size() - 1];
					new_vertices = polygon.Affine_Transf({0, 0, 0}, {0, 0, x_length}, {y_length, 0, 0} );
					
					int end = edge_weights[shelling_order[2]].size();
					int neighbor;
					int prev;
					for(neighbor = 0, prev = end-1; neighbor < end; prev = neighbor, neighbor++){
						new_vertex_coordinates[{shelling_order[2], adjacencies[shelling_order[2]][neighbor], adjacencies[shelling_order[2]][prev]}] = new_vertices[neighbor];
					}
					print_dictionary(new_vertex_coordinates);
					cout << "Finished with the first three faces" << "\n";
					build_polytopes(new_vertex_coordinates, shelling_num + 1);
				}

			}
		}
		else{
			for(auto& polygon:Smooth_Polygon_DB){
				if(edge_weights[shelling_num] == polygon.edge_lengths){
					new_vertex_coordinates = {};
					int end = polygon.number_vertices - 1;
					vector<int> origin_destination = vertex_coordinates[{shelling_order[shelling_num], adjacencies[shelling_order[shelling_num]][0], adjacencies[shelling_order[shelling_num]][end]}];
					vector<int> x_destination = vertex_coordinates[{shelling_order[shelling_num], adjacencies[shelling_order[shelling_num]][end], adjacencies[shelling_order[shelling_num]][end-1]}];
					vector<int> y_destination = vertex_coordinates[{shelling_order[shelling_num], adjacencies[shelling_order[shelling_num]][0], adjacencies[shelling_order[shelling_num]][1]}];
					new_vertices = polygon.Affine_Transf(origin_destination, x_destination, y_destination);
					//print_matrix(new_vertices);
					for(int i = 0, prev = end; i < polygon.number_vertices; prev = i, i++){
						new_vertex_coordinates[{shelling_order[shelling_num], adjacencies[shelling_order[shelling_num]][i], adjacencies[shelling_order[shelling_num]][prev]}] = new_vertices[i];
					}
					//print_dictionary(new_vertex_coordinates);
					if(mergable(new_vertex_coordinates, vertex_coordinates)){
						new_vertex_coordinates.merge(vertex_coordinates);
						print_dictionary(new_vertex_coordinates);
						build_polytopes(new_vertex_coordinates, shelling_num + 1);
					}
				}
			}
		}
	}
};

void unimodular3simplexexample(){
	//Initialization
	Triangulation K_4(4, {{1, 3, 2}, {2, 3, 0}, {0, 3, 1}, {1, 2, 0}}); 
	K_4.compute_a_shelling();
	K_4.rotate_adjacencies();
	K_4.print();
	//Smooth3Polytope(K_4, {0, 1, 2, 3}); //only input ""fixed"" triangulations!!

	K_4.build_polytopes({}, 0);
}

void cubeexample(){
	Triangulation Octahedron(6, {{1, 3, 4, 2}, {2, 5, 3, 0}, {0, 4, 5, 1}, {1, 5, 4, 0}, {3, 5, 2, 0}, {4, 3, 1, 2}});
	Octahedron.compute_a_shelling();
	Octahedron.rotate_adjacencies();
	//Octahedron.shelling_order = {0, 1, 2, 3, 4, 5};
	Octahedron.print();
	Octahedron.build_polytopes({}, 0); 
}

void haaseexample(){
	Triangulation Octahedron(6, {{1, 3, 4, 2}, {2, 5, 3, 0}, {0, 4, 5, 1}, {1, 5, 4, 0}, {3, 5, 2, 0}, {4, 3, 1, 2}});
	Octahedron.compute_a_shelling();
	Octahedron.rotate_adjacencies();
	//Octahedron.shelling_order = {0, 1, 2, 3, 4, 5};
	Octahedron.print();
	Octahedron.edge_weights = {{2, 2, 2, 2}}; 
	Octahedron.print();
}

void read_polygon_DB(string input_file_name="Smooth2Polytopesfixed2.txt"){
	//Reads the smooth polygons from text file with filename "input_file_name"
	//The input is assumed to be of the form
	//	3: 0 0 1 0 0 1
	//in clockwise order, and with #vertices as the start of the line

	//This function reads the polygons and puts them into standard form
	cout << "Reading Smooth Polygon Database from " << input_file_name << "..." << "\n";
	ifstream fin(input_file_name);
	int number_vertices;
	while(fin >> number_vertices){
		//Read a line of the file, putting the coordinates into "vertex_coordinates"
		char colon;
		int x_coordinate, y_coordinate;
		fin >> colon;
		vector<vector<int>> vertex_coordinates;
		int translation_x, translation_y;
		vertex_coordinates.push_back({0,0});
		fin >> translation_x >> translation_y;
		for(int i = 1; i<number_vertices; i++){
			fin >> x_coordinate >> y_coordinate;

			vertex_coordinates.push_back({x_coordinate - translation_x, y_coordinate - translation_y});
		}
		//After finishing with a line, process it into the Smooth_Polygon_DB
		//First put it into standard position
		int x_length, y_length;
		vector<vector<int>> standard_position_matrix;
		y_length = gcd(vertex_coordinates[1][0], vertex_coordinates[1][1]);
		x_length = gcd(vertex_coordinates[number_vertices - 1][0], vertex_coordinates[number_vertices-1][1]);
		standard_position_matrix = {{vertex_coordinates[number_vertices-1][0] / x_length, vertex_coordinates[number_vertices-1][1] / x_length},{vertex_coordinates[1][0] / y_length, vertex_coordinates[1][1] / y_length}};
		standard_position_matrix = matrix_inverse(standard_position_matrix);
		for(int i = 0; i < number_vertices; i++){
			vertex_coordinates[i] = matrix_multiply(standard_position_matrix, vertex_coordinates[i]);
		}
		//Compute its edge_lengths
		vector<int> edge_lengths = compute_edge_lengths(vertex_coordinates);
		//Define the incoming new Smooth Polygon
		Smooth_Polygon new_Polygon = Smooth_Polygon(number_vertices, 0, edge_lengths, vertex_coordinates);
		Smooth_Polygon_DB.insert(new_Polygon);
		//Now rotate its embeddings
		for(int i = 1; i < number_vertices; i++){
			Smooth_Polygon_DB.insert(new_Polygon.rotate_embedding())
		}
	}
}

int main(){
	/*
		Algorithm Steps:
			Construct/Import a library of smooth 2-polytopes (ie from Balletti)
			Construct/Import a database of triangulations of the sphere (ie from plantri)
				Massage the Data
					-Duplicate rotated non-identical embeddings
					-Fix Shelling order so that the first three are in the correct order
					-and rotate adjacencies so the first and last are previous numbers in the shelling order
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

	clock_t tStart = clock();
	read_polygon_DB();
	/*Smooth_Polygon Simplex_2 = Smooth_Polygon(3, 0, {1, 1, 1}, {{0, 0}, {0, 1}, {1, 0}});
	Smooth_Polygon_DB.push_back(Simplex_2); //initializing the smooth polytope database
	unimodular3simplexexample(); //running the example for constructing a unimodular simplex

	Smooth_Polygon Square = Smooth_Polygon(4, 0, {1, 1, 1, 1}, {{0, 0}, {0, 1}, {1, 1}, {1, 0}});
	Smooth_Polygon_DB.push_back(Square);
	cubeexample();

	Smooth_Polygon Hirzebruch_2 = Smooth_Polygon(4, 0, {1, 1, 1, 2}, {{0, 0}, {0, 1}, {1, 1}, {2, 0}});
	Smooth_Polygon_DB.push_back(Hirzebruch_2);
	Hirzebruch_2 = Smooth_Polygon(4, 0, {2, 1, 1, 1}, {{0, 0}, {0, 2}, {1, 1}, {1, 0}});
	Smooth_Polygon_DB.push_back(Hirzebruch_2);
	Hirzebruch_2 = Smooth_Polygon(4, 0, {1, 2, 1, 1}, {{0, 0}, {0, 1}, {2, 1}, {1, 0}});
	Smooth_Polygon_DB.push_back(Hirzebruch_2);
	Hirzebruch_2 = Smooth_Polygon(4, 0, {1, 1, 1, 2}, {{0, 0}, {0, 1}, {1, 2}, {1, 0}});
	Smooth_Polygon_DB.push_back(Hirzebruch_2);

	Smooth_Polygon Dilated_Square = Smooth_Polygon(4, 1, {2, 2, 2, 2}, {{0, 0}, {0, 2}, {2, 2}, {2, 0}});
	Smooth_Polygon_DB.push_back(Dilated_Square);*/

	//haaseexample();
	cubeexample();
	cout << "Time taken: \n" << (double)(clock()-tStart)/CLOCKS_PER_SEC << "\n";

}