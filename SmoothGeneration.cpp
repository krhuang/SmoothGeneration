#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <time.h>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <map>
#include <chrono>
#include <cassert>
#include "helper_matrix_functions.h" //Some name-wise self-explanatory functions, for printing, subtracting, multiplying, adding, matrices and vectors and dictionaries
using namespace std;

double affine_transformation_time;
double dictionary_merge_check_time;
double edge_length_allocation_time;
//The maximum # of lattice points in the 3-polytopes we generate. Previous work of Lundman has gone up to 16
const int MAX_LATTICE_POINTS = 24;
int polytopes_produced;
int affine_transformations_done;
int triangulations; 

//Outputs the produced polytopes to a file
void output_polytope(int number_vertices, map<set<int>, vector<int>>){
	cout << "outputting polytope" << endl;
}

//Returns all possible partitions of #balls into #boxes, with possibility of not using all the boxes
void balls_and_boxes_helper(int balls, int boxes, vector<int>& current, vector<pair<vector<int>, int>>& result, int used_weight){
	if ( (int) current.size() == boxes){
		result.push_back(make_pair(current, used_weight));
		return;
	}

	for (int i = 0; i <= balls; i++){
		current.push_back(i);
		balls_and_boxes_helper(balls - i, boxes, current, result, used_weight+i);
		current.pop_back();
	}
}

vector<pair<vector<int>, int>> balls_and_boxes(int balls, int boxes){
	vector<pair<vector<int>, int>> result;
	vector<int> current;
	balls_and_boxes_helper(balls, boxes, current, result, 0);
	return result;
}

//Given the vertex coordinates of a Smooth Polygon, computes its edge lengths in clockwise order, starting from the origin
vector<int> compute_edge_lengths(const vector<vector<int>>& vertex_coordinates){
	//Computes the clockwise lattice edge-lengths of Smooth Polygon, for its initialization. 
	//Example:
	//>print_vector({{0, 0}, {0, 3}, {2,0}})
	//>3 1 2 
	vector<int> edge_lengths;
	vector<int> difference;
	for(int i=0; i< (int) vertex_coordinates.size()-1; i++){
		difference = subtract_vector(vertex_coordinates[i+1], vertex_coordinates[i]);
		edge_lengths.push_back(gcd(difference[0], difference[1]));
	}
	difference = vertex_coordinates[vertex_coordinates.size()-1];
	edge_lengths.push_back(gcd(difference[0], difference[1]));
	return edge_lengths;
}

//TODO: This function could be without passing a copy, ie with & notation. Could be slightly faster, and better in general
//Takes the vertices of a Smooth Polygon with the first vector being the origin, and puts it in standard position (rays from the origin are along the x, y axes)
vector<vector<int>> standard_position(vector<vector<int>> vertex_coordinates){
	int x_length, y_length;
	int number_vertices = vertex_coordinates.size();
	vector<vector<int>> standard_position_matrix;
	y_length = gcd(vertex_coordinates[1][0], vertex_coordinates[1][1]);
	x_length = gcd(vertex_coordinates[number_vertices - 1][0], vertex_coordinates[number_vertices-1][1]);
	standard_position_matrix = {{vertex_coordinates[number_vertices-1][0] / x_length, vertex_coordinates[number_vertices-1][1] / x_length},{vertex_coordinates[1][0] / y_length, vertex_coordinates[1][1] / y_length}};
	standard_position_matrix = matrix_inverse(standard_position_matrix);
	for(int i = 0; i < number_vertices; i++){
		vertex_coordinates[i] = matrix_multiply(standard_position_matrix, vertex_coordinates[i]);
	}
	return vertex_coordinates;
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
		affine_transformations_done++;
		auto start_time = std::chrono::high_resolution_clock::now();
		vector <int> translation_vector = origin_destination;
		vector<vector<int>> new_vertices;
		int y_length = edge_lengths[0];
		int x_length = edge_lengths[edge_lengths.size() - 1];
		vector<int> first_col = divide_vector(subtract_vector(x_destination, translation_vector), x_length);
		vector<int> second_col = divide_vector(subtract_vector(y_destination, translation_vector), y_length);
		vector<vector<int>> lin_transform_matrix = {first_col, second_col}; 
		for(int i=0; i< (int) edge_lengths.size(); i++){
			new_vertices.push_back(add_vector(matrix_multiply(lin_transform_matrix, vertex_coordinates[i]), translation_vector));
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		 	// Calculate duration in seconds
    	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    	affine_transformation_time = affine_transformation_time + duration.count();

		return new_vertices; 


    	// Stop measuring time
    	
	}

	//Takes a Smooth_Polygon object, and rotates its embedding, so that a new vertex is the origin. 
	//TODO: this could be simpler after we fix the standard_position function
	void rotate_embedding(){
		//Rotate the edge_length vector
		//e.g. - 1 2 1 1 becomes 2 1 1 1
		rotate(edge_lengths.begin(), edge_lengths.begin()+1, edge_lengths.end()); 
		//Make the first vertex the new origin
		vector<int> difference_vector = vertex_coordinates[1];
		for(int i = 0; i < number_vertices; i++){
			vertex_coordinates[i] = subtract_vector(vertex_coordinates[i], difference_vector);
		}
		//Rotate it so the origin is the first vector
		rotate(vertex_coordinates.begin(), vertex_coordinates.begin()+1, vertex_coordinates.end());
		//Then put it into standard position
		vertex_coordinates = standard_position(vertex_coordinates); 
	}
};
//Database of Smooth polygons. 
//This should have multiple entries for different embeddings of the same smooth polygon with various vertices being the origin
//TODO it should also have mirrored smooth polygons
map<vector<int>, set<Smooth_Polygon>> Smooth_Polygon_DB;

//Triangulations, usually given by plantri in text form. 
//Each triangulation has [number_vertices] vertices and its edges are given clockwise
class Triangulation{
	public:
		int number_vertices{ 0 };
		int number_edges{ 0 };
		vector<vector<int>> adjacencies{}; //plantri-form
		//vector<vector<int>> edge_weights{}; //same as plantri-form but now with weights?
		int total_edge_weight{ 0 }; 
		//vector<vector<int>> adjacency_matrix{}; //0 1 symmetric matrix of incidences
		//vector<vector<int>> edge_weights_matrix{}; //initialized to have weight one on edges
		vector<int> shelling_order{};
		//Shelling Inverse, an inverse function to the shelling order which returns the index in which the input appears
		map<int, int> shelling_order_inverse{};
		int smooth_polytope_vertex_count;
		//Triangulation constructor
		Triangulation(int input_number_vertices, vector<vector<int>> input_adjacencies) 
			: number_vertices(input_number_vertices), number_edges(0), adjacencies(input_adjacencies)
		{
			//edge_weights = adjacencies;
			for(int vertex = 0; vertex < number_vertices; vertex++){
				for(int adjacency = 0; adjacency< (int) adjacencies[vertex].size(); adjacency++){
					//edge_weights[vertex][adjacency] = 0;
					number_edges++;
				}
			}
			number_edges = number_edges / 2;
			total_edge_weight = number_edges; 
			smooth_polytope_vertex_count = 2 + number_edges - number_vertices;
		}
	//Computes an arbitrary shelling order on the triangulation
	//Here a shelling requires that the first three vertices form a triangle, and that every new vertex thereafter must form a triangle with two of the previous vertices of the shelling
	/*
	vector<vector<int>> initial_edge_lengths(){
		vector<vector<int>> edge_weights = adjacencies; 
		for(int vertex = 0; vertex < number_vertices; vertex++){
			for(int adjacency = 0; adjacency< (int) adjacencies[vertex].size(); adjacency++){
				edge_weights[vertex][adjacency] = 0;	
			}
		}
		return edge_weights;
	}
	*/

	void compute_a_shelling(){
		shelling_order.push_back(0); 
		for(int i = 0; (int) shelling_order.size() < number_vertices; i++){
			int current_vertex = shelling_order[i];
			//int j = adjacencies[current_vertex].size()-1;
			for(int j = adjacencies[current_vertex].size()-1; j >= 0; j--){ //going clockwise through the adjacencies of the current vertex
				if(element_of_vector(adjacencies[current_vertex][j], shelling_order) == false){
					shelling_order.push_back(adjacencies[current_vertex][j]);
				}
			}
		}
	}

	void compute_shelling_inverse(){
		for(int i = 0; i < (int) shelling_order.size(); i++){
			shelling_order_inverse[shelling_order[i]] = i;
		}
	}

	//Rotates the various adjacency vectors so that the first and last entry are previously inside the shelling_order
	//Should only be run *after* computing a shelling order! 
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
			int current_vertex = shelling_order[rebuilt_shelling.size()];
			while((element_of_vector(adjacencies[current_vertex][0], rebuilt_shelling) && element_of_vector(adjacencies[current_vertex][adjacencies[current_vertex].size()-1], rebuilt_shelling)) == false){
				rotate(adjacencies[current_vertex].begin(), adjacencies[current_vertex].begin()+1, adjacencies[current_vertex].end()); 
			}
			rebuilt_shelling.push_back(current_vertex);
			current_vertex++;
		}
	}

	void invert_shelling(){
		vector<vector<int>> new_adjacencies;
		for(int shelling_num = 0; shelling_num < number_vertices; shelling_num++){
			vector<int> new_adjacency = adjacencies[shelling_order[shelling_num]];
			for(int i = 0; i < (int) new_adjacency.size(); i++){
				new_adjacency[i] = shelling_order_inverse[new_adjacency[i]];
			}
			new_adjacencies.push_back(new_adjacency);
		}
		adjacencies = new_adjacencies;
		for(int shelling_num = 0; shelling_num < number_vertices; shelling_num++){
			shelling_order[shelling_num] = shelling_num;
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
		cout << "It is (potentially) the dual graph of a smooth polytope with " << smooth_polytope_vertex_count << " vertices. \n"; 
		cout << "Its Adjacencies are given by" << "\n";
		print_matrix(adjacencies);
		//cout << "Its Edge Weights are given by" << "\n";
		//print_matrix(edge_weights);
		cout << "and it has total edge weight " << total_edge_weight << "\n";
		cout << "Its Shelling Order is ";
		print_vector(shelling_order);
	}

	//Returns possible edge-length allocations along with the total weight used, as the second of the pair
	vector <pair<vector<int>, int>> edge_length_allocations(int shelling_num, int remaining_weight, const vector<vector<int>>& edge_weights){
		//Fill in edge weights
		//Iterate through neighbors of the current vertex
		auto start_time = std::chrono::high_resolution_clock::now();
		vector<int> previous_weight;
		int num_previous_neighbors = 0; 
		int prev_neighbor;
		map<int, int> previous_edge_weight_map;
		for(int adjacency_index = 0; adjacency_index < (int) adjacencies[shelling_num].size(); adjacency_index++){
			prev_neighbor = adjacencies[shelling_num][adjacency_index];
			if(prev_neighbor < shelling_num){
				//Increment the counter
				num_previous_neighbors++;
				//Record the already-assigned edge-weight, and its index, into a dictionary
				//First find the edge weight
				for(int i = 0; i < (int) edge_weights[prev_neighbor].size(); i++){
					if(adjacencies[prev_neighbor][i] == shelling_num){
						previous_edge_weight_map[adjacency_index] = edge_weights[prev_neighbor][i];
						//cout << "Currently on shelling num " << shelling_num << endl;
						//cout << adjacency_index << " => " << edge_weights[prev_neighbor][i] << endl;
					}
				}
			}
		}
		
		vector <pair<vector<int>, int>> result;
		for (auto& [partition, used_weight]:balls_and_boxes(remaining_weight,adjacencies[shelling_num].size() - num_previous_neighbors)){
			//Increments weights by one (eliminates 0s)
			increment_one(partition);
			//Then insert previous weights which have already been decided
			//print_vector(partition);
			map<int,int>::iterator it;
			for(it = previous_edge_weight_map.begin(); it != previous_edge_weight_map.end(); it++){
				//cout << it->first << it->second << endl;
				partition.insert(partition.begin() + it->first, it->second);
				//print_vector(partition);
				//partition.push_back(1);
			}
			//print_vector(partition);
			result.push_back(make_pair(partition, used_weight));
		}
		//cout << "Found " << num_previous_neighbors << " previous neighbors" << endl;
		auto end_time = std::chrono::high_resolution_clock::now();

    	// Calculate duration in seconds
    	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    	edge_length_allocation_time += duration.count();
		return result;
	}

	void build_all_polytopes(){
		compute_a_shelling();
		rotate_adjacencies();
		compute_shelling_inverse();
		invert_shelling();
		//build_polytopes(); with many different edge lengths
		//vector<vector<int>> initialized_edge_lengths = initial_edge_lengths();
		//print();
		auto start_time = std::chrono::high_resolution_clock::now();
		build_polytopes_edge_weights_test({}, 0, MAX_LATTICE_POINTS - smooth_polytope_vertex_count, {});
		auto end_time = std::chrono::high_resolution_clock::now();
    	// Calculate duration in seconds
    	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    	//cout << "Time taken on building all polytopes: " << duration.count() << " seconds" << endl;
	}

	//A recursive function that builds 3-polytopes and appends them to the global variable
	//This should also have as input the edge-weights map which tells the recursion what edge weights have been used so far
	//As well as a remaining_weight counter to know when we've built too large of a polytope
	void build_polytopes_edge_weights_test(map<set<int>, vector<int>> vertex_coordinates, int shelling_num, int remaining_weight, vector<vector<int>> edge_weights){ 
		vector<vector<int>> new_vertices = {};
		map<set<int>, vector<int>> new_vertex_coordinates = {};
		vector<int> curent_vertex_edge_lengths;
	
		int used_weight = 0;
		if(remaining_weight < 0){
			cout << "Went below 0 remaining_weight!" << endl;
		}
		else if(shelling_num == number_vertices){
			cout << "Finished iterating through the triangulation" << endl;
			print_dictionary(vertex_coordinates);
			output_polytope(smooth_polytope_vertex_count, vertex_coordinates);
			polytopes_produced++;
		}
		else if(shelling_num == 0){
			for(auto& [current_vertex_edge_lengths,used_weight] : edge_length_allocations(shelling_num, remaining_weight, edge_weights)){
				for(auto& polygon:Smooth_Polygon_DB[current_vertex_edge_lengths]){
					assert(current_vertex_edge_lengths == polygon.edge_lengths);
					new_vertices = polygon.vertex_coordinates; 
					new_vertex_coordinates = vertex_coordinates; 
					vector<vector<int>> new_edge_weights = edge_weights;
					new_edge_weights.push_back(current_vertex_edge_lengths);
					for(auto& coordinate:new_vertices){
						coordinate.push_back(0); 
					}
					int end = adjacencies[0].size();
					int neighbor;
					int prev;
					for(neighbor = 0, prev = end - 1; neighbor < end; prev = neighbor, neighbor++){
						new_vertex_coordinates[{shelling_order[0], adjacencies[shelling_order[0]][neighbor], adjacencies[shelling_order[0]][prev]}] = new_vertices[neighbor];
					}
					//print_dictionary(new_vertex_coordinates);
					build_polytopes_edge_weights_test(new_vertex_coordinates, shelling_num+1, remaining_weight - used_weight, new_edge_weights); 
				}
			}
		}
		else if(shelling_num == 1){
			for(auto& [current_vertex_edge_lengths,used_weight]:edge_length_allocations(shelling_num, remaining_weight, edge_weights)){
				for(auto& polygon:Smooth_Polygon_DB[current_vertex_edge_lengths]){
					assert(current_vertex_edge_lengths == polygon.edge_lengths);
					new_vertices = polygon.vertex_coordinates;
					new_vertex_coordinates = vertex_coordinates;
					vector<vector<int>> new_edge_weights = edge_weights;
					new_edge_weights.push_back(current_vertex_edge_lengths);
					for(auto& coordinate:new_vertices){
						coordinate.insert(coordinate.begin(),0); 
					}
					int end = adjacencies[1].size();
					int neighbor;
					int prev;
					for(neighbor = 0, prev = end - 1; neighbor < end; prev = neighbor, neighbor++){
						new_vertex_coordinates[{shelling_order[1], adjacencies[shelling_order[1]][neighbor], adjacencies[shelling_order[1]][prev]}] = new_vertices[neighbor];
					}
					//print_dictionary(new_vertex_coordinates);
					build_polytopes_edge_weights_test(new_vertex_coordinates, shelling_num + 1, remaining_weight-used_weight, new_edge_weights);
				}
			}
		}
		else if(shelling_num == 2){
			for(auto& [current_vertex_edge_lengths,used_weight]:edge_length_allocations(shelling_num, remaining_weight-used_weight, edge_weights)){
				for(auto& polygon:Smooth_Polygon_DB[current_vertex_edge_lengths]){
					assert(current_vertex_edge_lengths == polygon.edge_lengths);
					new_vertex_coordinates = vertex_coordinates;
					vector<vector<int>> new_edge_weights = edge_weights;
					new_edge_weights.push_back(current_vertex_edge_lengths);
					int y_length = polygon.edge_lengths[0];
					int x_length = polygon.edge_lengths[polygon.edge_lengths.size() - 1];
					new_vertices = polygon.Affine_Transf({0, 0, 0}, {0, 0, x_length}, {y_length, 0, 0} );
					int end = adjacencies[2].size();
					int neighbor;
					int prev;
					for(neighbor = 0, prev = end-1; neighbor < end; prev = neighbor, neighbor++){
						new_vertex_coordinates[{shelling_order[2], adjacencies[shelling_order[2]][neighbor], adjacencies[shelling_order[2]][prev]}] = new_vertices[neighbor];
					}
					//print_dictionary(new_vertex_coordinates);
					//cout << "Finished with the first three faces" << "\n";
					build_polytopes_edge_weights_test(new_vertex_coordinates, shelling_num + 1, remaining_weight - used_weight, new_edge_weights);
				}
			}
		}
		else{
			for(auto& [current_vertex_edge_lengths,used_weight]:edge_length_allocations(shelling_num, remaining_weight, edge_weights)){
				for(auto& polygon:Smooth_Polygon_DB[current_vertex_edge_lengths]){
					assert(current_vertex_edge_lengths == polygon.edge_lengths);
					new_vertex_coordinates = {};
					vector<vector<int>> new_edge_weights = edge_weights;
					new_edge_weights.push_back(current_vertex_edge_lengths);
					map<set<int>, vector<int>> vertex_coordinates_copy = vertex_coordinates;
					int end = polygon.number_vertices - 1;
					vector<int> origin_destination = vertex_coordinates[{shelling_num, adjacencies[shelling_num][0], adjacencies[shelling_num][end]}];
					vector<int> x_destination = vertex_coordinates[{shelling_num, adjacencies[shelling_num][end], adjacencies[shelling_num][end-1]}];
					vector<int> y_destination = vertex_coordinates[{shelling_num, adjacencies[shelling_num][0], adjacencies[shelling_num][1]}];
					new_vertices = polygon.Affine_Transf(origin_destination, x_destination, y_destination);
					//Record the new_vertices into the dictionary
					for(int i = 0, prev = end; i < polygon.number_vertices; prev = i, i++){
						new_vertex_coordinates[{shelling_num, adjacencies[shelling_num][i], adjacencies[shelling_num][prev]}] = new_vertices[i];
					}
					//print_dictionary(new_vertex_coordinates);
					auto start_time = std::chrono::high_resolution_clock::now();
					if(mergable(new_vertex_coordinates, vertex_coordinates_copy)){
						new_vertex_coordinates.merge(vertex_coordinates_copy);
						// Stop measuring time
    					auto end_time = std::chrono::high_resolution_clock::now();
   							// Calculate duration in seconds
   						auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
						dictionary_merge_check_time += duration.count();
						build_polytopes_edge_weights_test(new_vertex_coordinates, shelling_num + 1, remaining_weight - used_weight, new_edge_weights);
					}
					else{
							// Stop measuring time
   						auto end_time = std::chrono::high_resolution_clock::now();
   							// Calculate duration in seconds
   						auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
						dictionary_merge_check_time += duration.count();
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

	//K_4.build_polytopes_edge_weights_test({}, 0);
}

void cubeexample(){
	Triangulation Octahedron(6, {{1, 3, 4, 2}, {2, 5, 3, 0}, {0, 4, 5, 1}, {1, 5, 4, 0}, {3, 5, 2, 0}, {4, 3, 1, 2}});
	Octahedron.build_all_polytopes();
}

/*void haaseexample(){
	Triangulation Octahedron(6, {{1, 3, 4, 2}, {2, 5, 3, 0}, {0, 4, 5, 1}, {1, 5, 4, 0}, {3, 5, 2, 0}, {4, 3, 1, 2}});
	Octahedron.compute_a_shelling();
	Octahedron.rotate_adjacencies();
	Octahedron.compute_shelling_inverse();
	Octahedron.edge_weights[0] = {2, 2, 2, 2};
	Octahedron.edge_weights[1][3] = 2;
	Octahedron.edge_weights[2][3] = 2;
	Octahedron.edge_weights[3][0] = 2;
	Octahedron.edge_weights[4][0] = 2;
	Octahedron.print();
	Octahedron.build_polytopes({}, 0);
}	*/

void read_plantri_triangulation(string input_file_name){
	cout << "Reading Plantri PLANAR CODE-format planar triangulations from " << input_file_name << "..." << "\n";
	//Check that the first 15 characters is >>planar_code<<
	ifstream plantri_in(input_file_name);
	char input_char;
	string planar_code_header;
	for(int i = 0; i < 15; i++){
		plantri_in >> input_char;
		planar_code_header += input_char;
	}
	assert(planar_code_header == ">>planar_code<<");


	int number_vertices;
	//Reads the next integer as a char, then casts it into an integer. I tried to read it as a hex but this didn't work, and I'm not sure why
	while(plantri_in.peek() != EOF){
		int current_vertex = 0;
		input_char = plantri_in.get();
		number_vertices = static_cast<int>(input_char); 
		//Read a PLANTRI CODE-style triangulation of the file, putting the adjacencies into "adjacencies"
		//Warning! PLANTRI CODE indexes from 1, see the below example from the appendix 
		//5  3 4 0  3 5 0  1 4 5 2 0  1 5 3 0  2 3 4 0
		vector<vector<int>> adjacencies(number_vertices);
		int neighbor;
		while(current_vertex < number_vertices){
			input_char = plantri_in.get();
			neighbor = static_cast<int>(input_char);
			if(neighbor == 0){
				current_vertex++; //Hitting a 0 means it's the end of the adjacencies of that vertex
			}
			else{
				adjacencies[current_vertex].push_back(neighbor - 1); //Add the number to the adjacencies, but translating conventions (PLANAR CODE indexes from 1)
			}
		}
		Triangulation new_Triangulation = Triangulation(number_vertices, adjacencies);
		new_Triangulation.build_all_polytopes();
		triangulations++;
		cout << triangulations << endl;
	}
}

//Reads the smooth polygons from text file with filename "input_file_name"
//The input is assumed to be of the form
//	3: 0 0 1 0 0 1
//in clockwise order, and with #vertices as the start of the line
//This function reads the polygons and puts them into standard form
void read_polygon_DB(string input_file_name="Smooth2Polytopesfixed2.txt"){
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
		vertex_coordinates = standard_position(vertex_coordinates);

		//Everything below could be consolidated as a function? 
		//Compute its edge_lengths
		vector<int> edge_lengths = compute_edge_lengths(vertex_coordinates);
		//Define the incoming new Smooth Polygon
		Smooth_Polygon new_Polygon = Smooth_Polygon(number_vertices, 0, edge_lengths, vertex_coordinates);
		Smooth_Polygon_DB[new_Polygon.edge_lengths].insert(new_Polygon);
		//Now rotate its embeddings
		for(int i = 1; i < number_vertices; i++){
			new_Polygon.rotate_embedding();
			Smooth_Polygon_DB[new_Polygon.edge_lengths].insert(new_Polygon);
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
					-And rotate adjacencies so the first and last are previous numbers in the shelling order
					-Relabel everything so that the shelling order is 0, 1, 2, ...
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
	string input_plantri_file;
	input_plantri_file = "plantri_output14"; 
    
	read_polygon_DB();
	//cout << Smooth_Polygon_DB.size() << " Smooth Polygons in the Database... \n";
	auto start_time = std::chrono::high_resolution_clock::now();
	read_plantri_triangulation(input_plantri_file); 
	//cubeexample();
	// Stop measuring time
    auto end_time = std::chrono::high_resolution_clock::now();

   	// Calculate duration in seconds
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
	cout << "Time taken: " << duration.count() << "\n";
	cout << "Time taken on affine transformations: " << affine_transformation_time << " seconds" << endl;
	cout << "Time taken on checking dictionary mergability: " << dictionary_merge_check_time << " seconds" << endl;
	cout << "Time spent on edge length allocations: " << edge_length_allocation_time << endl;
	cout << polytopes_produced << " smooth 3-polytopes were produced" << endl;
	cout << affine_transformations_done << " affine transformations" << endl;
}