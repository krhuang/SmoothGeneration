// helper_functions.h
#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <numeric>
#include <chrono>
#include <cassert>
#include <algorithm>

using namespace std;
using MyMap = map<set<int>, vector<int>>;

//vector & matrix operations
bool element_of_vector(int num, const vector<int>& vect);
void print_vector(const vector<int>& input_vect);
void print_matrix(const vector<vector<int>>& input_matrix);
vector<int> matrix_multiply(vector<vector<int>> matrix, vector<int> vect);
void increment_one(vector<int>& vect);
vector<vector<int>> matrix_inverse(vector<vector<int>> matrix);
vector<int> divide_vector(vector<int> input_vector, int mod_factor);
vector<int> add_vector(vector<int> first_summand_vector, vector<int> second_summand_vector);
vector<int> subtract_vector(vector<int> minuend_vector, vector<int> subtrahend_vector);
vector<vector<int>> flip_x_y_coordinates(const vector<vector<int>>& vertex_coordinates);

//dictionary functions
bool mergable(MyMap const& a, MyMap const& b);
void print_dictionary(const map<set<int>, vector<int>>& dictionary);

//printing polytopes
void print_polytope(int number_vertices, map<set<int>, vector<int>> vertex_dictionary);

//weight assignment partitioning
void balls_and_boxes_helper(int balls, int boxes, vector<int>& current, vector<pair<vector<int>, int>>& result, int used_weight);
vector<pair<vector<int>, int>> balls_and_boxes(int balls, int boxes);

//polytope operations
vector<int> compute_edge_lengths(const vector<vector<int>>& vertex_coordinates);
vector<vector<int>> standard_position(vector<vector<int>> vertex_coordinates);

#endif // HELPER_FUNCTIONS_H
