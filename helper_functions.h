// helper_functions.h
#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <cassert>

using namespace std;
using MyMap = map<set<int>, vector<int>>;

bool element_of_vector(int num, const vector<int>& vect);
void print_vector(const vector<int>& input_vect);
void print_matrix(const vector<vector<int>>& input_matrix);
void print_dictionary(const map<set<int>, vector<int>>& dictionary);
vector<int> matrix_multiply(vector<vector<int>> matrix, vector<int> vect);
void increment_one(vector<int>& vect);
vector<vector<int>> matrix_inverse(vector<vector<int>> matrix);
vector<int> divide_vector(vector<int> input_vector, int mod_factor);
vector<int> add_vector(vector<int> first_summand_vector, vector<int> second_summand_vector);
vector<int> subtract_vector(vector<int> minuend_vector, vector<int> subtrahend_vector);
vector<vector<int>> flip_x_y_coordinates(const vector<vector<int>>& vertex_coordinates);
bool mergable(MyMap const& a, MyMap const& b);


#endif // HELPER_FUNCTIONS_H
