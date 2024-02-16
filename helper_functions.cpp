// helper_functions.cpp
#include "helper_functions.h"

bool element_of_vector(int num, const vector<int>& vect) {
    return find(vect.begin(), vect.end(), num) != vect.end();
}

void print_vector(const vector<int>& input_vect) {
    for (const auto& element : input_vect) {
        cout << element << " ";
    }
    cout << "\n";
}

void print_matrix(const vector<vector<int>>& input_matrix) {
    for (const auto& row : input_matrix) {
        for (const auto& element : row) {
            cout << element << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

void print_dictionary(const map<set<int>, vector<int>>& dictionary) {
    cout << "A dictionary with " << "\n";
    cout << "keys => entries" << "\n";

    for (const auto& [key, value] : dictionary) {
        for (const auto& iter : key) {
            cout << iter << " ";
        }
        cout << "=> ";
        for (const auto& iter : value) {
            cout << iter << " ";
        }
        cout << "\n";
    }
}

vector<int> matrix_multiply(vector<vector<int>> matrix, vector<int> vect) {
    if (matrix[0].size() == 3) {
        return {matrix[0][0] * vect[0] + matrix[1][0] * vect[1],
                matrix[0][1] * vect[0] + matrix[1][1] * vect[1],
                matrix[0][2] * vect[0] + matrix[1][2] * vect[1]};
    } else {
        return {matrix[0][0] * vect[0] + matrix[1][0] * vect[1],
                matrix[0][1] * vect[0] + matrix[1][1] * vect[1]};
    }
}

vector<vector<int>> matrix_inverse(vector<vector<int>> matrix) {
    return {{matrix[1][1], -matrix[0][1]}, {-matrix[1][0], matrix[0][0]}};
}

void increment_one(vector<int>& input_vector) {
    for(int i = 0; i < (int) input_vector.size(); i++){
        input_vector[i]++;
    }
}

vector<int> divide_vector(vector<int> input_vector, int mod_factor) {
    for (int i = 0; i < (int) input_vector.size(); i++) {
        //Checks for the possibility of a weird division occurring 
        assert(input_vector[i] == 0 || input_vector[i] % mod_factor == 0);
        input_vector[i] = input_vector[i] / mod_factor;
    }
    return input_vector;
}

vector<int> add_vector(vector<int> first_summand_vector, vector<int> second_summand_vector) {
    if (first_summand_vector.size() != second_summand_vector.size()) {
        cout << "Adding vectors of different size!" << "\n";
    }
    vector<int> sum_vector;
    for (int i = 0; i < (int) first_summand_vector.size(); i++) {
        sum_vector.push_back(first_summand_vector[i] + second_summand_vector[i]);
    }
    return sum_vector;
}

vector<int> subtract_vector(vector<int> minuend_vector, vector<int> subtrahend_vector) {
    if (minuend_vector.size() != subtrahend_vector.size()) {
        cout << "Subtracting vectors of different size!" << "\n";
    }
    vector<int> difference_vector;
    for (int i = 0; i < (int) minuend_vector.size(); i++) {
        difference_vector.push_back(minuend_vector[i] - subtrahend_vector[i]);
    }
    return difference_vector;
}

//This implementation in O(max(n, m)) is from this SE thread: https://stackoverflow.com/questions/77809818/checking-two-stdmap-with-some-shared-keys-if-overlapping-keys-map-to-the-same/77809920?noredirect=1#comment137173196_77809920
bool mergable(MyMap const& a, MyMap const& b) {
    auto aIt = a.begin();
    auto bIt = b.begin();

    while (aIt != a.end() && bIt != b.end()) {
        if (aIt->first < bIt->first) {
            // a's key is less than b's key -> move to a's next entry
            ++aIt;
        } else if (bIt->first < aIt->first) {
            // b's key is less than a's key -> move to b's next entry
            ++bIt;
        } else {
            // the keys are equal -> make sure the values are equal...
            if (aIt->second != bIt->second) {
                return false; // value mismatch -> not mergable
            }
            // ...and advance to next entry of a and b
            ++aIt;
            ++bIt;
        }
    }

    return true; // no mismatches found -> maps are mergable
}


vector<vector<int>> flip_x_y_coordinates(const vector<vector<int>>& vertex_coordinates){
    vector<vector<int>> result;
    for(int vector_index = 0; vector_index < (int) vertex_coordinates.size(); vector_index++){
        result.push_back({vertex_coordinates[vector_index][1], vertex_coordinates[vector_index][0]});
    }
    return result;
}