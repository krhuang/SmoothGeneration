// helper_matrix_functions.cpp
#include "helper_matrix_functions.h"

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
            cout << iter;
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
    for(int i = 0; i < input_vector.size(); i++){
        input_vector[i]++;
    }
}

vector<int> divide_vector(vector<int> input_vector, int mod_factor) {
    for (int i = 0; i < input_vector.size(); i++) {
        if (input_vector[i] % mod_factor != 0) {
            cout << "Unclean division warning!" << "\n";
        }
        input_vector[i] = input_vector[i] / mod_factor;
    }
    return input_vector;
}

vector<int> add_vector(vector<int> first_summand_vector, vector<int> second_summand_vector) {
    if (first_summand_vector.size() != second_summand_vector.size()) {
        cout << "Adding vectors of different size!" << "\n";
    }
    vector<int> sum_vector;
    for (int i = 0; i < first_summand_vector.size(); i++) {
        sum_vector.push_back(first_summand_vector[i] + second_summand_vector[i]);
    }
    return sum_vector;
}

vector<int> subtract_vector(vector<int> minuend_vector, vector<int> subtrahend_vector) {
    if (minuend_vector.size() != subtrahend_vector.size()) {
        cout << "Subtracting vectors of different size!" << "\n";
    }
    vector<int> difference_vector;
    for (int i = 0; i < minuend_vector.size(); i++) {
        difference_vector.push_back(minuend_vector[i] - subtrahend_vector[i]);
    }
    return difference_vector;
}

//TODO:fix this without the copy functionality
bool mergable(map<set<int>, vector<int>> map1, map<set<int>, vector<int>> map2) {
    // Start measuring time
    map1.insert(map2.begin(), map2.end());
    map2.insert(map1.begin(), map1.end());
    if (map1 == map2) {
        return true;
    }
    cout << "Dictionaries not combinable: vertex assignment does not line up (not an error)" << "\n";
    return false;
}
