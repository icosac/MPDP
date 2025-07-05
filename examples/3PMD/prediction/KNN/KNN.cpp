#include <iostream>
#include <vector>
#include "annoylib.h"
#include "kissrandom.h"
#include <fstream>
#include <chrono>
#include <algorithm>

typedef float real_t;
using AnnoyIndex = Annoy::AnnoyIndex<int, real_t, Annoy::Euclidean, Annoy::Kiss32Random, Annoy::AnnoyIndexSingleThreadedBuildPolicy>;

#define RUN_TESTS 2

bool kill = false;

int rec_binary_search(const std::vector<std::tuple<int,int,int>> labels, int index, int low, int high) {
    int i = (int)(low + high) / 2;
    if (index >= std::get<1>(labels[i]) && index <= std::get<2>(labels[i])) {
        return std::get<0>(labels[i]);
    }
    else if (index <= std::get<1>(labels[i])) {
        return rec_binary_search(labels, index, low, i);
    }
    return rec_binary_search(labels, index, i+1, high);
}

int binary_search(const std::vector<std::tuple<int,int,int>> labels, int index){
    if (labels.size() == 0) {
        return -1;
    }
    else if (labels.size() == 1) {
        return std::get<0>(labels[0]);
    }
    else if (labels.size() == 2) {
        return index < std::get<1>(labels[1]) ? std::get<0>(labels[0]) : std::get<0>(labels[1]);
    }
    auto rec_res = rec_binary_search(labels, index, 0, labels.size()-1);
    return rec_res;
}


std::vector<int> predict_query(
    const AnnoyIndex& annoy_index,
    const std::vector<real_t>& query,
    const std::vector<std::tuple<int,int,int>>& labels,
    int k = 3
){
    if (k==3){
        std::cout << "Using default value for k: 3" << std::endl;
    }
    // Number of nearest neighbors
    std::vector<int> nearest_neighbors;
    std::vector<real_t> distances;

    // Get the k nearest neighbors
    annoy_index.get_nns_by_vector(query.data(), k, -1, &nearest_neighbors, &distances);

    std::vector<int> ret;
    for (int i = 0; i < nearest_neighbors.size(); i++) {
        ret.push_back(binary_search(labels, nearest_neighbors[i]));
    }

    return ret;
}


void test_set(
    const std::string& testset_name, 
    const AnnoyIndex& annoy_index, 
    const std::vector<std::tuple<int,int,int>>& labels,
    double & avg_time, double & min_time, double & max_time,
    int k = 3,
    int file_length = 0
){
    std::ifstream testset(testset_name);
    if (!testset.is_open()) {
        std::cerr << "Failed to open test set file!" << std::endl;
        return;
    }

    int n_features = 5;

    real_t kmax, theta_i, theta_f, alpha_m, alpha_f, th_m, len;
    int id_man_comb;

    // Skip first line of file
    std::string line;
    std::getline(testset, line);

    int num_correct = 0;
    int num_total = 0;

    while (testset >> kmax >> theta_i >> theta_f >> alpha_m >> alpha_f >> th_m >> id_man_comb >> len) {
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<real_t> query = {kmax, theta_i, theta_f, alpha_m, alpha_f};
        
        auto predicted_labels = predict_query(annoy_index, query, labels, k);
        
        if (predicted_labels.empty()) {
            continue;
        }

        if (std::find(predicted_labels.begin(), predicted_labels.end(), id_man_comb) != predicted_labels.end()) {
            num_correct ++;
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        double time = duration.count();
        avg_time += time;
        min_time = std::min(min_time, time);
        max_time = std::max(max_time, time);

        num_total ++;
        if (num_total % 10000 == 0) {
            if (file_length == 0) {
                std::cout << num_correct << "/" << num_total << std::endl;
            }
            else if (file_length != 0){
                printf("\r%.2f%% %.4f", num_total * 100.0 / file_length, avg_time / num_total);
                std::flush(std::cout);
            }
        }
    }

    std::cout << std::endl;
    std::cout << "Accuracy: " << 1.0 * num_correct / num_total << std::endl; 
    avg_time /= num_total;
}

int main() {
    // Define the number of dimensions (must match the Python model)
    int num_features = 5;
    AnnoyIndex annoy_index(num_features);

    // Load the saved Annoy index
    if (!annoy_index.load("knn.ann")) {
        std::cerr << "Failed to load Annoy index!" << std::endl;
        return 1;
    }

    std::cout << "Annoy index loaded successfully!" << std::endl;

    // Load labels intervals from file
    std::ifstream labels_file("y_labels_intervals.csv");
    if (!labels_file.is_open()) {
        std::cerr << "Failed to open labels file!" << std::endl;
        return 1;
    }

    std::vector<std::tuple<int, int, int>> labels;

    int label, lb, ub;
    while (labels_file >> label >> lb >> ub) {
        labels.push_back(std::make_tuple(label, lb, ub));
    }
    
    // Example query vector (same dimensionality as training data)
    std::vector<real_t> query = {0.5, 0.2, 0.8, 0.3, 0.9};

    // Perform KNN search
    int k = 21;
    auto predicted_labels = predict_query(annoy_index, query, labels, k);
    std::cout << "Predicted indices: ";
    for (int i = 0; i < predicted_labels.size(); i++) {
        std::cout << predicted_labels[i] << " ";
    }
    std::cout << std::endl;

    double avg_time = 0, min_time = 1e9, max_time = 0;

    switch(RUN_TESTS){
        case 1:
        case 3:
            avg_time = 0; min_time = 1e9; max_time = 0;
            test_set("/Users/enrico/Projects/mpdp/examples/3PMD/prediction/datasets/testset1.csv", annoy_index, labels, avg_time, min_time, max_time, k, 63504001);
            if (RUN_TESTS == 1) break;
        case 2:
            avg_time = 0; min_time = 1e9; max_time = 0;
            test_set("/Users/enrico/Projects/mpdp/examples/3PMD/prediction/datasets/testset2.csv", annoy_index, labels, avg_time, min_time, max_time, k, 2528173);
        default:
            break;
    }
    return 0;
}
