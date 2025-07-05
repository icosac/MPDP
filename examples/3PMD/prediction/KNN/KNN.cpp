#include <iostream>
#include <vector>
#include "annoylib.h"
#include "kissrandom.h"
#include <fstream>
#include <chrono>
#include <algorithm>
#include <chrono>
#include <map>

//! TODO Create table with variation of MAX_ELEMENT and accuracy

typedef float real_t;
using AnnoyIndex = Annoy::AnnoyIndex<int, real_t, Annoy::Euclidean, Annoy::Kiss32Random, Annoy::AnnoyIndexSingleThreadedBuildPolicy>;

#define RUN_TESTS 1

size_t not_enough = 0;
bool kill = false;
size_t max_number = 0;

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


std::vector<std::pair<int, int>> predict_query(
    const AnnoyIndex& annoy_index,
    const std::vector<real_t>& query,
    const std::vector<std::tuple<int,int,int>>& labels,
    int max_labels,
    int k 
){
    if (k==3){
        std::cout << "Using default value for k: 3" << std::endl;
    }
    // Number of nearest neighbors
    std::vector<int> nearest_neighbors;
    std::vector<real_t> distances;

    // Get the k nearest neighbors
    annoy_index.get_nns_by_vector(query.data(), k, -1, &nearest_neighbors, &distances);

    // Get map with unique labels and their frequencies
    std::map<int, int> label_map;
    for (int i = 0; i < nearest_neighbors.size(); i++) {
        int label = binary_search(labels, nearest_neighbors[i]);
        if (label != -1) {
            label_map[label]++;
        }
    }
    // Sort the labels by frequency
    std::vector<std::pair<int, int>> sorted_labels(label_map.begin(), label_map.end());
    std::sort(sorted_labels.begin(), sorted_labels.end(), [](const std::pair<int,int>& a, const std::pair<int,int>& b) {
        return a.second > b.second;
    });

    std::vector<std::pair<int, int>> ret;
    for (auto& pair : sorted_labels) {
        ret.push_back(std::make_pair(pair.first, pair.second));
        if (ret.size() >= max_labels) {
            break;
        }
    }

    if (ret.size() < max_labels){
        not_enough ++;
    }

    max_number = std::max((size_t)max_number, ret.size());

    return ret;
}


void test_set(
    const std::string& testset_name, 
    const AnnoyIndex& annoy_index, 
    const std::vector<std::tuple<int,int,int>>& labels,
    double & avg_time, double & min_time, double & max_time,
    int k,
    size_t max_labels,
    int file_length
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
    int num_exec = 0;
    int num_total = 0;

    while (testset >> kmax >> theta_i >> theta_f >> alpha_m >> alpha_f >> th_m >> id_man_comb >> len) {
        if (true || (rand() % 1000 == 0)) {
            num_exec ++;
            auto start = std::chrono::high_resolution_clock::now();
            std::vector<real_t> query = {kmax, theta_i, theta_f, alpha_m, alpha_f};
            
            std::vector<std::pair<int,int>> predicted_labels = predict_query(annoy_index, query, labels, max_labels, k);
            
            if (predicted_labels.empty()) {
                continue;
            }

            bool found = false;
            for (auto& pair : predicted_labels){
                if (pair.first == id_man_comb){
                    found = true;
                    num_correct++;
                    break;
                }
            }

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

            double time = duration.count();
            avg_time += time;
            min_time = std::min(min_time, time);
            max_time = std::max(max_time, time);

        }
        num_total ++;
        if (num_total % 10000 == 0) {
            if (file_length == 0) {
                std::cout << num_correct << "/" << num_exec << std::endl;
            }
            else if (file_length != 0){
                printf("\r%.2f%% %.4f", num_total * 100.0 / file_length, avg_time / num_exec);
                std::flush(std::cout);
            }
        }
    }

    std::cout << std::endl;
    std::cout << "Accuracy: " << 100.0 * num_correct / num_exec << std::endl; 
    avg_time /= num_exec;
}

int main() {
    srand(0);
    // Define the number of dimensions (must match the Python model)
    int num_features = 5;
    AnnoyIndex annoy_index(num_features);

    // Load the saved Annoy index
    if (!annoy_index.load("knn_new.ann")) {
        std::cerr << "Failed to load Annoy index!" << std::endl;
        return 1;
    }

    std::cout << "Annoy index loaded successfully!" << std::endl;

    // Load labels intervals from file
    std::ifstream labels_file("y_labels_intervals_new.csv");
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
    std::vector<real_t> query = {2.0615, -1.3258176636680323, -4.4674103172578254, -1.5707963267948966, -2.6516353273360647};

    // Perform KNN search
    int k = 10;
    auto start = std::chrono::high_resolution_clock::now();
    auto predicted_labels = predict_query(annoy_index, query, labels, k, 3);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time = duration.count();
    std::cout << "Time taken for KNN search: " << time << " microseconds" << std::endl;

    std::cout << "Predicted labels and frequencies:" << std::endl;
    for (const auto& pair : predicted_labels) {
        std::cout << "Label: " << pair.first << ", Frequency: " << pair.second << " " << pair.second*100.0/k << std::endl;
    }

    // std::cout << "Predicted indices: ";
    // for (int i = 0; i < predicted_labels.size(); i++) {
    //     std::cout << predicted_labels[i] << " ";
    // }
    // std::cout << std::endl;

    double avg_time = 0, min_time = 1e9, max_time = 0;

    switch(RUN_TESTS){
        case 1:
        case 3:
            avg_time = 0; min_time = 1e9; max_time = 0;
            // for (size_t max_labels = 1; max_labels < 5; max_labels++){
            test_set("/Users/enrico/Projects/mpdp/new_ds_out_4.csv", annoy_index, labels, avg_time, min_time, max_time, 20, 1, 1200000);
            std::cout << "The number of neighbors was not enough in " << not_enough << " cases" << std::endl;
            not_enough = 0;
            test_set("/Users/enrico/Projects/mpdp/new_ds_out_4.csv", annoy_index, labels, avg_time, min_time, max_time, 500, 2, 1200000);
            std::cout << "The number of neighbors was not enough in " << not_enough << " cases" << std::endl;
            not_enough = 0;
            test_set("/Users/enrico/Projects/mpdp/new_ds_out_4.csv", annoy_index, labels, avg_time, min_time, max_time, 800, 3, 1200000);
            std::cout << "The number of neighbors was not enough in " << not_enough << " cases" << std::endl;
            not_enough = 0;
            test_set("/Users/enrico/Projects/mpdp/new_ds_out_4.csv", annoy_index, labels, avg_time, min_time, max_time, 1000, 4, 1200000);
            std::cout << "The number of neighbors was not enough in " << not_enough << " cases" << std::endl;
            not_enough = 0;
            // }
            if (RUN_TESTS == 1) break;
        case 2:
            avg_time = 0; min_time = 1e9; max_time = 0;
            test_set("/Users/enrico/Projects/mpdp/examples/3PMD/prediction/datasets/testset2_new.csv", annoy_index, labels, avg_time, min_time, max_time, k, 3, 2104189);
        default:
            break;
    }
    return 0;
}
