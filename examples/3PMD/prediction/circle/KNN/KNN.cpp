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
    int k, 
    int knn_indexes 
){
    // Number of nearest neighbors
    std::vector<int> nearest_neighbors;
    std::vector<real_t> distances;

    bool enough = false;
    int k_size = 0;
    std::map<int, int> label_map;

    while (!enough) {
        // Get the k nearest neighbors
        annoy_index.get_nns_by_vector(query.data(), knn_indexes, -1, &nearest_neighbors, &distances);
        
        // Merge neighbors with their distances
        std::vector<std::pair<int, real_t>> neighbors_with_distances;
        for (size_t i = 0; i < nearest_neighbors.size(); i++) {
            neighbors_with_distances.push_back(std::make_pair(nearest_neighbors[i], distances[i]));
        }
        // Sort by distance
        std::sort(neighbors_with_distances.begin(), neighbors_with_distances.end(), [](const std::pair<int, real_t>& a, const std::pair<int, real_t>& b) {
            return a.second < b.second;
        });

        label_map.clear();
        for (int i = 0; i < neighbors_with_distances.size() && label_map.size() < k; i++) {
            int label = binary_search(labels, neighbors_with_distances[i].first);
            if (label != -1) {
                if (label_map.find(label) == label_map.end()) {
                    label_map[label] = 0;
                }
                label_map[label]++;
            }
        }

        k_size = label_map.size();
        if (k_size < k){
            knn_indexes *= 10;
        }
        else {
            enough = true;
        }
    }
    // std::cout << "Label map (top " << k_size << "):" << std::endl;
    // for (const auto& pair : label_map) {
    //     std::cout << "Label: " << pair.first << ", Frequency: " << pair.second << std::endl;
    // }
    if (label_map.size() == 0){
        throw std::runtime_error("No labels found for the nearest neighbors.");
    }
    // Sort the labels by frequency
    std::vector<std::pair<int, int>> sorted_labels(label_map.begin(), label_map.end());
    std::sort(sorted_labels.begin(), sorted_labels.end(), [](const std::pair<int,int>& a, const std::pair<int,int>& b) {
        return a.second > b.second;
    });

    std::vector<std::pair<int, int>> ret;
    for (auto& pair : sorted_labels) {
        ret.push_back(std::make_pair(pair.first, pair.second));
        if (ret.size() >= k) {
            break;
        }
    }

    if (ret.size() < k){
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
    int knn_index,
    size_t k,
    int file_length
){
    std::cout << "Testing on file: " << testset_name << " with knn_index=" << knn_index << " and k=" << k << std::endl;
    std::ifstream testset(testset_name);
    if (!testset.is_open()) {
        std::cerr << "Failed to open test set file!" << std::endl;
        return;
    }

    int n_features = 5;

    real_t kmax, theta_i, theta_f, alpha_m, alpha_f, th_m, len;
    real_t id_man_comb_f;

    // Skip first line of file
    std::string line;
    std::getline(testset, line);

    int num_correct = 0;
    int num_exec = 0;
    int num_total = 0;

    double avg_query_time = 0;
    size_t num_queries = 0;

    while (testset >> kmax >> theta_i >> theta_f >> alpha_m >> alpha_f >> th_m >> id_man_comb_f >> len) {
        int id_man_comb = static_cast<int>(id_man_comb_f);
        // std::cout << kmax << " " << theta_i << " " << theta_f << " " << alpha_m << " " << alpha_f << " " << th_m << " " << id_man_comb << " " << len << std::endl;
        if (true || (rand() % 10000 == 0)) {
            // num_exec ++;
            // auto start = std::chrono::high_resolution_clock::now();
            // std::vector<real_t> query = {kmax, theta_i, theta_f, alpha_m, alpha_f};

            // auto start_query_time = std::chrono::high_resolution_clock::now();
            // std::vector<std::pair<int,int>> predicted_labels = predict_query(annoy_index, query, labels, k, knn_index);
            // auto end_query_time = std::chrono::high_resolution_clock::now();

            // avg_query_time += std::chrono::duration_cast<std::chrono::microseconds>(end_query_time - start_query_time).count();
            // num_queries++;

            // if (predicted_labels.empty()) {
            //     continue;
            // }

            // for (auto& pair : predicted_labels){
            //     if (pair.first == id_man_comb){
            //         num_correct++;
            //         break;
            //     }
            // }

            // auto end = std::chrono::high_resolution_clock::now();
            // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

            // double time = duration.count();
            // avg_time += time;
            // min_time = std::min(min_time, time);
            // max_time = std::max(max_time, time);

        }
        num_total ++;
        if (num_total % 10000 == 0) {
            if (file_length == 0) {
                std::cout << num_correct << "/" << num_exec << std::endl;
            }
            else if (file_length != 0){
                printf("\r%.2f%% with avg time: %.4f", num_total * 100.0 / file_length, avg_time / num_exec);
                std::flush(std::cout);
            }
        }
    }

    std::cout << std::endl;
    std::cout << "The average query time for k=" << k << " knn_index=" << knn_index << " is: " << (avg_query_time / num_queries) << " microseconds over " << num_queries << " queries." << std::endl;
    std::cout << "Accuracy: " << 100.0 * num_correct / num_exec << std::endl; 
    std::cout << num_total << " " << file_length << std::endl;
    avg_time /= num_exec;
}

int main() {
    srand(0);
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

    float label, lb, ub;
    while (labels_file >> label >> lb >> ub) {
        labels.push_back(std::make_tuple((int)label, (int)lb, (int)ub));
    }
    std::cout << "Labels loaded successfully! Total labels: " << labels.size() << std::endl;
    for (const auto& tup : labels) {
        std::cout << "Label: " << std::get<0>(tup) << ", Interval: [" << std::get<1>(tup) << ", " << std::get<2>(tup) << "]" << std::endl;
    }
    
    // Example query vector (same dimensionality as training data)
    std::vector<real_t> query = {2.0615, -1.3258176636680323, -4.4674103172578254, -1.5707963267948966, -2.6516353273360647};

    // Perform KNN search
    auto start = std::chrono::high_resolution_clock::now();
    auto predicted_labels = predict_query(annoy_index, query, labels, 3, int(1e0));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time = duration.count();
    std::cout << "Time taken for KNN search: " << time << " microseconds" << std::endl;

    std::cout << "Predicted labels and frequencies:" << std::endl;
    for (const auto& pair : predicted_labels) {
        std::cout << "Label: " << pair.first << ", Frequency: " << pair.second << " " << pair.second*100.0/3 << std::endl;
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
            // for (size_t max_labels = 1; max_labels < 5; max_labels++){
            std::cout << "################" << std::endl;
            test_set("/Users/enrico/Projects/mpdp/circ_new_test_cc.csv", annoy_index, labels, avg_time, min_time, max_time, 1e0, 1, 1305244);
            avg_time = 0; min_time = 1e9; max_time = 0; not_enough = 0; 
            std::cout << "################" << std::endl;
            test_set("/Users/enrico/Projects/mpdp/circ_new_test_cc.csv", annoy_index, labels, avg_time, min_time, max_time, 1e2, 2, 1305244);
            avg_time = 0; min_time = 1e9; max_time = 0; not_enough = 0; 
            std::cout << "################" << std::endl;
            test_set("/Users/enrico/Projects/mpdp/circ_new_test_cc.csv", annoy_index, labels, avg_time, min_time, max_time, 1e2, 3, 1305244);
            avg_time = 0; min_time = 1e9; max_time = 0; not_enough = 0; 
            std::cout << "################" << std::endl;
            test_set("/Users/enrico/Projects/mpdp/circ_new_test_cc.csv", annoy_index, labels, avg_time, min_time, max_time, 1e3, 4, 1305244);
            avg_time = 0; min_time = 1e9; max_time = 0; not_enough = 0; 
            std::cout << "################" << std::endl;
            test_set("/Users/enrico/Projects/mpdp/circ_new_test_cc.csv", annoy_index, labels, avg_time, min_time, max_time, 1e5, 5, 1305244);
            avg_time = 0; min_time = 1e9; max_time = 0; not_enough = 0; 
            std::cout << "################" << std::endl;
            // }
            if (RUN_TESTS == 1) break;
        case 2:
            avg_time = 0; min_time = 1e9; max_time = 0;
            test_set("/Users/enrico/Projects/mpdp/examples/3PMD/prediction/datasets/testset2_new.csv", annoy_index, labels, avg_time, min_time, max_time, 1e5, 3, 2104189);
        default:
            break;
    }

    return 0;
}
