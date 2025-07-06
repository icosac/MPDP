#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include "hnswlib/hnswlib.h"

// #include "/Users/enrico/Projects/mpdp/examples/3PMD/prediction/KNN/HNSWlib/build/_deps/eigen-src/Eigen/Dense"
#include "Eigen/Dense"


// Function to shuffle and split dataset
void shuffleAndSplitEigen(  Eigen::MatrixXf& X, 
                            Eigen::VectorXi& y, 
                            float train_ratio, 
                            Eigen::MatrixXf& X_train, 
                            Eigen::VectorXi& y_train, 
                            Eigen::MatrixXf& X_test, 
                            Eigen::VectorXi& y_test) 
{
    int num_samples = X.rows();

    // Generate shuffled indices
    std::vector<int> indices(num_samples);
    for (int i = 0; i < num_samples; i++) indices[i] = i;

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(indices.begin(), indices.end(), g);

    // Compute train/test split index
    int train_size = static_cast<int>(train_ratio * num_samples);

    // Create train and test sets using the shuffled indices
    X_train.resize(train_size, X.cols());
    y_train.resize(train_size);
    X_test.resize(num_samples - train_size, X.cols());
    y_test.resize(num_samples - train_size);

    for (int i = 0; i < num_samples; i++) {
        if (i < train_size) {
            X_train.row(i) = X.row(indices[i]);
            y_train(i) = y(indices[i]);
        } 
        else {
            X_test.row(i - train_size) = X.row(indices[i]);
            y_test(i - train_size) = y(indices[i]);
        }
    }
}

int main() {
    std::ifstream X_file("/Users/enrico/Projects/mpdp/examples/3PMD/prediction/datasets/big_smaller.csv");

    if (!X_file.is_open()) {
        std::cerr << "Failed to open input file!" << std::endl;
        return 1;
    }

    int n_features = 5;

    double kmax, theta_i, theta_f, alpha_m, alpha_f, th_m, len;
    int id_man_comb;

    Eigen::MatrixXf X;
    Eigen::VectorXi y;

    // Skip first line of file
    std::string line;
    std::getline(X_file, line);

    while (X_file >> kmax >> theta_i >> theta_f >> alpha_m >> alpha_f >> th_m >> id_man_comb >> len) {
        X.conservativeResize(X.rows() + 1, n_features);
        X.row(X.rows() - 1) << kmax, theta_i, theta_f, alpha_m, alpha_f;
        y.conservativeResize(y.rows() + 1);
        y(y.rows() - 1) = id_man_comb;
    }

    // Split the dataset into training and testing sets
    float train_ratio = 0.8;
    Eigen::MatrixXf X_train, X_test;
    Eigen::VectorXi y_train, y_test;
    shuffleAndSplitEigen(X, y, train_ratio, X_train, y_train, X_test, y_test);

    std::cout << "Training set size: " << X_train.rows() << std::endl;
    std::cout << "Training Labels: " << y_train.rows() << std::endl;
    std::cout << "Testing set size: " << X_test.rows() << std::endl;
    std::cout << "Testing labels: " << y_test.rows() << std::endl;

    int num_elements = X_train.rows(); // Number of data points
    int M = 16;              // HNSW parameter: number of edges per node
    int ef_construction = 100; // HNSW parameter: controls accuracy/speed

    if (ef_construction >= num_elements) {
        std::cerr << "Error: ef_construction should be less than the number of data points" << std::endl;
        ef_construction = num_elements - 1;
    }

    // Initialize space and index
    // hnswlib::L2Space space(n_features);
    hnswlib::InnerProductSpace space(n_features);
    hnswlib::HierarchicalNSW<float> index(&space, num_elements, M, ef_construction);

    // Add data points (random example)
    for (int i = 0; i < num_elements; i++) {
        auto data = X_train.row(i);
        index.addPoint(data.data(), y_train(i)); // Index is used as the label
    }

    // Save the trained model
    index.saveIndex("hnsw_index.bin");
    std::cout << "HNSW model saved!" << std::endl;

    // Let's test the accuracy
    int k = 3; // Number of nearest neighbors
    int correct = 0;

    for (int i = 0; i < X_test.rows(); i++) {
        auto query = X_test.row(i);
        
        auto result = index.searchKnnCloserFirst(query.data(), k);

        for (int j = 0; j < k; j++) {
            if (result[j].second == y_test(i)) {
                correct++;
                break;
            }
        }

        // auto result = index.searchKnn(query.data(), k);

        // std::vector<int> predicted_labels(k);
        // for (int j = 0; j < k; j++) {
        //     auto res = result.top();
        //     predicted_labels[j] = res.second; // Index is used as the label
        //     result.pop();
        // }

        // // Get the most common label
        // std::map<int, int> label_counts;
        // for (int j = 0; j < k; j++) {
        //     label_counts[predicted_labels[j]]++;
        // }

        // int max_count = 0;
        // int predicted_label = -1;
        // for (const auto& pair : label_counts) {
        //     if (pair.second > max_count) {
        //         max_count = pair.second;
        //         predicted_label = pair.first;
        //     }
        // }

        // if (predicted_label == y_test(i)) {
        //     correct++;
        // }
    }

    std::cout << "Accuracy: " << 1.0*correct/X_test.rows() << std::endl;
    

    return 0;
}
