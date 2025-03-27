#include <iostream>
#include <vector>
#include "hnswlib/hnswlib.h"

int main() {
    int dim = 5;
    int k = 3;

    // Load the index
    hnswlib::L2Space space(dim);
    hnswlib::HierarchicalNSW<float> index(&space, "hnsw_index.bin");

    std::cout << "HNSW index loaded successfully!" << std::endl;

    // Query vector
    std::vector<float> query = {0.5, 0.2, 0.8, 0.3, 0.9};

    // Perform KNN search
    auto result = index.searchKnn(query.data(), k);

    std::cout << "Predicted indices: ";
    while (!result.empty()) {
        auto res = result.top();
        std::cout << res.first << " "; // Output nearest neighbor indices
        result.pop();
    }
    std::cout << std::endl;

    return 0;
}
