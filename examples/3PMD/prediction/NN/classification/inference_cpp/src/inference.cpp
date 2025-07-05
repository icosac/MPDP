#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <onnxruntime_cxx_api.h>

#include "model.hpp"

#include <chrono>


// Softmax function for normalizing output
std::vector<float> softmax(const std::vector<float>& input) {
    std::vector<float> result(input.size());
    float max_val = *std::max_element(input.begin(), input.end());
    float sum = 0.0f;
    
    // Compute exp(x_i - max_val) for numerical stability
    for (size_t i = 0; i < input.size(); i++) {
        result[i] = std::exp(input[i] - max_val);
        sum += result[i];
    }
    
    // Normalize
    for (size_t i = 0; i < input.size(); i++) {
        result[i] /= sum;
    }
    
    return result;
}

// Calculate angle from sin and cos values
float calculate_angle(float sin_val, float cos_val) {
    float angle_rad = std::atan2(sin_val, cos_val);
    float angle_deg = angle_rad * 180.0f / M_PI;
    return angle_deg;
}

void test_single(std::string model_path){
    std::cout << "Loading ONNX model: " << model_path << std::endl;
        
    // Load the model
    OnnxModel model(model_path);
    
    // Example input data - for a model with input shape [-1, 5]
    // Single sample with 5 features
    std::vector<float> input_data = {
        -1.0f, 1.0f, 2.3562f, 2.3562f, 1.5708f, 2.3562f, 1.0f, 2.3562f, 2.3562f, 1.5708f, 2.3562f
    };
    
    std::cout << "Input data: ";
    for (auto val : input_data) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    // Run inference based on model type
    if (model.isMultitaskModel()) {
        // Multi-task model
        std::cout << "Running multi-task inference..." << std::endl;
        auto [class_output, reg_output] = model.run_multitask(input_data);
        
        // Process classification output
        std::vector<float> class_probs = softmax(class_output);
        
        std::cout << "Classification probabilities:" << std::endl;
        for (size_t i = 0; i < class_probs.size(); ++i) {
            std::cout << "  Class " << i << ": " << std::fixed << std::setprecision(6) << class_probs[i] << std::endl;
        }
        
        // Find most likely class
        int predicted_class = std::distance(class_probs.begin(), 
                                          std::max_element(class_probs.begin(), class_probs.end()));
        float confidence = class_probs[predicted_class];
        
        std::cout << "Predicted class: " << predicted_class << " with confidence: " 
                  << std::fixed << std::setprecision(4) << confidence << std::endl;
        
        // Process regression output
        if (reg_output.size() >= 2) {
            float sin_val = reg_output[0];
            float cos_val = reg_output[1];
            float angle_deg = calculate_angle(sin_val, cos_val);
            
            std::cout << "Regression values - sin: " << sin_val << ", cos: " << cos_val << std::endl;
            std::cout << "Predicted angle: " << angle_deg << " degrees" << std::endl;
        }
    } else {
        // Single-task (classification) model
        std::cout << "Running classification inference..." << std::endl;
        std::vector<float> output = model.run(input_data);
        
        // Apply softmax to get probabilities
        std::vector<float> probabilities = softmax(output);
        
        // Print results
        std::cout << "Classification probabilities:" << std::endl;
        for (size_t i = 0; i < probabilities.size(); ++i) {
            std::cout << "  Class " << i << ": " << std::fixed << std::setprecision(6) << probabilities[i] << std::endl;
        }
        
        // Find most likely class
        int predicted_class = std::distance(probabilities.begin(), 
                                          std::max_element(probabilities.begin(), probabilities.end()));
        float confidence = probabilities[predicted_class];
        
        std::cout << "Predicted class: " << predicted_class << " with confidence: " 
                  << std::fixed << std::setprecision(4) << confidence << std::endl;
    }
}

void test_dataset(const std::string & model_path, const std::string & testset_path, int n_samples, bool skip_first_line = true){
    std::ifstream testset_file(testset_path);
    if (!testset_file.is_open()) {
        std::cerr << "Failed to open testset file: " << testset_path << std::endl;
        return;
    }
    std::cout << "Loading ONNX model: " << model_path << std::endl;
    // Load the model
    OnnxModel model(model_path);
    std::cout << "Model loaded successfully!" << std::endl;
    std::cout << "Running inference on testset..." << std::endl;
    
    // Skip first line
    if (skip_first_line) {
        std::string line;
        std::getline(testset_file, line);
    }

    size_t n_correct = 0, n_total = 0;
    double avg_time = 0, min_time = 1e9, max_time = 0;

    // Read testset data
    float kmax, th_i, th_f, alpha_m, alpha_f, th_m, len;
    int id_man;
    while(testset_file >> kmax >> th_i >> th_f >> alpha_m >> alpha_f >> th_m >> id_man >> len){
        // Prepare input data
        std::vector<float> input_data = {kmax, th_i, th_f, alpha_m, alpha_f};

        auto start = std::chrono::high_resolution_clock::now();
        
        // Single-task (classification) model
        std::vector<float> output = model.run(input_data);
        
        // Apply softmax to get probabilities
        std::vector<float> probabilities = softmax(output);
        std::vector<std::pair<int, float>> sorted_probabilities;
        for (size_t i = 0; i < probabilities.size(); ++i) {
            sorted_probabilities.push_back({i+1, probabilities[i]});
        }
        std::sort(sorted_probabilities.begin(), sorted_probabilities.end(), [](const auto& a, const auto& b) {
            return a.second > b.second;
        });

        auto it = std::find_if(sorted_probabilities.begin(), sorted_probabilities.begin()+n_samples, 
            [&id_man](const std::pair<int, float>& el){ return el.first == id_man; });

        if (it != sorted_probabilities.begin()+n_samples) {
            n_correct++;
        }
        else {
            // if (n_samples == 3 || n_samples == 4 || n_samples == 5) {
            //     std::cout << n_total << " ";
            // }
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double time = duration.count();
        avg_time += time;
        min_time = std::min(min_time, time);
        max_time = std::max(max_time, time);
        n_total++;

        if (n_total % static_cast<int>(1e5) == 0) {
            std::cout << "\rProcessed " << std::fixed << std::setprecision(2) << static_cast<float>(n_total)*100.0/1200000.0 << "%" << " ";
            std::cout << "Finishing in " << (avg_time / n_total) * (1200000.0 - n_total) / 1e6 << " seconds" << " ";
            std::cout << "Average time: " << avg_time / n_total << " microseconds";
            std::flush(std::cout);
        }

        if (n_total == 1000000){
            break;
        }
        
        // std::cout << "Predicted class: " << predicted_class << " with confidence: " 
        //             << std::fixed << std::setprecision(4) << confidence << std::endl;
    }

    std::cout << std::endl << "###########################\n" << "n_samples: " << n_samples << std::endl;
    std::cout << "Accuracy: " << static_cast<float>(n_correct) / static_cast<float>(n_total) * 100.0f << "%" << std::endl;
    std::cout << "Average time: " << avg_time / n_total << " microseconds" << std::endl;
    std::cout << "Minimum time: " << min_time << " microseconds" << std::endl;
    std::cout << "Maximum time: " << max_time << " microseconds" << std::endl;
    std::cout << "Total samples: " << n_total << std::endl;
    std::cout << "Correct predictions: " << n_correct << std::endl;
    std::cout << "Incorrect predictions: " << n_total - n_correct << std::endl;
    std::cout << std::endl << "###########################" << std::endl;
}


int main(int argc, char* argv[]) {
    try {
        // Check if model path is provided
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " <path_to_model.onnx> [path_to_testset]" << std::endl;
            std::cerr << "If not testset is provided, it will run a single test." << std::endl;
            return 1;
        }

        else if (argc == 2) {
            std::string model_path = argv[1];
            test_single(model_path);
        }
        else if (argc == 3) {
            std::string model_path = argv[1];
            std::string testset_path = argv[2];
            for (int n_samples = 1; n_samples < 5; n_samples++) {
                // std::cout << "Running test with " << n_samples << " samples..." << std::endl;
                test_dataset(model_path, testset_path, n_samples);
            }
            // test_dataset(model_path, testset_path, 3);
            // test_dataset(model_path, testset_path, 4);
            // test_dataset(model_path, testset_path, 5);
        }        
        return 0;
    }
    catch (const Ort::Exception& e) {
        std::cerr << "ONNX Runtime error: " << e.what() << std::endl;
        return 1;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}