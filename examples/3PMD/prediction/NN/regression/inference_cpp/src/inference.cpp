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
    OnnxRegressionModel model(model_path);
    
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
  

    // Single-task (regression) model
    std::cout << "Running regression inference..." << std::endl;
    std::vector<float> output = model.run(input_data);

    double pred_th = std::atan2(output[0], output[1]);
            
    std::cout << "Predicted: " << pred_th << " in radians, " << pred_th*180.0/M_PI << " degrees" << std::endl; 
}

void test_dataset(const std::string & model_path, const std::string & testset_path, bool skip_first_line = true){
    std::ifstream testset_file(testset_path);
    if (!testset_file.is_open()) {
        std::cerr << "Failed to open testset file: " << testset_path << std::endl;
        return;
    }
    std::cout << "Loading ONNX model: " << model_path << std::endl;
    // Load the model
    OnnxRegressionModel model(model_path);
    std::cout << "Model loaded successfully!" << std::endl;
    std::cout << "Running inference on testset..." << std::endl;
    
    // Skip first line
    if (skip_first_line) {
        std::string line;
        std::getline(testset_file, line);
    }

    double avg_error = 0.0, min_error = 1e9, max_error = 0.0;
    double avg_time = 0.0, min_time = 1e9, max_time = 0.0;
    size_t n_total = 0;

    // Read testset data
    float kmax, th_i, th_f, alpha_m, alpha_f, th_m, len;
    int id_man;
    while(testset_file >> kmax >> th_i >> th_f >> alpha_m >> alpha_f >> th_m >> id_man >> len){
        // Prepare input data
        std::vector<float> input_data = {kmax, th_i, th_f, alpha_m, alpha_f};

        auto start = std::chrono::high_resolution_clock::now();
        
        // Single-task (classification) model
        std::vector<float> output = model.run(input_data);
        double pred_angle = std::atan2(output[0], output[1]);
        
        double error = pred_angle - th_m;
        // if (error < 0) error += 2.0*M_PI;
        // else if (error > 2.0*M_PI) error -= 2.0*M_PI;

        if (error < -M_PI) error += 2.0*M_PI;
        else if (error > M_PI) error -= 2.0*M_PI;

        avg_error += error*error;
        max_error = std::max(max_error, error);
        min_error = std::min(min_error, error);

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double time = duration.count();
        avg_time += time;
        min_time = std::min(min_time, time);
        max_time = std::max(max_time, time);
        n_total++;

        if (n_total % static_cast<int>(1e5) == 0) {
            std::cout << "\rProcessed " << std::fixed << std::setprecision(2) << static_cast<float>(n_total)*100.0/4074835.0 << "%" << " ";
            std::cout << "Finishing in " << (avg_time / n_total) * (4074835.0 - n_total) / 1e6 << " seconds" << " ";
            std::cout << "Average time: " << avg_time / n_total << " microseconds";
            std::flush(std::cout);
        }
        
        // std::cout << "Predicted class: " << predicted_class << " with confidence: " 
        //             << std::fixed << std::setprecision(4) << confidence << std::endl;
    }

    std::cout << std::endl << "###########################\n" << std::endl;
    std::cout << "Average error: " << avg_error/static_cast<double>(n_total) << std::endl;
    std::cout << "Min error: " << min_error << std::endl; 
    std::cout << "Max error: " << max_error << std::endl; 
    std::cout << "Average time: " << avg_time / n_total << " microseconds" << std::endl;
    std::cout << "Minimum time: " << min_time << " microseconds" << std::endl;
    std::cout << "Maximum time: " << max_time << " microseconds" << std::endl;
    std::cout << "Total samples: " << n_total << std::endl;
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
            test_dataset(model_path, testset_path);
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