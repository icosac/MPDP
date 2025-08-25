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


enum ModelType {
    NONE,
    REGRESSION,
    CLASSIFICATION
};


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

void test_single(std::string model_path, ModelType model_type){
    std::cout << "Loading ONNX model: " << model_path << std::endl;
        
    // Load the model
    OnnxModel model(model_path);
    
    // Example input data - for a model with input shape [-1, 5]
    // Single sample with 5 features
    std::vector<float> input_data = {
        2.0f, 1.0f, 0.25f, 0.75f, (float)std::sin(M_PI/2.0), (float)std::cos(M_PI/2.0), (float)std::sin(-M_PI/2.0), (float)std::cos(-M_PI/2.0),
    };
    
    std::cout << "Input data: ";
    for (auto val : input_data) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    // Run inference based on model type
    if (model.isMultitaskModel()) {
    //     // Multi-task model
        std::cout << "Running multi-task inference..." << std::endl;

    //     auto [class_output, reg_output] = model.run_multitask(input_data);
        
    //     // Process classification output
    //     std::vector<float> class_probs = softmax(class_output);
        
    //     std::cout << "Classification probabilities:" << std::endl;
    //     for (size_t i = 0; i < class_probs.size(); ++i) {
    //         std::cout << "  Class " << i << ": " << std::fixed << std::setprecision(6) << class_probs[i] << std::endl;
    //     }
        
    //     // Find most likely class
    //     int predicted_class = std::distance(class_probs.begin(), 
    //                                       std::max_element(class_probs.begin(), class_probs.end()));
    //     float confidence = class_probs[predicted_class];
        
    //     std::cout << "Predicted class: " << predicted_class << " with confidence: " 
    //               << std::fixed << std::setprecision(4) << confidence << std::endl;
        
    //     // Process regression output
    //     if (reg_output.size() >= 2) {
    //         float sin_val = reg_output[0];
    //         float cos_val = reg_output[1];
    //         float angle_deg = calculate_angle(sin_val, cos_val);
            
    //         std::cout << "Regression values - sin: " << sin_val << ", cos: " << cos_val << std::endl;
    //         std::cout << "Predicted angle: " << angle_deg << " degrees" << std::endl;
    //     }
    } 
    // else {

    // Single-task (classification) model
    if (model_type == CLASSIFICATION) {
        std::cout << "Running classification..." << std::endl;
        auto now = std::chrono::high_resolution_clock::now();
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
    else if (model_type == REGRESSION) {
        std::cout << "Running regression..." << std::endl;
        auto now = std::chrono::high_resolution_clock::now();
        std::vector<float> output = model.run(input_data);

        std::cout << "Inference completed in " 
                << std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - now).count() 
                << " microseconds" << std::endl;
        
        for (auto val : output) {
            std::cout << val << " ";
        }

        float angle = atan2f(output[0], output[1]);
        if (angle < 0) angle += 2.0f * M_PI;
        float angle_deg = angle * 180.0f / M_PI;

        std::cout << std::endl << "Angle: " << angle << " " << angle_deg << " degrees" << std::endl;
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

    size_t n_total = 0;
    double avg_time = 0, min_time = 1e9, max_time = 0, error = 0.0;

    // Read testset data
    float kmax, c, xm, ym, th_i, th_f, th_m, len;
    int id_man;
    while(testset_file >> kmax >> c >> xm >> ym >> th_i >> th_f >> th_m >> id_man >> len){
        // Prepare input data
        std::vector<float> input_data = {kmax, c, xm, ym, (float)std::sin(th_i), (float)std::cos(th_i), (float)std::sin(th_f), (float)std::cos(th_f)};

        auto start = std::chrono::high_resolution_clock::now();
        
        // Single-task (classification) model
        std::vector<float> output = model.run(input_data);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        float pred_thm = atan2f(output[0], output[1]);
        if (pred_thm < 0) pred_thm += 2.0f * M_PI;
        float pred_thm_deg = pred_thm * 180.0f / M_PI;

        float local_error = std::abs(pred_thm - th_m);
        error += local_error;

        double time = duration.count();
        avg_time += time;
        min_time = std::min(min_time, time);
        max_time = std::max(max_time, time);
        n_total++;

        if (n_total % static_cast<int>(1e3) == 0) {
            std::cout << "\rProcessed " << std::fixed << std::setprecision(2) << static_cast<float>(n_total)*100.0/178942.0 << "%" << " ";
            std::cout << "Finishing in " << (avg_time / (float)n_total) * (178942.0 - n_total) / 1e6 << " seconds" << " ";
            std::cout << "Average time: " << avg_time / (float)n_total << " microseconds";
            std::flush(std::cout);
        }

        // std::cout << "Inference completed in " << duration.count()/1000.0 << " milliseconds " << local_error << std::endl;

        // std::cout << "Predicted class: " << predicted_class << " with confidence: "
        //             << std::fixed << std::setprecision(4) << confidence << std::endl;
    }

    std::cout << std::endl << "###########################\n" << std::endl;
    std::cout << "Mean Error: " << error / (float)n_total << std::endl;
    std::cout << "Average time: " << avg_time / (float)n_total << " microseconds" << std::endl;
    std::cout << "Minimum time: " << min_time << " microseconds" << std::endl;
    std::cout << "Maximum time: " << max_time << " microseconds" << std::endl;
    std::cout << "Total samples: " << n_total << std::endl;
    std::cout << std::endl << "###########################" << std::endl;
}


int main(int argc, char* argv[]) {
    try {
        // Check if model path is provided
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " <path_to_model.onnx> [--regression,--classification] [path_to_testset]" << std::endl;
            std::cerr << "If not testset is provided, it will run a single test." << std::endl;
            return 1;
        }

        else {
            std::string model_path = argv[1];
            std::string testset_path = "";
            ModelType model_type = NONE;

            for (int i = 2; i < argc; i++) {
                std::string arg = argv[i];
                if (arg == "--regression") {
                    model_type = REGRESSION;
                } else if (arg == "--classification") {
                    model_type = CLASSIFICATION;
                } else {
                    testset_path = arg;
                }
            }

            if (model_type == REGRESSION && model_type == CLASSIFICATION) {
                std::cerr << "Error: Cannot specify both --regression and --classification." << std::endl;
                return 1;
            }
            else if (model_type == NONE) {
                std::cerr << "Error: Must specify either --regression or --classification." << std::endl;
                return 1;
            }

            if (testset_path != "") {
                for (int n_samples = 1; n_samples < 2; n_samples++) {
                    test_dataset(model_path, testset_path, n_samples);
                }
            } else {
                test_single(model_path, model_type);
            }

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