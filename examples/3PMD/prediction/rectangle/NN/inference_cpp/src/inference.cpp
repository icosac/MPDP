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

void test_single(std::string model_path, ModelType model_type, bool use_gpu){
    std::cout << "Loading ONNX model: " << model_path << std::endl;

    auto time_model = std::chrono::high_resolution_clock::now();
    // Load the model
    OnnxModel model(model_path, use_gpu);

    std::cout << "Loading the model took"
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     std::chrono::high_resolution_clock::now() - time_model).count()
              << " microseconds" << std::endl;
    
    // Input data - for a model with input shape [-1, 5]
    std::vector<float> input_data = {
        1.7f, 1.0f, 0.1f, 0.1f,
        (float)std::sin(5.0*M_PI/12.0), (float)std::cos(5.0*M_PI/12.0),
        (float)std::sin(-M_PI/3.0), (float)std::cos(-M_PI/3.0)
    };
    
    std::cout << "Input data: ";
    for (auto val : input_data) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    // Single-task (classification) model
    if (model_type == CLASSIFICATION) {
        std::cout << "Running classification..." << std::endl;
        auto now = std::chrono::high_resolution_clock::now();
        std::vector<float> output = {};
        std::vector<float> probabilities = {};
        for (size_t i = 0; i < 1000; ++i) {
            output = model.run(input_data);

            // Apply softmax to get probabilities
            probabilities = softmax(output);
        }

        std::cout << "Solved in "
                << std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - now).count()/1000.0
                << " microseconds" << std::endl;
        
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
        std::vector<float> output = {};
        for (size_t i = 0; i < 1000; ++i) {
            output = model.run(input_data);
        }

        std::cout << "Inference completed in " 
                << std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - now).count() / 1000.0
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

void test_dataset_regression(const std::string & model_path, const std::string & testset_path, int n_samples, bool use_gpu, bool skip_first_line = true){
    std::ifstream testset_file(testset_path);
    if (!testset_file.is_open()) {
        std::cerr << "Failed to open testset file: " << testset_path << std::endl;
        return;
    }
    std::cout << "Loading ONNX model: " << model_path << std::endl;
    // Load the model
    OnnxModel model(model_path, use_gpu);
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
    std::ofstream out_file("inference_pred.csv");
    while(testset_file >> kmax >> c >> xm >> ym >> th_i >> th_f >> th_m >> id_man >> len){
        // Prepare input data
        std::vector<float> input_data = {kmax, c, xm, ym, (float)std::sin(th_i), (float)std::cos(th_i), (float)std::sin(th_f), (float)std::cos(th_f)};

        auto start = std::chrono::high_resolution_clock::now();
        std::vector<float> output = {};
        // Single-task (classification) model
        size_t n_iterations = 10;
        for (size_t i = 0; i < n_iterations; i++){
            output = model.run(input_data);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        float pred_thm = atan2f(output[0], output[1]);
        if (pred_thm < 0) pred_thm += 2.0f * M_PI;
        float pred_thm_deg = pred_thm * 180.0f / M_PI;

        float local_error = std::abs(pred_thm - th_m);
        error += local_error;

        double time = duration.count() / (n_iterations*1.0);  // average time per inference
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

        out_file << kmax << " " << c << " " << xm << " " << ym << " "
                 << th_i << " " << th_f << " " << th_m << " " << time << std::endl;

        // std::cout << "Inference completed in " << duration.count()/1000.0 << " milliseconds " << local_error << std::endl;

        // std::cout << "Predicted class: " << predicted_class << " with confidence: "
        //             << std::fixed << std::setprecision(4) << confidence << std::endl;
    }

    out_file.close();

    std::cout << std::endl << "###########################\n" << std::endl;
    std::cout << "Mean Error: " << error / (float)n_total << std::endl;
    std::cout << "Average time: " << avg_time / (float)n_total << " microseconds" << std::endl;
    std::cout << "Minimum time: " << min_time << " microseconds" << std::endl;
    std::cout << "Maximum time: " << max_time << " microseconds" << std::endl;
    std::cout << "Total samples: " << n_total << std::endl;
    std::cout << std::endl << "###########################" << std::endl;
}


void test_dataset_classification(
    const std::string & model_path,
    const std::string & testset_path,
    bool use_gpu,
    bool skip_first_line = true
){
    std::ifstream testset_file(testset_path);
    if (!testset_file.is_open()) {
        std::cerr << "Failed to open testset file: " << testset_path << std::endl;
        return;
    }
    std::cout << "Loading ONNX model: " << model_path << std::endl;
    // Load the model
    OnnxModel model(model_path, use_gpu);
    std::cout << "Model loaded successfully!" << std::endl;
    std::cout << "Running inference on testset..." << std::endl;

    // Skip first line
    if (skip_first_line) {
        std::string line;
        std::getline(testset_file, line);
    }

    size_t n_total = 0;
    double avg_time = 0, min_time = 1e9, max_time = 0;

    // Read testset data
    float kmax, c, xm, ym, th_i, th_f, th_m, len;
    int id_man;

    while(testset_file >> kmax >> c >> xm >> ym >> th_i >> th_f >> th_m >> id_man >> len){
        // Prepare input data
        std::vector<float> input_data = {kmax, c, xm, ym, (float)std::sin(th_i), (float)std::cos(th_i), (float)std::sin(th_f), (float)std::cos(th_f)};

        std::vector<float> output = {};
        std::vector<float> probabilities = {};
        // Single-task (classification) model
        auto start = std::chrono::high_resolution_clock::now();
        size_t n_iterations = 1;
        for (size_t i = 0; i < n_iterations; i++){
            output = model.run(input_data);
            probabilities = softmax(output);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        double time = duration.count() / (n_iterations*1.0);  // average time per inference
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
    }

    std::cout << std::endl << "###########################\n" << std::endl;
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
            std::cerr << "Usage: " << argv[0] << " <path_to_model.onnx> [--regression,--classification] [--gpu] [path_to_testset]" << std::endl;
            std::cerr << "If no testset is provided, it will run a single test." << std::endl;
            return 1;
        }

        else {
            std::string model_path = argv[1];
            std::string testset_path = "";
            ModelType model_type = NONE;
            bool use_gpu = false;

            for (int i = 2; i < argc; i++) {
                std::string arg = argv[i];
                if (arg == "--regression") {
                    model_type = REGRESSION;
                } else if (arg == "--classification") {
                    model_type = CLASSIFICATION;
                } else if (arg == "--gpu") {
                    use_gpu = true;
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
                if (model_type == CLASSIFICATION) {
                    test_dataset_classification(model_path, testset_path, use_gpu);
                }
                else if (model_type == REGRESSION) {
                    for (int n_samples = 1; n_samples < 2; n_samples++) {
                        test_dataset_regression(model_path, testset_path, n_samples, use_gpu);
                    }
                }
            } else {
                test_single(model_path, model_type, use_gpu);
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
