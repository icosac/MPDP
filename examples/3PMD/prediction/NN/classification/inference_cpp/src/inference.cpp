#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <onnxruntime_cxx_api.h>

class OnnxModel {
private:
    Ort::Env env;
    Ort::SessionOptions session_options;
    Ort::Session session{nullptr};
    std::vector<const char*> input_names;
    std::vector<const char*> output_names;
    std::vector<std::vector<int64_t>> input_node_dims;
    std::vector<std::vector<int64_t>> output_node_dims;
    bool is_multitask_model = false;

public:
    OnnxModel(const std::string& model_path) : env(ORT_LOGGING_LEVEL_WARNING, "onnx-model") {
        // Set graph optimization level
        session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);
        
        // Create session
        session = Ort::Session(env, model_path.c_str(), session_options);

        // Get input and output information
        Ort::AllocatorWithDefaultOptions allocator;
        
        // Get number of inputs and outputs
        size_t num_input_nodes = session.GetInputCount();
        size_t num_output_nodes = session.GetOutputCount();
        
        std::cout << "Model has " << num_input_nodes << " inputs and " 
                  << num_output_nodes << " outputs" << std::endl;
        
        if (num_output_nodes > 1) {
            is_multitask_model = true;
            std::cout << "Detected multi-task model with " << num_output_nodes << " outputs" << std::endl;
        }
        
        input_names.resize(num_input_nodes);
        output_names.resize(num_output_nodes);
        input_node_dims.resize(num_input_nodes);
        output_node_dims.resize(num_output_nodes);

        // Get input information
        for (size_t i = 0; i < num_input_nodes; i++) {
            // Get input name
            auto input_name = session.GetInputNameAllocated(i, allocator);
            input_names[i] = strdup(input_name.get());  // We need to duplicate the string
            
            // Get input dimensions
            auto type_info = session.GetInputTypeInfo(i);
            auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
            input_node_dims[i] = tensor_info.GetShape();
            
            // Print input information
            std::cout << "Input " << i << " : name=" << input_names[i] << std::endl;
            std::cout << "  Dimensions: ";
            for (auto dim : input_node_dims[i]) {
                std::cout << dim << " ";
            }
            std::cout << std::endl;
        }

        // Get output information
        for (size_t i = 0; i < num_output_nodes; i++) {
            // Get output name
            auto output_name = session.GetOutputNameAllocated(i, allocator);
            output_names[i] = strdup(output_name.get());  // We need to duplicate the string
            
            // Get output dimensions
            auto type_info = session.GetOutputTypeInfo(i);
            auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
            output_node_dims[i] = tensor_info.GetShape();
            
            // Print output information
            std::cout << "Output " << i << " : name=" << output_names[i] << std::endl;
            std::cout << "  Dimensions: ";
            for (auto dim : output_node_dims[i]) {
                std::cout << dim << " ";
            }
            std::cout << std::endl;
        }
    }

    ~OnnxModel() {
        // Free allocated strings
        for (auto name : input_names) {
            free((void*)name);
        }
        for (auto name : output_names) {
            free((void*)name);
        }
    }

    // Run inference with float input data for classification (single-task) model
    std::vector<float> run(const std::vector<float>& input_data) {
        // Create input tensor
        auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
        
        // Handle batch dimension
        std::vector<int64_t> input_dims = input_node_dims[0];
        
        // Set batch size for dynamic dimensions
        if (input_dims[0] == -1) {
            // For a single sample with 5 features
            input_dims[0] = 1;
        }
        
        // Create input tensor
        Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
            memory_info,
            const_cast<float*>(input_data.data()),
            input_data.size(),
            input_dims.data(),
            input_dims.size()
        );
        
        // Run inference
        auto output_tensors = session.Run(
            Ort::RunOptions{nullptr},
            input_names.data(),
            &input_tensor,
            1,
            output_names.data(),
            output_names.size()
        );
        
        // Get output data from the first output (classification)
        float* output_data = output_tensors[0].GetTensorMutableData<float>();
        
        // Get actual output shape
        auto output_tensor_info = output_tensors[0].GetTensorTypeAndShapeInfo();
        auto output_dims = output_tensor_info.GetShape();
        
        std::cout << "Output shape: ";
        for (auto dim : output_dims) {
            std::cout << dim << " ";
        }
        std::cout << std::endl;
        
        // Calculate total output size
        size_t output_size = 1;
        for (auto dim : output_dims) {
            output_size *= dim > 0 ? dim : 1; // Skip negative dimensions
        }
        
        // Copy output data to a vector
        std::vector<float> result(output_data, output_data + output_size);
        return result;
    }

    // Run inference with float input data for multi-task model
    std::pair<std::vector<float>, std::vector<float>> run_multitask(const std::vector<float>& input_data) {
        if (!is_multitask_model) {
            throw std::runtime_error("This is not a multi-task model");
        }
        
        // Create input tensor
        auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
        
        // Handle batch dimension
        std::vector<int64_t> input_dims = input_node_dims[0];
        
        // Set batch size for dynamic dimensions
        if (input_dims[0] == -1) {
            // For a single sample with 5 features
            input_dims[0] = 1;
        }
        
        // Create input tensor
        Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
            memory_info,
            const_cast<float*>(input_data.data()),
            input_data.size(),
            input_dims.data(),
            input_dims.size()
        );
        
        // Run inference
        auto output_tensors = session.Run(
            Ort::RunOptions{nullptr},
            input_names.data(),
            &input_tensor,
            1,
            output_names.data(),
            output_names.size()
        );
        
        // Process classification output
        float* class_output_data = output_tensors[0].GetTensorMutableData<float>();
        auto class_tensor_info = output_tensors[0].GetTensorTypeAndShapeInfo();
        auto class_dims = class_tensor_info.GetShape();
        
        size_t class_size = 1;
        for (auto dim : class_dims) {
            class_size *= dim > 0 ? dim : 1;
        }
        
        std::vector<float> class_result(class_output_data, class_output_data + class_size);
        
        // Process regression output
        float* reg_output_data = output_tensors[1].GetTensorMutableData<float>();
        auto reg_tensor_info = output_tensors[1].GetTensorTypeAndShapeInfo();
        auto reg_dims = reg_tensor_info.GetShape();
        
        size_t reg_size = 1;
        for (auto dim : reg_dims) {
            reg_size *= dim > 0 ? dim : 1;
        }
        
        std::vector<float> reg_result(reg_output_data, reg_output_data + reg_size);
        
        return {class_result, reg_result};
    }
    
    bool isMultitaskModel() const {
        return is_multitask_model;
    }
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

int main(int argc, char* argv[]) {
    try {
        // Check if model path is provided
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " <path_to_model.onnx>" << std::endl;
            return 1;
        }
        
        std::string model_path = argv[1];
        std::cout << "Loading ONNX model: " << model_path << std::endl;
        
        // Load the model
        OnnxModel model(model_path);
        
        // Example input data - for a model with input shape [-1, 5]
        // Single sample with 5 features
        std::vector<float> input_data = {
            1.0f, 2.3562f, 2.3562f, 1.5708f, 2.3562f
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