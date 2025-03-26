#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <onnxruntime_cxx_api.h>

class OnnxModel {
private:
    Ort::Env env;
    Ort::SessionOptions session_options;
    Ort::Session session{nullptr};
    std::vector<std::vector<int64_t>> input_node_dims;
    std::vector<std::vector<int64_t>> output_node_dims;
    const char* input_name = "input";
    const char* output_name = "output";

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
        
        input_node_dims.resize(num_input_nodes);
        output_node_dims.resize(num_output_nodes);

        // Get input dimensions
        for (size_t i = 0; i < num_input_nodes; i++) {
            // Get input name
            auto input_name_allocated = session.GetInputNameAllocated(i, allocator);
            std::string name_str = input_name_allocated.get();
            
            // Check if input name matches our hardcoded name
            std::cout << "Model input " << i << " name: " << name_str << std::endl;
            
            // Get input dimensions
            auto type_info = session.GetInputTypeInfo(i);
            auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
            input_node_dims[i] = tensor_info.GetShape();
            
            // Print input information
            std::cout << "Input " << i << " : name=" << name_str << std::endl;
            std::cout << "  Dimensions: ";
            for (auto dim : input_node_dims[i]) {
                std::cout << dim << " ";
            }
            std::cout << std::endl;
        }

        // Get output dimensions
        for (size_t i = 0; i < num_output_nodes; i++) {
            // Get output name
            auto output_name_allocated = session.GetOutputNameAllocated(i, allocator);
            std::string name_str = output_name_allocated.get();
            
            // Check if output name matches our hardcoded name
            std::cout << "Model output " << i << " name: " << name_str << std::endl;
            
            // Get output dimensions
            auto type_info = session.GetOutputTypeInfo(i);
            auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
            output_node_dims[i] = tensor_info.GetShape();
            
            // Print output information
            std::cout << "Output " << i << " : name=" << name_str << std::endl;
            std::cout << "  Dimensions: ";
            for (auto dim : output_node_dims[i]) {
                std::cout << dim << " ";
            }
            std::cout << std::endl;
        }
    }

    // Run inference with float input data
    std::vector<float> run(const std::vector<float>& input_data) {
        // Create input tensor
        auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
        
        // Handle dynamic dimensions (where dim == -1)
        std::vector<int64_t> input_dims = input_node_dims[0];
        
        // For this model, we have [-1, 5]
        // Let's set the batch size (first dimension) based on input size
        if (input_dims[0] == -1) {
            // Calculate batch size based on the total input size and remaining dimensions
            int64_t elements_per_batch = 1;
            for (size_t i = 1; i < input_dims.size(); i++) {
                elements_per_batch *= input_dims[i];
            }
            input_dims[0] = input_data.size() / elements_per_batch;
            
            std::cout << "Setting dynamic batch size to: " << input_dims[0] << std::endl;
            
            if (input_data.size() % elements_per_batch != 0) {
                throw std::runtime_error("Input size is not a multiple of the expected feature size");
            }
        }
        
        // Calculate total input size (product of all dimensions)
        size_t input_size = 1;
        for (auto dim : input_dims) {
            input_size *= dim;
        }
        
        if (input_data.size() != input_size) {
            throw std::runtime_error("Input data size does not match model input dimensions");
        }
        
        // Create input tensor with resolved dimensions
        Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
            memory_info,
            const_cast<float*>(input_data.data()),
            input_data.size(),
            input_dims.data(),
            input_dims.size()
        );
        
        // Using hardcoded input and output names
        const char* input_names[] = {input_name};
        const char* output_names[] = {output_name};
        
        std::cout << "Using input name: '" << input_name << "'" << std::endl;
        std::cout << "Using output name: '" << output_name << "'" << std::endl;
        
        // Run inference
        auto output_tensors = session.Run(
            Ort::RunOptions{nullptr},
            input_names,
            &input_tensor,
            1,
            output_names,
            1
        );
        
        // Get output data
        float* output_data = output_tensors[0].GetTensorMutableData<float>();
        
        // Get actual output shape (might have dynamic dimensions resolved)
        auto output_tensor_info = output_tensors[0].GetTensorTypeAndShapeInfo();
        auto output_dims = output_tensor_info.GetShape();
        
        std::cout << "Output shape after inference: ";
        for (auto dim : output_dims) {
            std::cout << dim << " ";
        }
        std::cout << std::endl;
        
        // Calculate total output size based on actual output dimensions
        size_t output_size = 1;
        for (auto dim : output_dims) {
            output_size *= dim;
        }
        
        // Copy output data to a vector
        std::vector<float> result(output_data, output_data + output_size);
        return result;
    }
};

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
        // Creating a batch of 2 samples, each with 5 features
        std::vector<float> input_data = {
            1.0f, 2.3562f, 2.3562f, 2.3562f, 1.5708,
        };
        
        // Run inference
        std::cout << "Running inference..." << std::endl;
        std::vector<float> output = model.run(input_data);
        
        // Print results
        std::cout << "Inference results:" << std::endl;
        for (size_t i = 0; i < output.size(); ++i) {
            std::cout << "  " << i << ": " << output[i] << std::endl;
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