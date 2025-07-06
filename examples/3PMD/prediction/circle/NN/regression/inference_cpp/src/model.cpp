#include "model.hpp"

OnnxRegressionModel::OnnxRegressionModel(const std::string& model_path) : env(ORT_LOGGING_LEVEL_WARNING, "onnx-model") {
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

OnnxRegressionModel::~OnnxRegressionModel() {
    // Free allocated strings
    for (auto name : input_names) {
        free((void*)name);
    }
    for (auto name : output_names) {
        free((void*)name);
    }
}


std::vector<float> OnnxRegressionModel::run(const std::vector<float>& input_data) {
    // Create input tensor
    auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
    
    std::vector<int64_t> input_dims = input_node_dims[0];

    // Set batch size dynamically if required
    if (input_dims[0] == -1) {
        input_dims[0] = 1;
    }
    
    Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
        memory_info,
        const_cast<float*>(input_data.data()),
        input_data.size()*4.0,
        input_dims.data(),
        input_dims.size()
    );
    
    auto output_tensors = session.Run(
        Ort::RunOptions{nullptr},
        input_names.data(),
        &input_tensor,
        1,
        output_names.data(),
        output_names.size()
    );

    // Ensure the model returns exactly two values
    float* output_data = output_tensors[0].GetTensorMutableData<float>();
    
    // Get output shape
    auto output_tensor_info = output_tensors[0].GetTensorTypeAndShapeInfo();
    auto output_dims = output_tensor_info.GetShape();

    // Ensure we have exactly two output values
    size_t output_size = 1;
    for (auto dim : output_dims) {
        output_size *= (dim > 0) ? dim : 1;
    }

    if (output_size != 2) {
        throw std::runtime_error("Expected output size of 2, but got " + std::to_string(output_size));
    }

    return std::vector<float>(output_data, output_data + output_size);
}

    
