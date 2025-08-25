#include "model.hpp"

OnnxModel::OnnxModel(const std::string& model_path) : env(ORT_LOGGING_LEVEL_WARNING, "onnx-model") {
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

OnnxModel::~OnnxModel() {
    // Free allocated strings
    for (auto name : input_names) {
        free((void*)name);
    }
    for (auto name : output_names) {
        free((void*)name);
    }
}