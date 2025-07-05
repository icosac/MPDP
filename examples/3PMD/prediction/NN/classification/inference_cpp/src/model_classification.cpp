#include "model.hpp"


// Run inference with float input data for classification (single-task) model
std::vector<float> OnnxModel::run(const std::vector<float>& input_data) {
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
    
    // std::cout << "Output shape: ";
    // for (auto dim : output_dims) {
    //     std::cout << dim << " ";
    // }
    // std::cout << std::endl;
    
    // Calculate total output size
    size_t output_size = 1;
    for (auto dim : output_dims) {
        output_size *= dim > 0 ? dim : 1; // Skip negative dimensions
    }
    
    // Copy output data to a vector
    std::vector<float> result(output_data, output_data + output_size);
    return result;
}
    
