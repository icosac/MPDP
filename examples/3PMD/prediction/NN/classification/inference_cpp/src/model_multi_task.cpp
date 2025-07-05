#include "model.hpp"



// Run inference with float input data for multi-task model
std::pair<std::vector<float>, std::vector<float>> OnnxModel::run_multitask(const std::vector<float>& input_data) {
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
