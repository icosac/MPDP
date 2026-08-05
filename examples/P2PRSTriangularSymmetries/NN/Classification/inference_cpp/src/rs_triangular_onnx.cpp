#include "rs_triangular_onnx.hpp"

#include <iostream>
#include <stdexcept>

namespace {

std::size_t static_width_from_shape(const std::vector<std::int64_t>& shape, const char* tensor_name)
{
    if (shape.size() != 2 || shape[1] <= 0) {
        throw std::runtime_error(std::string("Expected ") + tensor_name + " shape [batch, width].");
    }
    return static_cast<std::size_t>(shape[1]);
}

void print_shape(const char* label, const std::vector<std::int64_t>& shape)
{
    std::cout << label;
    for (std::int64_t dim : shape) {
        std::cout << ' ' << dim;
    }
    std::cout << '\n';
}

}  // namespace

RSTriangularOnnxModel::RSTriangularOnnxModel(
    const std::string& model_path,
    int intra_op_threads,
    bool verbose,
    bool use_cuda,
    int cuda_device_id)
    : env_(ORT_LOGGING_LEVEL_WARNING, "rs-triangular-classifier"),
      memory_info_(Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault))
{
    session_options_.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    session_options_.SetExecutionMode(ExecutionMode::ORT_SEQUENTIAL);
    session_options_.EnableCpuMemArena();
    session_options_.EnableMemPattern();
    if (intra_op_threads > 0) {
        session_options_.SetIntraOpNumThreads(intra_op_threads);
    }
    session_options_.SetInterOpNumThreads(1);

#ifdef RS_TRIANGULAR_USE_CUDA
    if (use_cuda) {
        OrtCUDAProviderOptions cuda_options{};
        cuda_options.device_id = cuda_device_id;
        session_options_.AppendExecutionProvider_CUDA(cuda_options);
        if (verbose) {
            std::cout << "Using ONNX Runtime CUDA execution provider on device "
                      << cuda_device_id << '\n';
        }
    }
#else
    if (use_cuda) {
        throw std::runtime_error(
            "This executable was built without CUDA support. "
            "Reconfigure with -DRS_TRIANGULAR_USE_CUDA=ON.");
    }
#endif

    session_ = Ort::Session(env_, model_path.c_str(), session_options_);

    Ort::AllocatorWithDefaultOptions allocator;
    const std::size_t input_count = session_.GetInputCount();
    const std::size_t output_count = session_.GetOutputCount();
    if (input_count != 1 || output_count != 1) {
        throw std::runtime_error("Expected a single-input, single-output classification ONNX model.");
    }

    input_name_storage_.reserve(input_count);
    output_name_storage_.reserve(output_count);
    input_names_.reserve(input_count);
    output_names_.reserve(output_count);

    for (std::size_t i = 0; i < input_count; ++i) {
        auto name = session_.GetInputNameAllocated(i, allocator);
        input_name_storage_.emplace_back(name.get());
        input_names_.push_back(input_name_storage_.back().c_str());
    }
    for (std::size_t i = 0; i < output_count; ++i) {
        auto name = session_.GetOutputNameAllocated(i, allocator);
        output_name_storage_.emplace_back(name.get());
        output_names_.push_back(output_name_storage_.back().c_str());
    }

    input_shape_ = session_.GetInputTypeInfo(0).GetTensorTypeAndShapeInfo().GetShape();
    output_shape_ = session_.GetOutputTypeInfo(0).GetTensorTypeAndShapeInfo().GetShape();
    input_width_ = static_width_from_shape(input_shape_, "input");
    output_width_ = static_width_from_shape(output_shape_, "output");
    current_input_shape_ = input_shape_;

    if (verbose) {
        std::cout << "Loaded " << model_path << '\n';
        std::cout << "Input name: " << input_names_[0] << '\n';
        print_shape("Input shape:", input_shape_);
        std::cout << "Output name: " << output_names_[0] << '\n';
        print_shape("Output shape:", output_shape_);
    }
}

std::vector<Ort::Value> RSTriangularOnnxModel::run(float* input_data, std::size_t sample_count)
{
    current_input_shape_[0] = static_cast<std::int64_t>(sample_count);
    const std::size_t input_value_count = sample_count * input_width_;

    Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
        memory_info_,
        input_data,
        input_value_count,
        current_input_shape_.data(),
        current_input_shape_.size());

    return session_.Run(
        run_options_,
        input_names_.data(),
        &input_tensor,
        1,
        output_names_.data(),
        output_names_.size());
}
