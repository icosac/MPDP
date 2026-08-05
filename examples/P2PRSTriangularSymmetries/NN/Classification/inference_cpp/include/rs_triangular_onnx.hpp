#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <onnxruntime_cxx_api.h>

class RSTriangularOnnxModel {
public:
    explicit RSTriangularOnnxModel(
        const std::string& model_path,
        int intra_op_threads = 0,
        bool verbose = false,
        bool use_cuda = false,
        int cuda_device_id = 0);

    std::vector<Ort::Value> run(float* input_data, std::size_t sample_count);

    std::size_t input_width() const { return input_width_; }
    std::size_t output_width() const { return output_width_; }
    const std::vector<std::int64_t>& input_shape() const { return input_shape_; }
    const std::vector<std::int64_t>& output_shape() const { return output_shape_; }

private:
    Ort::Env env_;
    Ort::SessionOptions session_options_;
    Ort::Session session_{nullptr};
    Ort::MemoryInfo memory_info_;
    Ort::RunOptions run_options_;

    std::vector<std::string> input_name_storage_;
    std::vector<std::string> output_name_storage_;
    std::vector<const char*> input_names_;
    std::vector<const char*> output_names_;

    std::vector<std::int64_t> input_shape_;
    std::vector<std::int64_t> output_shape_;
    std::vector<std::int64_t> current_input_shape_;

    std::size_t input_width_ = 0;
    std::size_t output_width_ = 0;
};
