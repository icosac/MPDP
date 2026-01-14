#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <onnxruntime_cxx_api.h>
#ifdef USE_CUDA
#include <onnxruntime_c_api.h>
#endif


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
        OnnxModel(const std::string& model_path, bool use_gpu);
        ~OnnxModel();
    
        // Run inference with float input data for classification (single-task) model
        std::vector<float> run(const std::vector<float>& input_data);
    
        // Run inference with float input data for multi-task model
        std::pair<std::vector<float>, std::vector<float>> run_multitask(const std::vector<float>& input_data);
        
        bool isMultitaskModel() const { return is_multitask_model; }
    };
