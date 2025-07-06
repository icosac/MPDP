#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <onnxruntime_cxx_api.h>


class OnnxRegressionModel {
    private:
        Ort::Env env;
        Ort::SessionOptions session_options;
        Ort::Session session{nullptr};
        std::vector<const char*> input_names;
        std::vector<const char*> output_names;
        std::vector<std::vector<int64_t>> input_node_dims;
        std::vector<std::vector<int64_t>> output_node_dims;
    
    public:
        OnnxRegressionModel(const std::string& model_path);
        ~OnnxRegressionModel();
    
        // Run inference with float input data for classification (single-task) model
        std::vector<float> run(const std::vector<float>& input_data);
    
        // Run inference with float input data for multi-task model
        std::pair<std::vector<float>, std::vector<float>> run_multitask(const std::vector<float>& input_data);
    };