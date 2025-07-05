#include <torch/script.h>
#include <torch/torch.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <Eigen/Dense>
#include <nlohmann/json.hpp> // Include a JSON library like nlohmann/json
#include <fstream>

torch::Device device(torch::cuda::is_available() ? torch::kCUDA : torch::kCPU);

std::vector<float> parseInputString(const std::string& input) {
    std::vector<float> values;
    std::stringstream ss(input);
    std::string item;
    while (std::getline(ss, item, ',')) {
        values.push_back(std::stof(item));
    }
    return values;
}

Eigen::VectorXf applyScaler(const Eigen::VectorXf& input, const Eigen::VectorXf& mean, const Eigen::VectorXf& scale) {
    return (input - mean).cwiseQuotient(scale);
}

std::pair<Eigen::VectorXf, Eigen::VectorXf> loadScaler(const std::string& scaler_json_path) {
    std::ifstream file(scaler_json_path);
    if (!file) {
        throw std::runtime_error("Error opening scaler JSON file: " + scaler_json_path);
    }

    nlohmann::json scaler_json;
    file >> scaler_json;

    std::vector<float> mean_vec = scaler_json["mean"];
    std::vector<float> scale_vec = scaler_json["scale"];

    Eigen::VectorXf mean = Eigen::Map<Eigen::VectorXf>(mean_vec.data(), mean_vec.size());
    Eigen::VectorXf scale = Eigen::Map<Eigen::VectorXf>(scale_vec.data(), scale_vec.size());

    return {mean, scale};
}

std::tuple<torch::Tensor, torch::Tensor> predict(torch::jit::script::Module& model, torch::Tensor input, int N) {
    model.to(device);
    model.eval();
    input = input.to(device);
    
    torch::Tensor output = model.forward({input}).toTensor();
    torch::Tensor probabilities = torch::softmax(output, 1);
    torch::Tensor top_probs, top_indices;
    std::tie(top_probs, top_indices) = torch::topk(probabilities, N, 1);

    return {top_indices.cpu(), top_probs.cpu()};
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_values_or_csv_file>" << std::endl;
        return 1;
    }

    std::string input_arg = argv[1];
    std::string scaler_json_path = "models/model1_scaler.json"; // Path to the scaler JSON file

    Eigen::VectorXf mean, scale;
    try {
        std::tie(mean, scale) = loadScaler(scaler_json_path);
    } catch (const std::exception& e) {
        std::cerr << "Error loading scaler: " << e.what() << std::endl;
        return -1;
    }

    torch::jit::script::Module model;

    try {
        model = torch::jit::load("models/model_cpp.pt", device);
    } catch (const c10::Error& e) {
        std::cerr << "Error loading the model." << std::endl;
        std::cerr << e.what() << std::endl;
        return -1;
    }

    if (input_arg.find(".csv") != std::string::npos) {
        std::ifstream file(input_arg);
        if (!file) {
            std::cerr << "Error opening file: " << input_arg << std::endl;
            return -1;
        }
        
        float kmax, theta_i, theta_f, alpha_m, alpha_f;
        int id_man_comb;

        std::vector<std::vector<float>> data;
        std::vector<int> labels;

        std::string tmp;
        std::getline(file, tmp); // Skip header

        int counter = 0;

        while(file >> kmax >> theta_i >> theta_f >> alpha_m >> alpha_f >> id_man_comb) {
            Eigen::VectorXf input(5);
            input << kmax, theta_i, theta_f, alpha_m, alpha_f;
            Eigen::VectorXf scaled_input = applyScaler(input, mean, scale);

            data.push_back(std::vector<float>(scaled_input.data(), scaled_input.data() + scaled_input.size()));
            labels.push_back(id_man_comb);
            // if (counter == 5) break;
            counter++;
        }

        torch::Tensor input_tensor = torch::from_blob(data.data(), {(long)data.size(), (long)data[0].size()}, torch::kFloat32);
        
        auto [top_classes, top_probs] = predict(model, input_tensor, 1);
        // std::cout << "Top-N Predicted Classes: " << std::endl << top_classes << std::endl;
        // std::cout << "Top-N Probabilities: " << std::endl << top_probs << std::endl;

        int correct = 0;
        int total = 0;
        for (int i = 0; i < labels.size(); i++) {
            for (int j = 0; j < top_classes.size(1); j++) {
                if (labels[i] == top_classes[i][j].item<int>()) {
                    correct++;
                    break;
                }
            }
        }
        std::cout << "Accuracy: " << (float)correct / labels.size() << std::endl;

    } else {
        std::vector<float> input_values = parseInputString(input_arg);
        Eigen::VectorXf input = Eigen::Map<Eigen::VectorXf>(input_values.data(), input_values.size());
        Eigen::VectorXf scaled_input = applyScaler(input, mean, scale);

        torch::Tensor input_tensor = torch::from_blob(scaled_input.data(), {(long)scaled_input.size()}, torch::kFloat32).unsqueeze(0);
        
        auto [top_classes, top_probs] = predict(model, input_tensor, 5);
        std::cout << "Top-N Predicted Classes: " << top_classes << std::endl;
        std::cout << "Top-N Probabilities: " << top_probs << std::endl;
    }

    return 0;
}
