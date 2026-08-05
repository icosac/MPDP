#include "rs_triangular_onnx.hpp"

#include <rs.hh>
#include <cmath>
#include <numeric>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr std::size_t kRawFeatureCount = 3;

const std::vector<int>& default_class_labels()
{
    static const std::vector<int> labels = {
        1, 3, 6, 10, 11, 12, 14, 15, 16, 20, 21, 23,
        27, 31, 32, 33, 35, 39, 40, 43, 44, 47, 48};
    return labels;
}

struct Options {
    std::string model_path;
    std::string input;
    std::string labels_path;
    std::size_t topk = 3;
    std::size_t batch_size = 8192;
    int threads = 0;
    int cuda_device_id = 0;
    bool use_cuda = false;
    bool verbose = false;
};

struct Metrics {
    std::size_t total = 0;
    std::size_t correct = 0;
    std::size_t batches = 0;
    double inference_seconds = 0.0;
    double min_batch_seconds = std::numeric_limits<double>::max();
    double max_batch_seconds = 0.0;
    
    std::vector<int> unsolvable_cases;
    std::vector<bool> completely_unsolvable_cases;
    std::vector<bool> opt_man_in_top_k_cases;
    std::vector<std::pair<double, double>> errors;
    std::vector<double> ml_times_us;
};

bool parse_size(const char* value, std::size_t& out)
{
    char* end = nullptr;
    const unsigned long long parsed = std::strtoull(value, &end, 10);
    if (end == value || *end != '\0') {
        return false;
    }
    out = static_cast<std::size_t>(parsed);
    return true;
}

bool parse_int(const char* value, int& out)
{
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value || *end != '\0') {
        return false;
    }
    out = static_cast<int>(parsed);
    return true;
}

void print_usage(const char* exe)
{
    std::cerr
        << "Usage:\n"
        << "  " << exe << " <model.onnx> <thi,thf,kmax> [--topk N] [--labels labels.txt] [--threads N] [--cuda] [--cuda-device N] [--verbose]\n"
        << "  " << exe << " <model.onnx> <dataset.csv> [--topk N] [--batch-size N] [--labels labels.txt] [--threads N] [--cuda] [--cuda-device N] [--verbose]\n\n"
        << "Dataset CSV must contain columns: thi thf kmax ... man ...\n"
        << "The exported ONNX model expects raw [thi, thf, kmax]; trigonometric features and scaling are inside the graph.\n";
}

Options parse_args(int argc, char** argv)
{
    if (argc < 3) {
        print_usage(argv[0]);
        throw std::runtime_error("Missing required arguments.");
    }

    Options options;
    options.model_path = argv[1];
    options.input = argv[2];

    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--topk" && i + 1 < argc) {
            if (!parse_size(argv[++i], options.topk)) {
                throw std::runtime_error("Invalid --topk value.");
            }
        } else if (arg == "--batch-size" && i + 1 < argc) {
            if (!parse_size(argv[++i], options.batch_size)) {
                throw std::runtime_error("Invalid --batch-size value.");
            }
        } else if (arg == "--labels" && i + 1 < argc) {
            options.labels_path = argv[++i];
        } else if (arg == "--threads" && i + 1 < argc) {
            if (!parse_int(argv[++i], options.threads)) {
                throw std::runtime_error("Invalid --threads value.");
            }
        } else if (arg == "--cuda") {
            options.use_cuda = true;
        } else if (arg == "--cuda-device" && i + 1 < argc) {
            if (!parse_int(argv[++i], options.cuda_device_id)) {
                throw std::runtime_error("Invalid --cuda-device value.");
            }
        } else if (arg == "--verbose") {
            options.verbose = true;
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown or incomplete option: " + arg);
        }
    }

    if (options.topk == 0) {
        throw std::runtime_error("--topk must be greater than zero.");
    }
    if (options.batch_size == 0) {
        throw std::runtime_error("--batch-size must be greater than zero.");
    }
    if (options.cuda_device_id < 0) {
        throw std::runtime_error("--cuda-device must be non-negative.");
    }
    return options;
}

bool file_exists(const std::string& path)
{
    std::ifstream file(path);
    return file.good();
}

std::vector<float> parse_single_input(const std::string& text)
{
    std::string normalized = text;
    std::replace(normalized.begin(), normalized.end(), ',', ' ');

    std::istringstream stream(normalized);
    std::vector<float> values;
    float value = 0.0f;
    while (stream >> value) {
        values.push_back(value);
    }

    if (values.size() != kRawFeatureCount) {
        throw std::runtime_error("Single input must contain exactly thi, thf, kmax.");
    }
    return values;
}

std::vector<int> load_labels(const std::string& path)
{
    if (path.empty()) {
        return default_class_labels();
    }

    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Failed to open labels file: " + path);
    }

    std::vector<int> labels;
    int label = 0;
    while (file >> label) {
        labels.push_back(label);
    }
    if (labels.empty()) {
        throw std::runtime_error("Labels file is empty: " + path);
    }
    return labels;
}

std::vector<int> build_label_to_class(const std::vector<int>& labels)
{
    int max_label = 0;
    for (int label : labels) {
        if (label > max_label) {
            max_label = label;
        }
    }
    std::vector<int> label_to_class(static_cast<std::size_t>(max_label) + 1, -1);
    for (std::size_t class_idx = 0; class_idx < labels.size(); ++class_idx) {
        if (labels[class_idx] >= 0) {
            label_to_class[static_cast<std::size_t>(labels[class_idx])] = static_cast<int>(class_idx);
        }
    }
    return label_to_class;
}

std::vector<std::pair<int, float>> topk_logits(
    const float* logits,
    std::size_t class_count,
    const std::vector<int>& labels,
    std::size_t topk)
{
    std::vector<std::pair<int, float>> ranked;
    ranked.reserve(class_count);
    for (std::size_t class_idx = 0; class_idx < class_count; ++class_idx) {
        ranked.emplace_back(labels[class_idx], logits[class_idx]);
    }

    const std::size_t keep = std::min(topk, ranked.size());
    std::partial_sort(
        ranked.begin(),
        ranked.begin() + static_cast<std::ptrdiff_t>(keep),
        ranked.end(),
        [](const auto& a, const auto& b) {
            return a.second > b.second;
        });
    ranked.resize(keep);
    return ranked;
}

void print_class_probabilities(
    std::ostream& out,
    const float* logits,
    std::size_t class_count,
    const std::vector<int>& labels)
{
    double max_logit = -std::numeric_limits<double>::infinity();
    for (std::size_t class_idx = 0; class_idx < class_count; ++class_idx) {
        max_logit = std::max(max_logit, static_cast<double>(logits[class_idx]));
    }

    std::vector<double> exp_logits;
    exp_logits.reserve(class_count);
    double exp_sum = 0.0;
    for (std::size_t class_idx = 0; class_idx < class_count; ++class_idx) {
        const double exp_logit = std::exp(static_cast<double>(logits[class_idx]) - max_logit);
        exp_logits.push_back(exp_logit);
        exp_sum += exp_logit;
    }

    const auto old_flags = out.flags();
    const auto old_precision = out.precision();
    std::vector<std::pair<int, double>> probabilities;
    probabilities.reserve(class_count);
    for (std::size_t class_idx = 0; class_idx < class_count; ++class_idx) {
        const double probability = exp_sum > 0.0 ? exp_logits[class_idx] / exp_sum : 0.0;
        probabilities.emplace_back(labels[class_idx], probability);
    }
    std::sort(
        probabilities.begin(),
        probabilities.end(),
        [](const auto& a, const auto& b) {
            return a.second > b.second;
        });

    out << "All class probabilities (maneuver, probability, descending): ";
    out << std::scientific << std::setprecision(16);
    for (const auto& item : probabilities) {
        out << "(" << item.first << ", " << item.second << ") ";
    }
    out << std::endl;
    out.flags(old_flags);
    out.precision(old_precision);
}

void evaluate_prediction(
    const float* logits,
    std::size_t class_count,
    const std::vector<int>& label_to_class,
    const std::vector<int>& class_labels,
    int expected,
    std::size_t topk,
    const float* raw_features,
    float optimal_len,
    Metrics& metrics)
{
    metrics.unsolvable_cases.push_back(0);

    // Check if expected label is in top-k
    bool expected_in_topk = false;
    if (expected >= 0 && static_cast<std::size_t>(expected) < label_to_class.size()) {
        const int expected_class = label_to_class[static_cast<std::size_t>(expected)];
        if (expected_class >= 0) {
            const float expected_logit = logits[expected_class];
            std::size_t greater_count = 0;
            expected_in_topk = true;
            for (std::size_t class_idx = 0; class_idx < class_count; ++class_idx) {
                if (logits[class_idx] > expected_logit && ++greater_count >= topk) {
                    expected_in_topk = false;
                    break;
                }
            }
        }
    }

    std::vector<std::pair<int, float>> ranked = topk_logits(logits, class_count, class_labels, topk);
    
    Configuration2 ci(-1.0, 0.0, static_cast<double>(raw_features[0])); // thi
    Configuration2 cf(1.0, 0.0, static_cast<double>(raw_features[1]));  // thf
    double kmax = static_cast<double>(raw_features[2]);

    double best_length = std::numeric_limits<double>::infinity();
    int best_man = -1;
    bool found_valid = false;
    bool found_optimal = expected_in_topk;

    for (const auto& r : ranked) {
        int predicted_maneuver = r.first;
        RS myRS = RS(ci, cf, { kmax, static_cast<double>(predicted_maneuver) });
        myRS.solve();
        if (!std::isfinite(myRS.l()) || myRS.l() >= 1e99) {
            metrics.unsolvable_cases.back() += 1;
        } else {
            if (myRS.l() < best_length) {
                best_length = myRS.l();
                best_man = predicted_maneuver;
            }
            found_valid = true;
        }
        if (myRS.getNman() == expected) {
            found_optimal = true;
        }
        if (std::isfinite(myRS.l()) && myRS.l() < 1e99 && std::abs(myRS.l() - static_cast<double>(optimal_len)) <= 1e-6) {
            found_optimal = true;
        }
    }

    if (!found_optimal){
        std::cerr << "Warning: Optimal maneuver not found in top-" << topk << " predictions. Expected: " << expected << ", Optimal length: " << optimal_len << ", Best predicted length: " << best_length << std::endl;
        std::cerr << "thi: " << raw_features[0] << ", thf: " << raw_features[1] << ", kmax: " << raw_features[2] << std::endl;
        std::cerr << "Top-" << topk << " predictions (maneuver, logit): ";
        for (const auto& r : ranked) {
            std::cerr << "(" << r.first << ", " << r.second << ") ";
        }        
        std::cerr << std::endl;
        std::cerr << "Optimal maneuver: " << expected << ", Optimal length: " << optimal_len << std::endl;
        std::cerr << "Best predicted length: " << best_length << " with maneuver " << best_man << std::endl;
        print_class_probabilities(std::cerr, logits, class_count, class_labels);
    }

    metrics.completely_unsolvable_cases.push_back(!found_valid);
    metrics.opt_man_in_top_k_cases.push_back(found_optimal);
    metrics.errors.push_back({optimal_len, best_length});

    if (found_optimal) {
        metrics.correct += 1;
    }
}

bool parse_dataset_line(const std::string& line, float* raw_features, float& len, int& label)
{
    if (line.empty()) {
        return false;
    }

    const char* cursor = line.c_str();
    char* next = nullptr;

    raw_features[0] = std::strtof(cursor, &next);
    if (next == cursor) {
        return false;
    }
    cursor = next;

    raw_features[1] = std::strtof(cursor, &next);
    if (next == cursor) {
        return false;
    }
    cursor = next;

    raw_features[2] = std::strtof(cursor, &next);
    if (next == cursor) {
        return false;
    }
    cursor = next;

    len = std::strtof(cursor, &next);
    if (next == cursor) {
        return false;
    }
    cursor = next;

    const long parsed_label = std::strtol(cursor, &next, 10);
    if (next == cursor) {
        return false;
    }
    label = static_cast<int>(parsed_label);
    return true;
}

void score_batch(
    RSTriangularOnnxModel& model,
    std::vector<float>& batch_features,
    const std::vector<float>& batch_lens,
    const std::vector<int>& batch_labels,
    std::size_t batch_count,
    const std::vector<int>& label_to_class,
    const std::vector<int>& class_labels,
    std::size_t topk,
    Metrics& metrics)
{
    const auto start = std::chrono::steady_clock::now();
    std::vector<Ort::Value> outputs = model.run(batch_features.data(), batch_count);
    const auto end = std::chrono::steady_clock::now();

    const double seconds = std::chrono::duration<double>(end - start).count();
    metrics.inference_seconds += seconds;
    metrics.min_batch_seconds = std::min(metrics.min_batch_seconds, seconds);
    metrics.max_batch_seconds = std::max(metrics.max_batch_seconds, seconds);
    ++metrics.batches;

    const float* logits = outputs[0].GetTensorData<float>();
    const std::size_t class_count = model.output_width();

    double ml_time_per_sample_us = (seconds * 1'000'000.0) / batch_count;
    
    for (std::size_t row = 0; row < batch_count; ++row) {
        metrics.ml_times_us.push_back(ml_time_per_sample_us);
        evaluate_prediction(
            logits + row * class_count,
            class_count,
            label_to_class,
            class_labels,
            batch_labels[row],
            topk,
            batch_features.data() + row * kRawFeatureCount,
            batch_lens[row],
            metrics);
    }
    metrics.total += batch_count;
}

void print_summary(const Metrics& metrics, std::size_t top_k) {
    if (metrics.total == 0) return;
    
    std::cout << std::fixed << std::setprecision(6);
    
    // Print the summary of unsolvable cases
    std::cout << "=========== Summary of unsolvable cases ===========\n";
    double avg_unsolvable = std::accumulate(metrics.unsolvable_cases.begin(), metrics.unsolvable_cases.end(), 0.0) / metrics.unsolvable_cases.size();
    std::cout << "Average unsolvable manouevres per test: " << avg_unsolvable << " out of " << top_k << std::endl;
    int total_unsolvable = std::accumulate(metrics.unsolvable_cases.begin(), metrics.unsolvable_cases.end(), 0);
    std::cout << "Unsolvable manouevres: " << total_unsolvable << " out of " << metrics.total * top_k << std::endl;
    int total_completely_unsolvable = std::accumulate(metrics.completely_unsolvable_cases.begin(), metrics.completely_unsolvable_cases.end(), 0);
    std::cout << "Completely unsolvable tests: " << total_completely_unsolvable << " out of " << metrics.total << std::endl;
    
    // Print the summary of errors, average error, standard devition and also average percentage of error
    std::cout << "=========== Summary of length error ===========\n";
    double total_error = 0.0;
    double total_percentage_error = 0.0;
    for (const auto& error_pair : metrics.errors) {
        double gt_length = error_pair.first;
        double predicted_length = error_pair.second;
        if (!std::isfinite(predicted_length) || predicted_length >= 1e99) {
            continue;
        }
        double error = std::abs(predicted_length - gt_length);
        double percentage_error = (error / gt_length) * 100.0;
        total_error += error;
        total_percentage_error += percentage_error;
    }
    double average_error = total_error / metrics.errors.size();
    double average_percentage_error = total_percentage_error / metrics.errors.size();
    double std_dev_error = 0.0;
    double std_dev_percentage_error = 0.0;
    for (const auto& error_pair : metrics.errors) {
        double gt_length = error_pair.first;
        double predicted_length = error_pair.second;
        if (!std::isfinite(predicted_length) || predicted_length >= 1e99) {
            continue;
        }
        double error = predicted_length - gt_length;
        double percentage_error = (error / gt_length) * 100.0;
        std_dev_error += (error - average_error) * (error - average_error);
        std_dev_percentage_error += (percentage_error - average_percentage_error) * (percentage_error - average_percentage_error);
    }

    int n_wrong_guesses = std::count_if(metrics.errors.begin(), metrics.errors.end(), [](const std::pair<double, double>& error_pair) {
        double gt_length = error_pair.first;
        double predicted_length = error_pair.second;
        if (!std::isfinite(predicted_length) || predicted_length >= 1e99) {
            return false;
        }
        double error = std::abs(predicted_length - gt_length);
        return error > 1e-6;
    });

    double max_error_percent = 0.0;
    for (const auto& error_pair : metrics.errors) {
        double gt_length = error_pair.first;
        double predicted_length = error_pair.second;
        if (!std::isfinite(predicted_length) || predicted_length >= 1e99) {
            continue;
        }
        double error = std::abs(predicted_length - gt_length);
        double percentage_error = (error / gt_length) * 100.0;
        if (percentage_error > max_error_percent) {
            max_error_percent = percentage_error;
        }
    }

    std_dev_error = std::sqrt(std_dev_error / metrics.errors.size());
    std_dev_percentage_error = std::sqrt(std_dev_percentage_error / metrics.errors.size());
    std::cout << "Average error: " << average_error << std::endl;
    std::cout << "Average percentage error: " << average_percentage_error << "%" << std::endl;
    std::cout << "Average error on wrong guesses: " << (n_wrong_guesses ? (average_error * metrics.total / n_wrong_guesses) : 0.0) << std::endl;
    std::cout << "Average percentage error on wrong guesses: " << (n_wrong_guesses ? (average_percentage_error * metrics.total / n_wrong_guesses) : 0.0) << "%" << std::endl;
    std::cout << "Standard deviation of error: " << std_dev_error << std::endl;
    std::cout << "Standard deviation of percentage error: " << std_dev_percentage_error << "%" << std::endl;
    std::cout << "Maximum percentage error: " << max_error_percent << "%" << std::endl;

    // Print the summary of times
    std::cout << "=========== Summary of computational times ===========\n";
    double total_ml_time = 0.0;
    for (const auto& t : metrics.ml_times_us) {
        total_ml_time += t;
    }
    double average_ml_time = total_ml_time / metrics.ml_times_us.size();
    double std_dev_ml = 0.0;
    for (const auto& t : metrics.ml_times_us) {
        std_dev_ml += (t - average_ml_time) * (t - average_ml_time);
    }
    std_dev_ml = std::sqrt(std_dev_ml / metrics.ml_times_us.size());
    
    std::cout << "Average ML time: " << average_ml_time << " microseconds, Standard deviation: " << std_dev_ml << std::endl;

    // Print the summary of optimal manoeuvre in top_k cases
    std::cout << "=========== Summary of optimal manoeuvre in top_k cases ===========\n";
    int total_optimal_in_top_k = std::accumulate(metrics.opt_man_in_top_k_cases.begin(), metrics.opt_man_in_top_k_cases.end(), 0);
    std::cout << "Optimal manoeuvre in top_k cases: " << total_optimal_in_top_k << " out of " << metrics.total << std::endl;
    std::cout << "Unoptimal manoeuvre in top_k cases: " << (metrics.total - total_optimal_in_top_k) << std::endl;
    double percentage_optimal_in_top_k = (static_cast<double>(total_optimal_in_top_k) / metrics.total) * 100.0;
    std::cout << "Percentage of optimal manoeuvre in top_k cases: " << percentage_optimal_in_top_k << "%" << std::endl;
}

Metrics run_dataset(
    RSTriangularOnnxModel& model,
    const std::string& dataset_path,
    const std::vector<int>& class_labels,
    std::size_t topk,
    std::size_t batch_size)
{
    std::ifstream file(dataset_path);
    if (!file) {
        throw std::runtime_error("Failed to open dataset: " + dataset_path);
    }

    std::string line;
    std::getline(file, line);

    std::vector<float> batch_features(batch_size * kRawFeatureCount);
    std::vector<float> batch_lens(batch_size);
    std::vector<int> batch_labels(batch_size);
    const std::vector<int> label_to_class = build_label_to_class(class_labels);
    Metrics metrics;
    std::size_t batch_count = 0;

    while (std::getline(file, line)) {
        int label = 0;
        float len = 0.0f;
        float raw[kRawFeatureCount] = {};
        if (!parse_dataset_line(line, raw, len, label)) {
            continue;
        }

        float* dst = batch_features.data() + batch_count * kRawFeatureCount;
        dst[0] = raw[0];
        dst[1] = raw[1];
        dst[2] = raw[2];
        batch_lens[batch_count] = len;
        batch_labels[batch_count] = label;
        ++batch_count;

        if (batch_count == batch_size) {
            score_batch(model, batch_features, batch_lens, batch_labels, batch_count, label_to_class, class_labels, topk, metrics);
            batch_count = 0;
        }
    }

    if (batch_count > 0) {
        score_batch(model, batch_features, batch_lens, batch_labels, batch_count, label_to_class, class_labels, topk, metrics);
    }
    return metrics;
}

void run_single(RSTriangularOnnxModel& model, const std::string& input, const std::vector<int>& class_labels, std::size_t topk)
{
    std::vector<float> raw_features = parse_single_input(input);
    const auto start = std::chrono::steady_clock::now();
    std::vector<Ort::Value> outputs = model.run(raw_features.data(), 1);
    const auto end = std::chrono::steady_clock::now();
    const double inference_seconds = std::chrono::duration<double>(end - start).count();
    const float* logits = outputs[0].GetTensorData<float>();

    const auto ranked = topk_logits(logits, model.output_width(), class_labels, topk);
    std::cout << "Top-" << ranked.size() << " maneuver logits:\n";
    for (const auto& item : ranked) {
        std::cout << "  man " << item.first << ": " << std::fixed << std::setprecision(6) << item.second << '\n';
    }
    print_class_probabilities(std::cout, logits, model.output_width(), class_labels);
    std::cout << std::fixed << std::setprecision(9);
    std::cout << "Inference time: " << inference_seconds << " s\n";
    std::cout << "Inference latency: " << inference_seconds * 1000.0 << " ms\n";
    std::cout << "Inference latency: " << inference_seconds * 1000000.0 << " us\n";
    std::cout << "Inference throughput: "
              << (inference_seconds > 0.0 ? 1.0 / inference_seconds : 0.0)
              << " samples/s\n";
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options options = parse_args(argc, argv);
        RSTriangularOnnxModel model(
            options.model_path,
            options.threads,
            options.verbose,
            options.use_cuda,
            options.cuda_device_id);
        std::vector<int> class_labels = load_labels(options.labels_path);

        if (model.input_width() != kRawFeatureCount) {
            throw std::runtime_error("Model input width is not 3. Expected raw [thi, thf, kmax].");
        }
        if (class_labels.size() != model.output_width()) {
            throw std::runtime_error("Class label count does not match model output width.");
        }

        const std::size_t topk = std::min(options.topk, class_labels.size());
        if (file_exists(options.input)) {
            const auto wall_start = std::chrono::steady_clock::now();
            const Metrics metrics = run_dataset(model, options.input, class_labels, topk, options.batch_size);
            const auto wall_end = std::chrono::steady_clock::now();
            const double wall_seconds = std::chrono::duration<double>(wall_end - wall_start).count();

            print_summary(metrics, topk);
        } else {
            run_single(model, options.input, class_labels, topk);
        }

        return 0;
    } catch (const Ort::Exception& error) {
        std::cerr << "ONNX Runtime error: " << error.what() << '\n';
        return 1;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
