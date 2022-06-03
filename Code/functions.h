// data.h
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <algorithm>
#include <chrono>
#include <ratio>
#include "omp.h"
#include "variance_est.h"
#include "hash_chain.h"
#include "relative_error.h"
#include "absolute_error.h"
#include "hash_chain_kwise.h"
#include "xoroshiro64star.h"
#include "pred_pair.h"

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

static inline unsigned my_sign(const double x)
{
    return (*(unsigned long long*)&x >> 63);
}

uint8_t fast_mod(const uint8_t input, const uint8_t ceil) {
    // apply the modulo operator only when needed
    // (i.e. when the input is greater than the ceiling)
    return input >= ceil ? input % ceil : input;
    // NB: the assumption here is that the numbers are positive
}

bool cmp_low(const Feature& a, const Feature& b)
{
    return a.low_ > b.low_;
}

bool cmp_up(const Feature& a, const Feature& b)
{
    return a.up_ > b.up_;
}

bool cmp_avg(const Feature& a, const Feature& b)
{
    return a.avg_ > b.avg_;
}

bool cmp_ratio(const Feature& a, const Feature& b)
{
    return a.ratio_ < b.ratio_;
}

bool cmp_diff(const Feature& a, const Feature& b)
{
    return a.diff_ < b.diff_;
}

bool cmp_vec(const int a, const int b) 
{
    return a < b;
}

bool cmp_idx_value_pair(const std::pair<int, double>& a, const std::pair<int, double>& b)
{
    return a.first < b.first;
}

void ReadFile(std::string file_name, std::vector<std::vector<double>>& file_array)
{
    std::vector<std::vector<double>> temp_file_array;
    std::ifstream in_file(file_name, std::ios::in);
    std::string line_str;
    getline(in_file, line_str);
    int count = 1;
    while (getline(in_file, line_str))
    {
        //printf("%d\n", count);
        count++;
        //cout << line_str << endl; 
        std::stringstream ss(line_str);
        std::string str;
        std::vector<double> line_array;  
        while (getline(ss, str, ','))
            line_array.push_back(stod(str));
        temp_file_array.push_back(line_array);
    }
    in_file.close();
    file_array.resize(temp_file_array[0].size());
    for (int i = 0; i < file_array.size(); ++i)
    {
        file_array[i].resize(temp_file_array.size());
    }
    for (int i = 0; i < file_array.size(); ++i)
    {
        for (int j = 0; j < file_array[i].size(); ++j)
        {
            file_array[i][j] = temp_file_array[j][i];
        } 
    }
}

void ReadFileHash(std::string file_name, std::vector<std::vector<double>>& file_array)
{
    std::vector<std::vector<double>> temp_file_array, hash_file_array;
    std::ifstream in_file(file_name, std::ios::in);
    std::string line_str;
    getline(in_file, line_str);
    int count = 1;
    while (getline(in_file, line_str))
    {
        //printf("%d\n", count);
        count++;
        //cout << line_str << endl; 
        std::stringstream ss(line_str);
        std::string str;
        std::vector<double> line_array;  
        while (getline(ss, str, ','))
            line_array.push_back(stod(str));
        temp_file_array.push_back(line_array);
    }
    in_file.close();
    int num_record = temp_file_array.size();
    std::vector<int> permutation_vec(num_record);
    for (int i = 0; i < permutation_vec.size(); i++)
    {
        permutation_vec[i] = i;
    }
    HashChain hash_chain(4, 25);
    for (int i = 0; i < permutation_vec.size(); i++)
    {
        hash_chain.Insert(permutation_vec[i]);
    }
    for (int i = 0; i < hash_chain.table_.size(); i++)
    {
        for (std::list<int>::const_iterator iter = hash_chain.table_[i].begin(); iter != hash_chain.table_[i].end(); iter++)
        {
            hash_file_array.push_back(temp_file_array[*iter]);
        }
    }
    printf("After Hashing NumofRecord %d\n", int(hash_file_array.size()));
    file_array.resize(hash_file_array[0].size());
    for (int i = 0; i < file_array.size(); ++i)
    {
        file_array[i].resize(hash_file_array.size());
    }
    for (int i = 0; i < file_array.size(); ++i)
    {
        for (int j = 0; j < file_array[0].size(); ++j)
        {
            file_array[i][j] = hash_file_array[j][i];
        } 
    }
}

void ReadFileHashUseIndex(std::string file_name, std::vector<std::vector<double>>& file_array)
{
    std::vector<std::vector<double>> temp_file_array;
    std::ifstream in_file(file_name, std::ios::in);
    std::string line_str;
    getline(in_file, line_str);
    int count = 1;
    while (getline(in_file, line_str))
    {
        //printf("%d\n", count);
        count++;
        //cout << line_str << endl; 
        std::stringstream ss(line_str);
        std::string str;
        std::vector<double> line_array;  
        while (getline(ss, str, ','))
            line_array.push_back(stod(str));
        temp_file_array.push_back(line_array);
    }
    in_file.close();
    int num_record = temp_file_array.size();
    std::vector<int> permutation_vec(num_record);
    for (int i = 0; i < permutation_vec.size(); i++)
    {
        permutation_vec[i] = i;
    }
    HashChain hash_chain(4, 25);
    for (int i = 0; i < permutation_vec.size(); i++)
    {
        hash_chain.Insert(permutation_vec[i]);
    }
    int permutation_idx = 0;
    for (int i = 0; i < hash_chain.table_.size(); i++)
    {
        for (std::list<int>::const_iterator iter = hash_chain.table_[i].begin(); iter != hash_chain.table_[i].end(); iter++)
        {
            permutation_vec[permutation_idx] = *iter;
            permutation_idx ++;
        }
    }
    printf("After Hashing NumofRecord %d\n", int(permutation_vec.size()));
    file_array.resize(temp_file_array[0].size());
    for (int i = 0; i < file_array.size(); ++i)
    {
        file_array[i].resize(temp_file_array.size());
    }
    for (int i = 0; i < file_array.size(); ++i)
    {
        for (int j = 0; j < file_array[i].size(); ++j)
        {
            file_array[i][j] = temp_file_array[permutation_vec[j]][i];
        } 
    }
}

void AlignedReadFileHashUseIndexAddPushBack(std::string file_name, std::vector<std::vector<double>>& file_array)
{
    std::vector<std::vector<double>> temp_file_array, hash_file_array;
    std::ifstream in_file(file_name, std::ios::in);
    std::string line_str;
    getline(in_file, line_str);
    int count = 1;
    while (getline(in_file, line_str))
    {
        //printf("%d\n", count);
        count++;
        //cout << line_str << endl; 
        std::stringstream ss(line_str);
        std::string str;
        std::vector<double> line_array;  
        while (getline(ss, str, ','))
            line_array.push_back(stod(str));
        temp_file_array.push_back(line_array);
    }
    in_file.close();
    int num_record = temp_file_array.size();
    std::vector<int> permutation_vec(num_record);
    for (int i = 0; i < permutation_vec.size(); i++)
    {
        permutation_vec[i] = i;
    }
    HashChain hash_chain(4, 25);
    for (int i = 0; i < permutation_vec.size(); i++)
    {
        hash_chain.Insert(permutation_vec[i]);
    }
    int permutation_idx = 0;
    for (int i = 0; i < hash_chain.table_.size(); i++)
    {
        for (std::list<int>::const_iterator iter = hash_chain.table_[i].begin(); iter != hash_chain.table_[i].end(); iter++)
        {
            permutation_vec[permutation_idx] = *iter;
            hash_file_array.push_back(temp_file_array[*iter]);
            permutation_idx ++;
        }
    }
    printf("After Hashing NumofRecord %d\n", int(permutation_vec.size()));
    file_array.resize(temp_file_array[0].size());
    for (int i = 0; i < file_array.size(); ++i)
    {
        file_array[i].resize(temp_file_array.size());
    }
    for (int i = 0; i < file_array.size(); ++i)
    {
        for (int j = 0; j < file_array[i].size(); ++j)
        {
            file_array[i][j] = temp_file_array[permutation_vec[j]][i];
        } 
    }
}

void AlignedReadFileHashStatic(std::string file_name, std::vector<std::vector<double>>& file_array)
{
    static std::vector<std::vector<double>> temp_file_array, hash_file_array;
    std::ifstream in_file(file_name, std::ios::in);
    std::string line_str;
    getline(in_file, line_str);
    int count = 1;
    while (getline(in_file, line_str))
    {
        //printf("%d\n", count);
        count++;
        //cout << line_str << endl; 
        std::stringstream ss(line_str);
        std::string str;
        std::vector<double> line_array;  
        while (getline(ss, str, ','))
            line_array.push_back(stod(str));
        temp_file_array.push_back(line_array);
    }
    in_file.close();
    int num_record = temp_file_array.size();
    std::vector<int> permutation_vec(num_record);
    for (int i = 0; i < permutation_vec.size(); i++)
    {
        permutation_vec[i] = i;
    }
    HashChain hash_chain(4, 25);
    for (int i = 0; i < permutation_vec.size(); i++)
    {
        hash_chain.Insert(permutation_vec[i]);
    }
    for (int i = 0; i < hash_chain.table_.size(); i++)
    {
        for (std::list<int>::const_iterator iter = hash_chain.table_[i].begin(); iter != hash_chain.table_[i].end(); iter++)
        {
            hash_file_array.push_back(temp_file_array[*iter]);
        }
    }
    printf("After Hashing NumofRecord %d\n", int(hash_file_array.size()));
    file_array.resize(hash_file_array[0].size());
    for (int i = 0; i < file_array.size(); ++i)
    {
        file_array[i].resize(hash_file_array.size());
    }
    for (int i = 0; i < file_array.size(); ++i)
    {
        for (int j = 0; j < file_array[0].size(); ++j)
        {
            file_array[i][j] = hash_file_array[j][i];
        } 
    }
}

void AlignedReadFileHashUseIndexStatic(std::string file_name, std::vector<std::vector<double>>& file_array)
{
    static std::vector<std::vector<double>> temp_file_array;
    std::ifstream in_file(file_name, std::ios::in);
    std::string line_str;
    getline(in_file, line_str);
    int count = 1;
    while (getline(in_file, line_str))
    {
        //printf("%d\n", count);
        count++;
        //cout << line_str << endl; 
        std::stringstream ss(line_str);
        std::string str;
        std::vector<double> line_array;  
        while (getline(ss, str, ','))
            line_array.push_back(stod(str));
        temp_file_array.push_back(line_array);
    }
    in_file.close();
    int num_record = temp_file_array.size();
    std::vector<int> permutation_vec(num_record);
    for (int i = 0; i < permutation_vec.size(); i++)
    {
        permutation_vec[i] = i;
    }
    HashChain hash_chain(4, 25);
    for (int i = 0; i < permutation_vec.size(); i++)
    {
        hash_chain.Insert(permutation_vec[i]);
    }
    int permutation_idx = 0;
    for (int i = 0; i < hash_chain.table_.size(); i++)
    {
        for (std::list<int>::const_iterator iter = hash_chain.table_[i].begin(); iter != hash_chain.table_[i].end(); iter++)
        {
            permutation_vec[permutation_idx] = *iter;
            permutation_idx ++;
        }
    }
    printf("After Hashing NumofRecord %d\n", int(permutation_vec.size()));
    file_array.resize(temp_file_array[0].size());
    for (int i = 0; i < file_array.size(); ++i)
    {
        file_array[i].resize(temp_file_array.size());
    }
    for (int i = 0; i < file_array.size(); ++i)
    {
        for (int j = 0; j < file_array[i].size(); ++j)
        {
            file_array[i][j] = temp_file_array[permutation_vec[j]][i];
        } 
    }
}

// void ReadSupportSize(std::string file_name, std::vector<int>& support_size_array)
// {
//     std::ifstream in_file(file_name, std::ios::in);
//     std::string line_str;
//     int count = 1;
//     while (getline(in_file, line_str))
//     {
//         //printf("%d ", count);
//         count++;
//         support_size_array.push_back(stoi(line_str));
//     }
// }

void ReadCumSize(std::string file_name, std::vector<uint32_t>& cum_size)
{
    std::ifstream in_file(file_name, std::ios::in);
    std::string line_str;
    int count = 1;
    while (getline(in_file, line_str))
    {
        std::stringstream ss(line_str);
        uint32_t tmp;
        ss >> tmp;
        cum_size.push_back(tmp);
    }
}

void GenrDistRandNum(std::vector<int> &random_permutation, int prev_num_sample, int num_sample)
{
    int num_element = random_permutation.size();
    for (int i = prev_num_sample; i < num_sample; ++i)
    {
        int target_pos = rand() % (num_element - i);
        std::swap(random_permutation[i], random_permutation[target_pos + i]);
    }
    std::sort(random_permutation.begin() + prev_num_sample, random_permutation.begin() + num_sample);
}

double AccuracyFilter(const int& num_pair, std::vector<int> exact, std::vector<int> est)
{
    std::vector<int> diff;
    sort(exact.begin(), exact.end());
    sort(est.begin(), est.end());
    set_symmetric_difference(exact.begin(), exact.end(), est.begin(), est.end(), inserter(diff, diff.end()));
    return 1.0 - double(diff.size()) / double(num_pair);
}

double Precision(std::vector<int> exact, std::vector<int> est)
{
    if (est.empty())
    {
        if (exact.empty())
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }
    std::vector<int> true_pos;
    sort(exact.begin(), exact.end());
    sort(est.begin(), est.end());
    set_intersection(exact.begin(), exact.end(), est.begin(), est.end(), inserter(true_pos, true_pos.end()));
    return double(true_pos.size()) / double(est.size());
}

double Recall(std::vector<int> exact, std::vector<int> est)
{
    if (exact.empty())
    {
        return 1.0;   
    }
    std::vector<int> true_pos;
    sort(exact.begin(), exact.end());
    sort(est.begin(), est.end());
    set_intersection(exact.begin(), exact.end(), est.begin(), est.end(), inserter(true_pos, true_pos.end()));
    return double(true_pos.size()) / double(exact.size());
}

double F1Measure(const std::vector<int>& exact, const std::vector<int>& est)
{
    double precision = Precision(exact, est), recall = Recall(exact, est);
    if (precision <= 0 && recall <= 0)
    {
        return 0.0;
    }
    else
    {
        return 2.0 * precision * recall / (precision + recall);
    }
}

double NDCG(const std::vector<double>& mut_info, const std::vector<int>& exact, const std::vector<int>& est)
{
    if (exact.size() != est.size())
    {
        printf("Exact: %d and estimated: %d sets are not matched.\n", int(exact.size()), int(est.size()));
        return 0.0;
    }
    double DCG = 0.0, IDCG = 0.0;
    for (int i = 0; i < exact.size(); ++i)
    {
        IDCG += (pow(2.0, mut_info[exact[i]]) - 1) / log2(i + 2.0);
    }
    //printf("IDCG: %f\n", IDCG);
    for (int i = 0; i < est.size(); ++i)
    {
        DCG += (pow(2.0, mut_info[est[i]]) - 1) / log2(i + 2.0);
    }
   // printf("DCG: %f\n", DCG);
    return DCG / IDCG;
}

RelativeError EpsilonDeltaRelativeError(const std::vector<double>& exact_topk_value, const std::vector<double>& approx_topk_value, const double& empirical_threshold)
{
    if (exact_topk_value.size() != approx_topk_value.size())
    {
        printf("Exact %d and approx %d sets are not matched.\n", int(exact_topk_value.size()), int(approx_topk_value.size()));
        return RelativeError(0.0, 0.0);
    }
    double max_relative_error = 0.0, total_relative_error = 0.0;
    for (int i = 0; i < exact_topk_value.size(); i++)
    {
        double relative_error;
        if (exact_topk_value[i] >= empirical_threshold)
        {
            relative_error = abs(exact_topk_value[i] - approx_topk_value[i]) / exact_topk_value[i];
        }
        else
        {
            relative_error = abs(exact_topk_value[i] - approx_topk_value[i]) / empirical_threshold;
        }
        if (relative_error > max_relative_error)
        {
            max_relative_error = relative_error;
        }
        total_relative_error += relative_error;
    }
    return RelativeError(max_relative_error, total_relative_error / exact_topk_value.size());
}

AbsoluteError SelectAbsoluteError(const std::vector<double>& exact_variance, const std::vector<double>& approx_variance)
{
    if (exact_variance.size() != approx_variance.size())
    {
        printf("Exact %d and approx %d sets are not matched.\n", int(exact_variance.size()), int(approx_variance.size()));
        return AbsoluteError(0.0, 0.0);
    }
    double max_absolute_error = 0.0, total_absolute_error = 0.0;
    for (int i = 0; i < exact_variance.size(); i++)
    {
        double absolute_error = abs(exact_variance[i] - approx_variance[i]);
        if (absolute_error > max_absolute_error)
        {
            max_absolute_error = absolute_error;
        }
        total_absolute_error += absolute_error;
    }
    return AbsoluteError(max_absolute_error, total_absolute_error / exact_variance.size());
}

RelativeError TopkRelativeError(const std::vector<double>& exact_topk_value, const std::vector<double>& approx_topk_value)
{
    if (exact_topk_value.size() != approx_topk_value.size())
    {
        printf("Exact %d and approx %d sets are not matched.\n", int(exact_topk_value.size()), int(approx_topk_value.size()));
        return RelativeError(0.0, 0.0);
    }
    double max_relative_error = 0.0, mean_relative_error = 0.0;
    for (int i = 0; i < exact_topk_value.size(); i++)
    {
        if (exact_topk_value[i] > 0.0)
        {
            if (abs(exact_topk_value[i] - approx_topk_value[i]) / exact_topk_value[i] > max_relative_error)
            {
                max_relative_error = abs(exact_topk_value[i] - approx_topk_value[i]) / exact_topk_value[i];
            }
            mean_relative_error += abs(exact_topk_value[i] - approx_topk_value[i]) / exact_topk_value[i];
        }
        else if (exact_topk_value[i] == 0.0)
        {
            if (exact_topk_value[i] - approx_topk_value[i] == 0.0)
            {
                mean_relative_error += 0.0;
            }
            else 
            {
                printf("Try to divide 0!\n");
                printf("i %d value %f\n", i, exact_topk_value[i]);
            }
        }
        else
        {
            printf("Try to divide negative value!\n");
        }
    }
    return RelativeError(max_relative_error, mean_relative_error / exact_topk_value.size());
}

AbsoluteError TopkAbsoluteError(const std::vector<double>& exact_topk_value, const std::vector<double>& approx_topk_value)
{
    if (exact_topk_value.size() != approx_topk_value.size())
    {
        printf("Exact %d and approx %d sets are not matched.\n", int(exact_topk_value.size()), int(approx_topk_value.size()));
        return AbsoluteError(0.0, 0.0);
    }
    double max_absolute_error = 0.0, total_absolute_error = 0.0;
    for (int i = 0; i < exact_topk_value.size(); i++)
    {
        double absolute_error = abs(exact_topk_value[i] - approx_topk_value[i]);
        if (absolute_error > max_absolute_error)
        {
            max_absolute_error = absolute_error;
        }
        total_absolute_error += absolute_error;
    }
    return AbsoluteError(max_absolute_error, total_absolute_error / exact_topk_value.size());
}

RelativeError FilterRelativeError(const std::vector<std::pair<int, double>>& exact_filter_pair, const std::vector<std::pair<int, double>>& approx_filter_pair)
{
    if (exact_filter_pair.size() != approx_filter_pair.size())
    {
        printf("Exact %d and approx %d sets are not matched.\n", int(exact_filter_pair.size()), int(approx_filter_pair.size()));
        return RelativeError(0.0, 0.0);
    }
    int cnt = 0;
    double max_relative_error = 0.0, mean_relative_error = 0.0;
    std::vector<std::pair<int, double>> temp_exact_pair(exact_filter_pair), temp_approx_pair(approx_filter_pair);
    sort(temp_exact_pair.begin(), temp_exact_pair.end(), cmp_idx_value_pair);
    sort(temp_approx_pair.begin(), temp_approx_pair.end(), cmp_idx_value_pair);
    for (int i = 0; i < exact_filter_pair.size(); i++)
    {
        if (temp_exact_pair[i].first == temp_approx_pair[i].first)
        {
            if (abs(temp_exact_pair[i].second - temp_approx_pair[i].second) / temp_exact_pair[i].second > max_relative_error)
            {
                max_relative_error = abs(temp_exact_pair[i].second - temp_approx_pair[i].second) / temp_exact_pair[i].second;
            }
            mean_relative_error += abs(temp_exact_pair[i].second - temp_approx_pair[i].second) / temp_exact_pair[i].second;
            cnt ++;
        }
        else
        {
            continue;
        }
    }
    return RelativeError(max_relative_error, mean_relative_error / cnt);
}

double ExactVarianceScore(const std::vector<double>& data_vec)
{
    int num_element = data_vec.size();
    double variable2_sum = 0.0, variable_sum = 0.0, variance_score = 0.0;
    for (int i = 0; i < num_element; i++)
    {
        variable2_sum += pow(data_vec[i], 2.0);
        variable_sum += data_vec[i];
    }
    variance_score = variable2_sum / num_element - pow(variable_sum / num_element, 2.0);
    return variance_score;
}

double ExactPredVarianceScore(const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec)
{
    int count = 0;
    double sum = 0.0, sum_square = 0.0, variance_score = 0.0;
    for (int i = 0; i < pred_vec.size(); i++)
    {
        if (abs(pred_vec[i] - pred_value) <= 1e-14)
        {
            sum += data_vec[i];
            sum_square += pow(data_vec[i], 2.0);
            count ++;
        }
    }
    variance_score = sum_square / count - pow(sum / count, 2.0);
    return variance_score;
}

double TestPredVarianceScore(const double pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, const int block_size)
{
    int count = 0, idx = 0, rand_start = rand(), test = 0;
    double sum = 0.0, sum_square = 0.0, variance_score = 0.0;
    while (count < block_size)
    {
        // int sample_index = (rand_start + idx) % (pred_vec.size());
        int sample_index = 100 + idx;
        idx ++;
        test += Next() % 4;
        if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
        {
            sum += data_vec[sample_index];
            sum_square += pow(data_vec[sample_index], 2.0);
            count++;           
        }
    }
    variance_score = sum_square / count - pow(sum / count, 2.0);
    printf("Count %d TotalScan %d Variance %f TestSum %d\n", count, idx, variance_score, test);
    return variance_score;
}

double TestPredVarianceScore2(const double pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, const int block_size)
{
    int count = 0, idx = 0, rand_start = rand(), test = 0;
    double sum = 0.0, sum_square = 0.0, variance_score = 0.0;
    // for (int i = 0; i < 10; i++)
    // {
    //     printf("%d\n", nextSmall());
    // }
    while (count < block_size)
    {
        // int sample_index = (rand_start + idx) % (pred_vec.size());
        int sample_index = 100 + idx;
        idx ++;
        test += NextShort() % 4;
        if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
        {
            sum += data_vec[sample_index];
            sum_square += pow(data_vec[sample_index], 2.0);
            count++;           
        }
    }
    variance_score = sum_square / count - pow(sum / count, 2.0);
    printf("Count %d TotalScan %d Variance %f TestSum %d\n", count, idx, variance_score, test);
    return variance_score;
}

// double TestPredVarianceScore3(const double pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, const int block_size)
// {
//     int count = 0, idx = 0, rand_start = rand(), test = 0;
//     double sum = 0.0, sum_square = 0.0, variance_score = 0.0;
//     // for (int i = 0; i < 10; i++)
//     // {
//     //     printf("%d\n", nextSmall());
//     // }
//     while (count < block_size)
//     {
//         // int sample_index = (rand_start + idx) % (pred_vec.size());
//         int sample_index = 100 + idx;
//         idx ++;
//         test += Xoshiro256plus::nextSmall() % 4;
//         if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
//         {
//             sum += data_vec[sample_index];
//             sum_square += pow(data_vec[sample_index], 2.0);
//             count++;           
//         }
//     }
//     variance_score = sum_square / count - pow(sum / count, 2.0);
//     printf("Count %d TotalScan %d Variance %f TestSum %d\n", count, idx, variance_score, test);
//     return variance_score;
// }

double ExactAverageScore(const std::vector<double>& data_vec)
{
    int num_element = data_vec.size();
    double variable_sum = 0.0, average_score = 0.0;
    for (int i = 0; i < num_element; i++)
    {
        variable_sum += data_vec[i];
    }
    average_score = variable_sum / num_element;
    return average_score;
}

void UpdateSingleBlock(const int& block_size, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const int& block_index, const std::vector<double>& data_vec, const std::vector<int>& sample_index_count, const int& start_slot_first, const int& end_slot_last, bool& exceed_final_slot, double& precise)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();

    // int count = 0;
    double sum = 0.0, sum_square = 0.0;
    if (likely(exceed_final_slot == false))
    {
        for (int i = start_slot_first; i < end_slot_last; i++)
        {
            sum += data_vec[i] * sample_index_count[i];
            sum_square += data_vec[i] * data_vec[i] * sample_index_count[i];
            // count ++;
        }
    }
    else
    {
        for (int i = start_slot_first; i < data_vec.size(); i++)
        {
            sum += data_vec[i] * sample_index_count[i];
            sum_square += data_vec[i] * data_vec[i] * sample_index_count[i];
            // count ++;
        }
        for (int i = 0; i < end_slot_last; i++)
        {
            sum += data_vec[i] * sample_index_count[i];
            sum_square += data_vec[i] * data_vec[i] * sample_index_count[i];
            // count ++;
        }
    }
    // printf("Block %d Count %d\n", block_index, count);
    MoM_pair.first[block_index] = sum / block_size;
    MoM_pair.second[block_index] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    precise += duration.count();
}

void UpdateSingleBlock2(const int& block_size, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const int& block_index, const std::vector<double>& data_vec, const std::vector<uint32_t>& sample_vec, const int& start_slot_first, const int& end_slot_last, bool& exceed_final_slot, double& precise)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();

    // int count = 0;
    double sum = 0.0, sum_square = 0.0;
    for (int i = 0; i < sample_vec.size(); i++)
    {
        double sample_value = data_vec[sample_vec[i]];
        sum += sample_value;
        sum_square += sample_value * sample_value;
    }
    // printf("Block %d Count %d\n", block_index, count);
    MoM_pair.first[block_index] = sum / block_size;
    MoM_pair.second[block_index] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    precise += duration.count();
}

void UpdateCumSeq(const int block_size, const int num_element, const std::vector<int>& sample_start_pos, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const std::vector<double>& data_vec, const std::vector<double>& data_slot_idx, const std::vector<uint32_t>& cum_size, double& precise)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    
    for (uint8_t i = 0; i < sample_start_pos.size(); i++)
    {
        int start_slot = data_slot_idx[sample_start_pos[i]], end_slot = data_slot_idx[(sample_start_pos[i] + block_size - 1) % num_element], sample_stop_pos = (sample_start_pos[i] + block_size - 1) % num_element, num_slot = cum_size.size() - 1;
        double sum = 0.0, sum_square = 0.0;
        // printf("start_slot %d end_slot %d sample_start_pos %d sample_stop_pos %d\n", start_slot, end_slot, sample_start_pos[i], sample_stop_pos);
        
        // printf("First slot end + 1 %d\n", cum_size[start_slot].second);
        uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i];
        for (uint8_t j = 0; j < start_slot_sample_size; j++)
        {
            double sample_value = data_vec[cum_size[start_slot] + Next() % (cum_size[start_slot + 1] - cum_size[start_slot])];
            // printf("j %d random index %d\n", j, random_index);
            sum += sample_value;
            sum_square += sample_value * sample_value;
        }
        
        // printf("Medium slot\n");
        if (likely(start_slot < end_slot))
        {
            for (uint32_t j = start_slot + 1; j < end_slot; j++)
            {
                uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                if (slot_size == 1)
                {
                    double sample_value = data_vec[cum_size[j]];
                    sum += sample_value;
                    sum_square += sample_value * sample_value;
                }
                else
                {
                    for (uint8_t k = 0; k < slot_size; k++)
                    {
                        double sample_value = data_vec[cum_size[j] + Next() % slot_size];
                        sum += sample_value;
                        sum_square += sample_value * sample_value;
                    }
                }
            }
        }
        else
        {
            for (uint32_t j = start_slot + 1; j < num_slot; j++)
            {
                uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                if (slot_size == 1)
                {
                    double sample_value = data_vec[cum_size[j]];
                    sum += sample_value;
                    sum_square += sample_value * sample_value;
                }
                else
                {
                    for (uint8_t k = 0; k < slot_size; k++)
                    {
                        double sample_value = data_vec[cum_size[j] + Next() % slot_size];
                        sum += sample_value;
                        sum_square += sample_value * sample_value;
                    }
                }
            }
            for (uint32_t j = 0; j < end_slot; j++)
            {
                uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                if (slot_size == 1)
                {
                    double sample_value = data_vec[cum_size[j]];
                    sum += sample_value;
                    sum_square += sample_value * sample_value;
                }
                else
                {
                    for (uint8_t k = 0; k < slot_size; k++)
                    {
                        double sample_value = data_vec[cum_size[j] + Next() % slot_size];
                        sum += sample_value;
                        sum_square += sample_value * sample_value;
                    }
                }
            }
        }
        
        // printf("Last slot\n");
        uint8_t end_slot_sample_size = sample_stop_pos - cum_size[end_slot] + 1;
        for (uint8_t j = 0; j < end_slot_sample_size; j++)
        {
            double sample_value = data_vec[cum_size[end_slot] + Next() % (cum_size[end_slot + 1] - cum_size[end_slot])];
            sum += sample_value;
            sum_square += sample_value * sample_value;
        }
        MoM_pair.first[i] = sum / block_size;
        MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    precise += duration.count();
}

void UpdateCumSeqPred(const int block_size, const std::vector<uint32_t>& sample_start_pos, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, const std::vector<double>& data_slot_idx, uint32_t& num_scan, const std::vector<uint32_t>& cum_size, double& access_data_time)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    
    int num_element = data_vec.size(), num_slot = cum_size.size() - 1;
    for (uint8_t i = 0; i < sample_start_pos.size(); i++)
    {
        int start_slot = data_slot_idx[sample_start_pos[i]];
        double sum = 0.0, sum_square = 0.0;
        uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i], slot_size = cum_size[start_slot + 1] - cum_size[start_slot];
        num_scan += start_slot_sample_size;
        uint32_t count = 0, current_slot = start_slot;
        for (uint8_t j = 0; j < start_slot_sample_size; j++)
        {
            uint32_t sample_index = cum_size[start_slot] + NextShort() % slot_size;
            if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
            {
                sum += data_vec[sample_index];
                sum_square += data_vec[sample_index] * data_vec[sample_index];
                count ++;
            }
        }
        bool flag = true;
        while (flag)
        {
            current_slot += 1;
            if (unlikely(current_slot >= num_slot))
            {
                current_slot -= num_slot;
            }
            slot_size = cum_size[current_slot + 1] - cum_size[current_slot];
            num_scan += slot_size;
            if (likely(slot_size == 1))
            {
                uint32_t sample_index = cum_size[current_slot];
                if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
                {
                    sum += data_vec[sample_index];
                    sum_square += data_vec[sample_index] * data_vec[sample_index];
                    count ++;
                    if (unlikely(count >= block_size))
                    {
                        flag = false;
                    }
                }
            }
            else
            {
                // uint8_t slot_rand = NextShort();
                for (uint8_t k = 0; k < slot_size; k++)
                {
                    uint32_t sample_index = cum_size[current_slot] + NextShort() % slot_size;
                    // slot_rand >>= 2;
                    if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
                    {
                        sum += data_vec[sample_index];
                        sum_square += data_vec[sample_index] * data_vec[sample_index];
                        count ++;
                        if (unlikely(count >= block_size))
                        {
                            k = slot_size;
                            flag = false;
                        }
                    }
                }
            }
        }
        MoM_pair.first[i] = sum / block_size;
        MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    access_data_time += duration.count();
}

void UpdateCumSeqPredBlock4(const int block_size, const std::vector<uint32_t>& sample_start_pos, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, const std::vector<double>& data_slot_idx, uint32_t& num_scan, const std::vector<uint32_t>& cum_size, double& access_data_time)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    
    int num_element = data_vec.size(), num_slot = cum_size.size() - 1;
    for (uint8_t i = 0; i < sample_start_pos.size(); i++)
    {
        int start_slot = data_slot_idx[sample_start_pos[i]];
        double sum = 0.0, sum_square = 0.0;
        uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i], slot_size = cum_size[start_slot + 1] - cum_size[start_slot];
        num_scan += start_slot_sample_size;
        uint32_t count = 0, current_slot = start_slot, sign;
        for (uint8_t j = 0; j < start_slot_sample_size; j++)
        {
            uint32_t sample_index = cum_size[start_slot] + NextShort() % slot_size;
            sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
            sum += sign * data_vec[sample_index];
            sum_square += sign * data_vec[sample_index] * data_vec[sample_index];
            count += sign;
        }
        // bool flag = true;
        while (count < block_size)
        {
            current_slot += 1;
            if (unlikely(current_slot >= num_slot))
            {
                current_slot -= num_slot;
            }
            slot_size = cum_size[current_slot + 1] - cum_size[current_slot];
            num_scan += slot_size;
            if (likely(slot_size == 1))
            {
                uint32_t sample_index = cum_size[current_slot];
                sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
                sum += sign * data_vec[sample_index];
                sum_square += sign * data_vec[sample_index] * data_vec[sample_index];
                count += sign;
            }
            else
            {
                uint32_t slot_rand = Next();
                for (uint8_t k = 0; k < slot_size; k++)
                {
                    uint32_t sample_index = cum_size[current_slot] + (slot_rand & 0xff) % slot_size;
                    slot_rand >>= 8;
                    sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
                    sum += sign * data_vec[sample_index];
                    sum_square += sign * data_vec[sample_index] * data_vec[sample_index];
                    count += sign;
                    if (unlikely(count >= block_size))
                    {
                        k = slot_size;
                    }
                }
            }
        }
        MoM_pair.first[i] = sum / block_size;
        MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    access_data_time += duration.count();
}

// void UpdateCumSeqPredBlock5(const int block_size, const std::vector<uint32_t>& sample_start_pos, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, const std::vector<double>& data_slot_idx, uint32_t& num_scan, const std::vector<uint32_t>& cum_size, double& access_data_time)
// {
//     std::chrono::high_resolution_clock::time_point start, end; 
//     std::chrono::duration<double> duration;
//     start = std::chrono::high_resolution_clock::now();
    
//     int num_element = data_vec.size(), num_slot = cum_size.size() - 1;
//     for (uint8_t i = 0; i < sample_start_pos.size(); i++)
//     {
//         int start_slot = data_slot_idx[sample_start_pos[i]];
//         double sum = 0.0, sum_square = 0.0;
//         uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i], slot_size = cum_size[start_slot + 1] - cum_size[start_slot];
//         num_scan += start_slot_sample_size;
//         uint32_t count = 0, current_slot = start_slot, sign;
//         for (uint8_t j = 0; j < start_slot_sample_size; j++)
//         {
//             uint32_t sample_index = cum_size[start_slot] + NextShort() % slot_size;
//             sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
//             sum += sign * data_vec[sample_index];
//             sum_square += sign *  data_vec[sample_index] * data_vec[sample_index];
//             count += sign;
//         }
//         while (count < block_size)
//         {
//             current_slot += 1;
//             if (unlikely(current_slot >= num_slot))
//             {
//                 current_slot -= num_slot;
//             }
//             slot_size = cum_size[current_slot + 1] - cum_size[current_slot];
//             num_scan += slot_size;
//             switch (slot_size)
//             {
//                 [[likely]] case 1:
//                 {
//                     uint32_t sample_index = cum_size[current_slot];
//                     sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
//                     sum += sign * data_vec[sample_index];
//                     sum_square += sign * data_vec[sample_index] * data_vec[sample_index];
//                     count += sign;
//                     break;
//                 }
//                 default:
//                 {
//                     uint32_t slot_rand = Next();
//                     for (uint8_t k = 0; k < slot_size; k++)
//                     {
//                         uint32_t sample_index = cum_size[current_slot] + (slot_rand & 0xff) % slot_size;
//                         slot_rand >>= 8;
//                         sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
//                         sum += sign * data_vec[sample_index];
//                         sum_square += sign * data_vec[sample_index] * data_vec[sample_index];
//                         count += sign;
//                         if (unlikely(count >= block_size))
//                         {
//                             k = slot_size;
//                         }
//                     }
//                     break;
//                 }                    
//             }
//         }
//         MoM_pair.first[i] = sum / block_size;
//         MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
//     }
//     end = std::chrono::high_resolution_clock::now();
//     duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//     access_data_time += duration.count();
// }

void UpdateCumSeqPredBlock6(const int block_size, const std::vector<uint32_t>& sample_start_pos, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, const std::vector<double>& data_slot_idx, uint32_t& num_scan, const std::vector<uint32_t>& cum_size, double& access_data_time)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    
    int num_element = data_vec.size(), num_slot = cum_size.size() - 1;
    for (uint8_t i = 0; i < sample_start_pos.size(); i++)
    {
        int start_slot = data_slot_idx[sample_start_pos[i]];
        double sum = 0.0, sum_square = 0.0;
        uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i], slot_size = cum_size[start_slot + 1] - cum_size[start_slot];
        num_scan += start_slot_sample_size;
        uint32_t count = 0, current_slot = start_slot, sign;
        for (uint8_t j = 0; j < start_slot_sample_size; j++)
        {
            uint32_t sample_index = cum_size[start_slot] + NextShort() % slot_size;
            sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
            sum += sign * data_vec[sample_index];
            sum_square += sign *  data_vec[sample_index] * data_vec[sample_index];
            count += sign;
        }
        // bool flag = true;
        while (count < block_size)
        {
            current_slot += 1;
            if (unlikely(current_slot >= num_slot))
            {
                current_slot -= num_slot;
            }
            slot_size = cum_size[current_slot + 1] - cum_size[current_slot];
            num_scan += slot_size;
            if (likely(slot_size == 1))
            {
                uint32_t sample_index = cum_size[current_slot];
                sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
                sum += sign * data_vec[sample_index];
                sum_square += sign * data_vec[sample_index] * data_vec[sample_index];
                count += sign;
            }
            else
            {
                // uint32_t slot_rand = Next();
                for (uint8_t k = 0; k < slot_size; k++)
                {
                    uint32_t sample_index = cum_size[current_slot] + NextShort() % slot_size;
                    // slot_rand >>= 8;
                    sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
                    sum += sign * data_vec[sample_index];
                    sum_square += sign * data_vec[sample_index] * data_vec[sample_index];
                    count += sign;
                    if (unlikely(count >= block_size))
                    {
                        k = slot_size;
                    }
                }
            }
        }
        MoM_pair.first[i] = sum / block_size;
        MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    access_data_time += duration.count();
}

// void UpdateCumSeqPred(const int block_size, const std::vector<uint32_t>& sample_start_pos, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, const std::vector<double>& data_slot_idx, uint32_t& num_scan, const std::vector<uint32_t>& cum_size, double& access_data_time)
// {
//     std::chrono::high_resolution_clock::time_point start, end; 
//     std::chrono::duration<double> duration;
//     start = std::chrono::high_resolution_clock::now();
    
//     int num_element = data_vec.size(), num_slot = cum_size.size() - 1;
//     for (uint8_t i = 0; i < sample_start_pos.size(); i++)
//     {
//         int start_slot = data_slot_idx[sample_start_pos[i]];
//         double sum = 0.0, sum_square = 0.0;
//         // printf("start_slot %d end_slot %d sample_start_pos %d sample_stop_pos %d\n", start_slot, end_slot, sample_start_pos[i], sample_stop_pos);
        
//         // printf("First slot end + 1 %d\n", cum_size[start_slot].second);
//         uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i], slot_size = cum_size[start_slot + 1] - cum_size[start_slot];
//         num_scan += start_slot_sample_size;
//         uint32_t count = 0, current_slot = start_slot;
//         for (uint8_t j = 0; j < start_slot_sample_size; j++)
//         {
//             uint32_t sample_index = cum_size[start_slot] + NextShort() % slot_size;
//             if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
//             {
//                 sum += data_vec[sample_index];
//                 sum_square += data_vec[sample_index] * data_vec[sample_index];
//                 count ++;
//             }
//         }
//         bool flag = true;
//         while (flag)
//         {
//             current_slot += 1;
//             if (current_slot >= num_slot)
//             {
//                 current_slot -= num_slot;
//             }
//             slot_size = cum_size[current_slot + 1] - cum_size[current_slot];
//             num_scan += slot_size;
//             if (slot_size == 1)
//             {
//                 uint32_t sample_index = cum_size[current_slot];
//                 if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
//                 {
//                     sum += data_vec[sample_index];
//                     sum_square += data_vec[sample_index] * data_vec[sample_index];
//                     count ++;
//                     if (count >= block_size)
//                     {
//                         flag = false;
//                     }
//                 }
//             }
//             else
//             {
//                 // uint8_t slot_rand = NextShort();
//                 for (uint8_t k = 0; k < slot_size; k++)
//                 {
//                     uint32_t sample_index = cum_size[current_slot] + NextShort() % slot_size;
//                     // slot_rand >>= 2;
//                     if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
//                     {
//                         sum += data_vec[sample_index];
//                         sum_square += data_vec[sample_index] * data_vec[sample_index];
//                         count ++;
//                         if (count >= block_size)
//                         {
//                             k = slot_size;
//                             flag = false;
//                         }
//                     }
//                 }
//             }
//         }
//         MoM_pair.first[i] = sum / block_size;
//         MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
//     }
//     end = std::chrono::high_resolution_clock::now();
//     duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//     access_data_time += duration.count();
// }

// void UpdateCumSeqPred(const int block_size, const std::vector<uint32_t>& sample_start_pos, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, const std::vector<double>& data_slot_idx, uint32_t& num_scan, const std::vector<uint32_t>& cum_size, double& access_data_time)
// {
//     std::chrono::high_resolution_clock::time_point start, end; 
//     std::chrono::duration<double> duration;
//     start = std::chrono::high_resolution_clock::now();
    
//     int num_element = data_vec.size(), num_slot = cum_size.size() - 1;
//     for (uint8_t i = 0; i < sample_start_pos.size(); i++)
//     {
//         int start_slot = data_slot_idx[sample_start_pos[i]];
//         double sum = 0.0, sum_square = 0.0;

//         uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i], slot_size = cum_size[start_slot + 1] - cum_size[start_slot];
//         num_scan += start_slot_sample_size;
//         uint32_t count = 0, current_slot = start_slot;
//         for (uint8_t j = 0; j < start_slot_sample_size; j++)
//         {
//             uint32_t sample_index = cum_size[start_slot] + NextShort() % slot_size;
//             if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
//             {
//                 sum += data_vec[sample_index];
//                 sum_square += data_vec[sample_index] * data_vec[sample_index];
//                 count ++;
//             }
//         }
//         bool flag = true;
//         while (flag) 
//         {
//             current_slot += 1;
//             if (current_slot >= num_slot) [[unlikely]]
//             {
//                 current_slot -= num_slot;
//             }
//             slot_size = cum_size[current_slot + 1] - cum_size[current_slot];
//             num_scan += slot_size;
//             if (slot_size == 1) [[likely]]
//             {
//                 uint32_t sample_index = cum_size[current_slot];
//                 if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
//                 {
//                     sum += data_vec[sample_index];
//                     sum_square += data_vec[sample_index] * data_vec[sample_index];
//                     count ++;
//                     if (count >= block_size) [[unlikely]]
//                     {
//                         flag = false;
//                     }
//                 }
//             }
//             else
//             {
//                 // uint8_t slot_rand = NextShort();
//                 for (uint8_t k = 0; k < slot_size; k++)
//                 {
//                     uint32_t sample_index = cum_size[current_slot] + NextShort() % slot_size;
//                     // slot_rand >>= 2;
//                     if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
//                     {
//                         sum += data_vec[sample_index];
//                         sum_square += data_vec[sample_index] * data_vec[sample_index];
//                         count ++;
//                         if (count >= block_size) [[unlikely]]
//                         {
//                             k = slot_size;
//                             flag = false;
//                         }
//                     }
//                 }
//             }
//         }
//         MoM_pair.first[i] = sum / block_size;
//         MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
//     }
//     end = std::chrono::high_resolution_clock::now();
//     duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//     access_data_time += duration.count();
// }

void UpdateCumSeqPredBlock2(const int block_size, const std::vector<int>& sample_start_pos, std::vector<uint32_t>& sample_vec, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, const std::vector<double>& data_slot_idx, const std::vector<uint32_t>& cum_size, uint32_t& num_total_scan, double& sample_index_time, double& access_data_time)
{
    sample_vec.resize((block_size + 4) * sample_start_pos.size());
    
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    
    int size = sample_vec.size(), num_element = data_vec.size(), num_slot = cum_size.size() - 1;
    printf("SampleIndex\n");
    for (uint8_t i = 0; i < sample_start_pos.size(); i++)
    {
    //     int count = 0, idx = 0;
    //     while (count < block_size)
    //     {
    //         int sample_index = sample_start_pos[i] + idx;
    //         if (sample_index >= num_element)
    //         {
    //             sample_index -= num_element;
    //         }
    //         idx ++;
    //         if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
    //         {
    //             sample_vec[(block_size + 4) * i + count] = sample_index;
    //             count++;           
    //         }
    //     }
    //     num_total_scan += idx;
    // }
        uint32_t start_slot = data_slot_idx[sample_start_pos[i]], count = 0, current_slot = start_slot;
        double sum = 0.0, sum_square = 0.0;
        // printf("start_slot %d end_slot %d sample_start_pos %d sample_stop_pos %d\n", start_slot, end_slot, sample_start_pos[i], sample_stop_pos);
        
        // printf("First slot end + 1 %d\n", cum_size[start_slot].second);
        uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i], slot_size = cum_size[start_slot + 1] - cum_size[start_slot];
        num_total_scan += start_slot_sample_size;
        for (uint8_t j = 0; j < start_slot_sample_size; j++)
        {
            uint32_t sample_index = cum_size[start_slot] + Next() % slot_size;
            if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
            {
                // if (sample_vec_index >= size)
                // {
                //     printf("index %d\n", sample_vec_index);
                // }
                sample_vec[(block_size + 4) * i + count] = sample_index;
                // sample_vec_index ++;
                count ++;
            }
        }
        // bool flag = true;
        while (count < block_size)
        {
            current_slot += 1;
            if (current_slot >= num_slot)
            {
                current_slot -= num_slot;
            }
            slot_size = cum_size[current_slot + 1] - cum_size[current_slot];
            num_total_scan += slot_size;
            if (slot_size == 1)
            {
                uint32_t sample_index = cum_size[current_slot];
                if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
                {
                    // if (sample_vec_index >= size)
                    // {
                    //     printf("index %d\n", sample_vec_index);
                    // }
                    sample_vec[(block_size + 4) * i + count] = sample_index;
                    // sample_vec_index ++;
                    count ++;
                    // if (unlikely(count >= block_size))
                    // {
                    //     flag = false;
                    // }
                }
            }
            else
            {
                for (uint8_t k = 0; k < slot_size; k++)
                {
                    int sample_index = cum_size[current_slot] + Next() % slot_size;
                    if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
                    {
                        // if (sample_vec_index >= size)
                        // {
                        //     printf("index %d\n", sample_vec_index);
                        // }
                        sample_vec[(block_size + 4) * i + count] = sample_index;
                        // sample_vec_index ++;
                        count ++;
                        // if (unlikely(count >= block_size))
                        // {
                        //     k = slot_size;
                        //     flag = false;
                        // }
                    }
                }
                // std::sort(sample_vec.begin() + (block_size + 4) * i, sample_vec.begin() + (block_size + 4) * i + slot_size);
            }
        }  
        // printf("i %d Count %d sample_index %d\n", i, count, sample_vec_index);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    sample_index_time += duration.count();
    // printf("RealSample %d\n", sample_vec_index);
    start = std::chrono::high_resolution_clock::now();
    for (uint8_t i = 0; i < sample_start_pos.size(); i++)
    {
        double sum = 0.0, sum_square = 0.0;
        for (int j = 0; j < block_size; j++)
        {
            double sample_value = data_vec[sample_vec[i * block_size + j]];
            sum += sample_value;
            sum_square += sample_value * sample_value;
        }
        MoM_pair.first[i] = sum / block_size;
        MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
    }
    // printf("FinishOneFeature\n");
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    access_data_time += duration.count();
}

void UpdateCumSeqPredBlock3(const int block_size, const std::vector<int>& sample_start_pos, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, const std::vector<double>& data_slot_idx, const std::vector<uint32_t>& cum_size, double& access_data_time)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    
    int num_element = data_vec.size(), num_slot = cum_size.size() - 1;
    for (uint8_t i = 0; i < sample_start_pos.size(); i++)
    {
        int start_slot = data_slot_idx[sample_start_pos[i]];
        double sum = 0.0, sum_square = 0.0;
        // printf("start_slot %d end_slot %d sample_start_pos %d sample_stop_pos %d\n", start_slot, end_slot, sample_start_pos[i], sample_stop_pos);
        
        // printf("First slot end + 1 %d\n", cum_size[start_slot].second);
        uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i];

        uint32_t count = 0, current_slot = start_slot;
        for (uint8_t j = 0; j < start_slot_sample_size; j++)
        {
            uint32_t sample_index = cum_size[start_slot] + Next() % (cum_size[start_slot + 1] - cum_size[start_slot]);
            // uint32_t sample_index = cum_size[start_slot] + j;
            if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
            {
                sum += data_vec[sample_index];
                sum_square += data_vec[sample_index] * data_vec[sample_index];
                count ++;
            }
        }
        // bool flag = true;
        while (count < block_size)
        {
            current_slot = (current_slot + 1) % num_slot;
            uint8_t slot_size = cum_size[current_slot + 1] - cum_size[current_slot];
            if (slot_size == 1)
            {
                uint32_t sample_index = cum_size[current_slot];
                if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
                {
                    sum += data_vec[sample_index];
                    sum_square += data_vec[sample_index] * data_vec[sample_index];
                    count ++;
                    // if (count >= block_size)
                    // {
                    //     flag = false;
                    // }
                }
            }
            else
            {
                for (uint8_t j = 0; j < slot_size; j++)
                {
                    uint32_t sample_index = cum_size[current_slot] + Next() % slot_size;
                    // uint32_t sample_index = cum_size[current_slot] + j;
                    if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
                    {
                        sum += data_vec[sample_index];
                        sum_square += data_vec[sample_index] * data_vec[sample_index];
                        count ++;
                        // if (count >= block_size)
                        // {
                        //     j = slot_size;
                        //     flag = false;
                        // }
                    }
                }
            }
        }
        MoM_pair.first[i] = sum / count;
        MoM_pair.second[i] = (sum_square / count - pow(sum / count, 2.0)) * count / (count - 1.0);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    access_data_time += duration.count();
}

void MomSamplePredIndex(const int block_size, const std::vector<u_int32_t>& sample_start_pos, std::vector<uint32_t>& sample_vec, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_slot_idx, const std::vector<uint32_t>& cum_size, double& medium)
{
    sample_vec.resize((block_size + 4) * sample_start_pos.size());
    
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    
    int size = sample_vec.size(), num_element = pred_vec.size(), num_slot = cum_size.size() - 1;
    for (uint8_t i = 0; i < sample_start_pos.size(); i++)
    {
        uint32_t start_slot = data_slot_idx[sample_start_pos[i]], count = 0, current_slot = start_slot;

        uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i];
        
        for (uint8_t j = 0; j < start_slot_sample_size; j++)
        {
            uint32_t sample_index = cum_size[start_slot] + Next() % (cum_size[start_slot + 1] - cum_size[start_slot]);
            if (likely(pred_vec[sample_index] == pred_value))
            {
                sample_vec[(block_size + 4) * i + count] = sample_index;
                count ++;
            }
        }

        while (count < block_size)
        {
            current_slot = (current_slot + 1) % num_slot;
            uint8_t slot_size = cum_size[current_slot + 1] - cum_size[current_slot];
            if (slot_size == 1)
            {
                uint32_t sample_index = cum_size[current_slot];
                if (likely(pred_vec[sample_index] == pred_value))
                {
                    sample_vec[(block_size + 4) * i + count] = sample_index;
                    count ++;
                }
            }
            else
            {
                for (uint8_t k = 0; k < slot_size; k++)
                {
                    int sample_index = cum_size[current_slot] + Next() % slot_size;
                    if (likely(pred_vec[sample_index] == pred_value))
                    {
                        sample_vec[(block_size + 4) * i + count] = sample_index;
                        count ++;
                    }
                }
            }
        }  
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    medium += duration.count();
}

void MomSamplePredIndexBlock2(const int block_size, const std::vector<u_int32_t>& sample_start_pos, std::vector<uint32_t>& sample_vec, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_slot_idx, uint32_t& num_scan, const std::vector<uint32_t>& cum_size, double& medium)
{
    sample_vec.resize((block_size + 4) * sample_start_pos.size());
    
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    
    int num_element = pred_vec.size(), num_slot = cum_size.size() - 1;
    bool sign;

    for (uint8_t i = 0; i < sample_start_pos.size(); i++)
    {
        uint32_t start_slot = data_slot_idx[sample_start_pos[i]], count = 0, current_slot = start_slot;

        uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i], slot_size = cum_size[start_slot + 1] - cum_size[start_slot];
        num_scan += start_slot_sample_size;
        for (uint8_t j = 0; j < start_slot_sample_size; j++)
        {
            uint32_t sample_index = cum_size[start_slot] + Next() % slot_size;
            sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
            sample_vec[(block_size + 4) * i + count + !sign] = sample_index;
            count += sign;
        }

        while (count < block_size)
        {
            current_slot += 1;
            if (unlikely(current_slot >= num_slot))
            {
                current_slot -= num_slot;
            }
            slot_size = cum_size[current_slot + 1] - cum_size[current_slot];
            num_scan += slot_size;
            if (likely(slot_size == 1))
            {
                uint32_t sample_index = cum_size[current_slot];
                sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
                sample_vec[(block_size + 4) * i + count + !sign] = sample_index;
                count += sign;
            }
            else
            {
                for (uint8_t k = 0; k < slot_size; k++)
                {
                    uint32_t sample_index = cum_size[current_slot] + Next() % slot_size;
                    sign = my_sign(abs(pred_vec[sample_index] - pred_value) - 1e-12);
                    sample_vec[(block_size + 4) * i + count + !sign] = sample_index;
                    count += sign;
                }
            }
        }  
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    medium += duration.count();
}

void SampleBlockRandomNum(std::vector<uint32_t>& sample_vec, int num_element, int num_var, int block_size, const std::vector<int>& sample_start_pos, const std::vector<double>& data_slot_idx, const std::vector<uint32_t>& cum_size, double& precise)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    sample_vec.resize(block_size * num_var);
    int num_slot = cum_size.size() - 1;
    for (uint8_t i = 0; i < sample_start_pos.size(); i++)
    {
        int start_slot = data_slot_idx[sample_start_pos[i]], end_slot = data_slot_idx[(sample_start_pos[i] + block_size - 1) % num_element], sample_stop_pos = (sample_start_pos[i] + block_size - 1) % num_element, sample_vector_index = 0;
        double sum = 0.0, sum_square = 0.0;
        // printf("start_slot %d end_slot %d sample_start_pos %d sample_stop_pos %d\n", start_slot, end_slot, sample_start_pos[i], sample_stop_pos);
        
        // printf("First slot end + 1 %d\n", cum_size[start_slot].second);
        uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i];
        for (uint8_t j = 0; j < start_slot_sample_size; j++)
        {
            sample_vec[block_size * i + sample_vector_index] = cum_size[start_slot] + Next() % (cum_size[start_slot + 1] - cum_size[start_slot]);
            // printf("j %d random index %d\n", j, random_index);
            sample_vector_index ++;
        }
        
        // printf("Medium slot\n");
        if (likely(start_slot < end_slot))
        {
            for (uint32_t j = start_slot + 1; j < end_slot; j++)
            {
                uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                if (slot_size == 1)
                {
                    sample_vec[block_size * i + sample_vector_index] = cum_size[j];
                    sample_vector_index ++;
                }
                else
                {
                    for (uint8_t k = 0; k < slot_size; k++)
                    {
                        sample_vec[block_size * i + sample_vector_index] = cum_size[j] + Next() % slot_size;
                        sample_vector_index ++;
                    }
                }
            }
        }
        else
        {
            for (uint32_t j = start_slot + 1; j < num_slot; j++)
            {
                uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                if (slot_size == 1)
                {
                    sample_vec[block_size * i + sample_vector_index] = cum_size[j];
                    sample_vector_index ++;
                }
                else
                {
                    for (uint8_t k = 0; k < slot_size; k++)
                    {
                        sample_vec[block_size * i + sample_vector_index] = cum_size[j] + Next() % slot_size;
                        sample_vector_index ++;
                    }
                }
            }
            for (uint32_t j = 0; j < end_slot; j++)
            {
                uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                if (slot_size == 1)
                {
                    sample_vec[block_size * i + sample_vector_index] = cum_size[j];
                    sample_vector_index ++;
                }
                else
                {
                    for (uint8_t k = 0; k < slot_size; k++)
                    {
                        sample_vec[block_size * i + sample_vector_index] = cum_size[j] + Next() % slot_size;
                        sample_vector_index ++;
                    }
                }
            }
        }
        
        // printf("Last slot\n");
        uint8_t end_slot_sample_size = sample_stop_pos - cum_size[end_slot] + 1;
        for (uint8_t j = 0; j < end_slot_sample_size; j++)
        {
            sample_vec[block_size * i + sample_vector_index] = cum_size[end_slot] + Next() % (cum_size[end_slot + 1] - cum_size[end_slot]);
            sample_vector_index ++;
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    precise += duration.count();
}

void UpdateCumRandom(const int& sample_size, std::pair<double, double>& cum_pair, const std::vector<double>& data_vec, double& precise)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    
    cum_pair = std::pair<double, double>(0.0, 0.0);
    for (int i = 0; i < sample_size; i++)
    {
        int sample_index = rand() % data_vec.size();
        cum_pair.first += data_vec[sample_index];
        cum_pair.second += pow(data_vec[sample_index], 2.0);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    precise += duration.count();
}

void UpdateCumRandPred(const int& sample_size, std::pair<double, double>& cum_pair, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, uint32_t& num_scan, double& precise)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    
    uint32_t count = 0;
    cum_pair = std::pair<double, double>(0.0, 0.0);
    // int sample = rand() % data_vec.size();
    // printf("size %lu sample %d\n", data_vec.size(), sample);
    while(count < sample_size)
    {
        int sample_index = rand() % data_vec.size();
        num_scan ++;
        if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
        {
            cum_pair.first += data_vec[sample_index];
            cum_pair.second += pow(data_vec[sample_index], 2.0);
            count ++;
        }
    }

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    precise += duration.count();
}

void UpdateCumRandPredBlock2(const int& sample_size, std::vector<int>& sample_vec, std::pair<double, double>& cum_pair, const double& pred_value, const std::vector<double>& pred_vec, const std::vector<double>& data_vec, double& sample_index_time, double& access_data_time)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
    uint32_t count = 0;
    cum_pair = std::pair<double, double>(0.0, 0.0);
    sample_vec.resize(0);
    while(count < sample_size)
    {
        int sample_index = rand() % data_vec.size();
        if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
        {
            sample_vec.push_back(sample_index);
            count ++;
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    sample_index_time += duration.count();

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < sample_vec.size(); i++)
    {
        cum_pair.first += data_vec[sample_vec[i]];
        cum_pair.second += pow(data_vec[sample_vec[i]], 2.0);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    access_data_time += duration.count();
}

void ApproxSamplePredIndex(const int& sample_size, std::vector<uint32_t>& sample_vec, const double& pred_value, const std::vector<double>& pred_vec, uint32_t& num_scan, double& medium)
{
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();
     
    printf("RandPredBlock2\n");
    sample_vec.resize(sample_size);
    uint32_t count = 0;
    while(count < sample_size)
    {
        int sample_index = rand() % pred_vec.size();
        num_scan ++;
        if (abs(pred_vec[sample_index] - pred_value) <= 1e-12)
        {
            sample_vec[count] = sample_index;
            count ++;
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    medium += duration.count();
}

// void UpdateCumSeqIndex(const int& block_size, const int& num_element, const std::vector<int>& sample_start_pos, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const std::vector<double>& data_vec, const std::vector<double>& data_slot_idx, const std::vector<uint32_t>& cum_size, double& precise)
// {
//     std::chrono::high_resolution_clock::time_point start, end; 
//     std::chrono::duration<double> duration;
//     start = std::chrono::high_resolution_clock::now();
//     // printf("data_slot_index\n");
//     // for (int i = 0; i < 100; i++)
//     // {
//     //     printf("i %d %lf\n", i, data_slot_idx[i]);
//     // }
    
//     for (int i = 0; i < sample_start_pos.size(); i++)
//     {
//         int start_slot = data_slot_idx[sample_start_pos[i]], end_slot = data_slot_idx[(sample_start_pos[i] + block_size - 1) % num_element], sample_stop_pos = (sample_start_pos[i] + block_size - 1) % num_element, num_slot = cum_size.size() - 1, count = 0;
//         double sum = 0.0, sum_square = 0.0;
//         std::vector<int32_t> sample_idx(block_size);
//         // printf("start_slot %d end_slot %d sample_start_pos %d sample_stop_pos %d\n", start_slot, end_slot, sample_start_pos[i], sample_stop_pos);
        
//         // printf("First slot end + 1 %d\n", cum_size[start_slot].second);
//         uint8_t first_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i];
//         for (uint8_t j = 0; j < first_slot_sample_size; j++)
//         {
//             sample_idx[count] = cum_size[start_slot] + Next() % (cum_size[start_slot + 1] - cum_size[start_slot]);
//             count ++;
//             // printf("j %d random index %d\n", j, random_index);
//             // sum += sample_value;
//             // sum_square += sample_value * sample_value;
//         }
        
//         // printf("Medium slot\n");
//         for (int j = (start_slot + 1) % num_slot; j < end_slot % num_slot; j++)
//         {
//             uint8_t slot_size = cum_size[j + 1] - cum_size[j];
//             if (slot_size == 1)
//             {
//                 sample_idx[count] = cum_size[j];
//                 // double sample_value = data_vec[cum_size[j]];
//                 // sum += sample_value;
//                 // sum_square += sample_value * sample_value;
//                 count ++;
//             }
//             else
//             {
//                 for (uint8_t k = 0; k < slot_size; k++)
//                 {
//                     sample_idx[count] = cum_size[j] + Next() % slot_size;
//                     count ++;
//                     // double sample_value = data_vec[cum_size[j] + Next() % slot_size];
//                     // sum += sample_value;
//                     // sum_square += sample_value * sample_value;
//                 }
//             }
//         }
//         // printf("Last slot\n");
//         for (int j = 0; j < sample_stop_pos - cum_size[end_slot] + 1; j++)
//         {
//             sample_idx[count] = cum_size[end_slot] + Next() % (cum_size[end_slot + 1] - cum_size[end_slot]);
//             count ++;
//             // double sample_value = data_vec[cum_size[end_slot] + Next() % (cum_size[end_slot + 1] - cum_size[end_slot])];
//             // sum += sample_value;
//             // sum_square += sample_value * sample_value;
//         }
        
//         for (uint j = 0; j < sample_idx.size(); j++)
//         {
//             double sample_value = data_vec[sample_idx[j]];
//             sum += sample_value;
//             sum_square += sample_value * sample_value;
//         }

//         // for (uint8_t j = 0; j < first_slot_sample_size; j++)
//         // {
//         //     double sample_value = data_vec[cum_size[start_slot] + Next() % (cum_size[start_slot + 1] - cum_size[start_slot])];
//         //     // printf("j %d random index %d\n", j, random_index);
//         //     sum += sample_value;
//         //     sum_square += sample_value * sample_value;
//         // }
        
//         // // printf("Medium slot\n");
//         // for (int j = (start_slot + 1) % num_slot; j < end_slot % num_slot; j++)
//         // {
//         //     uint8_t slot_size = cum_size[j + 1] - cum_size[j];
//         //     if (slot_size == 1)
//         //     {
//         //         double sample_value = data_vec[cum_size[j]];
//         //         sum += sample_value;
//         //         sum_square += sample_value * sample_value;
//         //     }
//         //     else
//         //     {
//         //         for (uint8_t k = 0; k < slot_size; k++)
//         //         {
//         //             double sample_value = data_vec[cum_size[j] + Next() % slot_size];
//         //             sum += sample_value;
//         //             sum_square += sample_value * sample_value;
//         //         }
//         //     }
//         // }
//         // // printf("Last slot\n");
//         // for (int j = 0; j < sample_stop_pos - cum_size[end_slot] + 1; j++)
//         // {
//         //     double sample_value = data_vec[cum_size[end_slot] + Next() % (cum_size[end_slot + 1] - cum_size[end_slot])];
//         //     sum += sample_value;
//         //     sum_square += sample_value * sample_value;
//         // }
//         MoM_pair.first[i] = sum / block_size;
//         MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
//     }
//     end = std::chrono::high_resolution_clock::now();
//     duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//     precise += duration.count();
// }

// void UpdateCumSeqIndex2(const int& block_size, const int& num_element, const std::vector<int>& sample_start_pos, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const std::vector<double>& data_vec, const std::vector<double>& data_slot_idx, const std::vector<uint32_t>& cum_size, double& precise, std::vector<uint32_t>& count_idx)
// {
//     std::chrono::high_resolution_clock::time_point start, end; 
//     std::chrono::duration<double> duration;
//     start = std::chrono::high_resolution_clock::now();
//     // printf("data_slot_index\n");
//     // for (int i = 0; i < 100; i++)
//     // {
//     //     printf("i %d %lf\n", i, data_slot_idx[i]);
//     // }
    
//     for (int i = 0; i < sample_start_pos.size(); i++)
//     {
//         int start_slot = data_slot_idx[sample_start_pos[i]], end_slot = data_slot_idx[(sample_start_pos[i] + block_size - 1) % num_element], start_slot_first = cum_size[start_slot], end_slot_last = cum_size[end_slot + 1], sample_stop_pos = (sample_start_pos[i] + block_size - 1) % num_element, num_slot = cum_size.size() - 1;
//         if (start_slot_first < end_slot_last)
//         {
//             std::fill(count_idx.begin() + start_slot_first, count_idx.begin() + end_slot_last, 0);
//         }
//         else
//         {
//             std::fill(count_idx.begin() + start_slot_first, count_idx.end(), 0);
//             std::fill(count_idx.begin(), count_idx.begin() + end_slot_last, 0);
//         }
//         double sum = 0.0, sum_square = 0.0;
//         // std::vector<int32_t> sample_idx(block_size);
//         // printf("start_slot %d end_slot %d sample_start_pos %d sample_stop_pos %d\n", start_slot, end_slot, sample_start_pos[i], sample_stop_pos);
        
//         // printf("First slot end + 1 %d\n", cum_size[start_slot].second);
//         uint8_t first_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i];
//         for (uint8_t j = 0; j < first_slot_sample_size; j++)
//         {
//             count_idx[cum_size[start_slot] + Next() % (cum_size[start_slot + 1] - cum_size[start_slot])] += 1;
//             // printf("j %d random index %d\n", j, random_index);
//             // sum += sample_value;
//             // sum_square += sample_value * sample_value;
//         }
        
//         // printf("Medium slot\n");
//         for (int j = (start_slot + 1) % num_slot; j < end_slot % num_slot; j++)
//         {
//             uint8_t slot_size = cum_size[j + 1] - cum_size[j];
//             if (slot_size == 1)
//             {
//                 count_idx[cum_size[j]] += 1;
//                 // double sample_value = data_vec[cum_size[j]];
//                 // sum += sample_value;
//                 // sum_square += sample_value * sample_value;
//                 // count ++;
//             }
//             else
//             {
//                 for (uint8_t k = 0; k < slot_size; k++)
//                 {
//                     count_idx[cum_size[j] + Next() % slot_size] += 1;
//                     // count ++;
//                     // double sample_value = data_vec[cum_size[j] + Next() % slot_size];
//                     // sum += sample_value;
//                     // sum_square += sample_value * sample_value;
//                 }
//             }
//         }
//         // printf("Last slot\n");
//         for (int j = 0; j < sample_stop_pos - cum_size[end_slot] + 1; j++)
//         {
//             count_idx[cum_size[end_slot] + Next() % (cum_size[end_slot + 1] - cum_size[end_slot])] += 1;
//             // count ++;
//             // double sample_value = data_vec[cum_size[end_slot] + Next() % (cum_size[end_slot + 1] - cum_size[end_slot])];
//             // sum += sample_value;
//             // sum_square += sample_value * sample_value;
//         }
        
//         if (start_slot < end_slot)
//         {
//             for (uint j = start_slot_first; j < end_slot_last; j++)
//             {
//                 sum += data_vec[j] * count_idx[j];
//                 sum_square += data_vec[j] * data_vec[j] * count_idx[j];
//             }
//         }
//         else
//         {
//             for (uint j = start_slot_first; j < num_element; j++)
//             {
//                 sum += data_vec[j] * count_idx[j];
//                 sum_square += data_vec[j] * data_vec[j] * count_idx[j];
//             }
//             for (uint j = 0; j < end_slot_last; j++)
//             {
//                 sum += data_vec[j] * count_idx[j];
//                 sum_square += data_vec[j] * data_vec[j] * count_idx[j];
//             }
//         }

//         // for (uint8_t j = 0; j < first_slot_sample_size; j++)
//         // {
//         //     double sample_value = data_vec[cum_size[start_slot] + next() % (cum_size[start_slot + 1] - cum_size[start_slot])];
//         //     // printf("j %d random index %d\n", j, random_index);
//         //     sum += sample_value;
//         //     sum_square += sample_value * sample_value;
//         // }
        
//         // // printf("Medium slot\n");
//         // for (int j = (start_slot + 1) % num_slot; j < end_slot % num_slot; j++)
//         // {
//         //     uint8_t slot_size = cum_size[j + 1] - cum_size[j];
//         //     if (slot_size == 1)
//         //     {
//         //         double sample_value = data_vec[cum_size[j]];
//         //         sum += sample_value;
//         //         sum_square += sample_value * sample_value;
//         //     }
//         //     else
//         //     {
//         //         for (uint8_t k = 0; k < slot_size; k++)
//         //         {
//         //             double sample_value = data_vec[cum_size[j] + next() % slot_size];
//         //             sum += sample_value;
//         //             sum_square += sample_value * sample_value;
//         //         }
//         //     }
//         // }
//         // // printf("Last slot\n");
//         // for (int j = 0; j < sample_stop_pos - cum_size[end_slot] + 1; j++)
//         // {
//         //     double sample_value = data_vec[cum_size[end_slot] + next() % (cum_size[end_slot + 1] - cum_size[end_slot])];
//         //     sum += sample_value;
//         //     sum_square += sample_value * sample_value;
//         // }
//         MoM_pair.first[i] = sum / block_size;
//         MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
//     }
//     end = std::chrono::high_resolution_clock::now();
//     duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//     precise += duration.count();
// }

// void UpdateOneBlockOld(const int& block_size, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const int& block_index, const std::vector<double>& data_vec, const std::vector<int>& sample_index_count, const int& start_slot_first, const int& end_slot_last, bool& exceed_final_slot, double& precise)
// {
//     std::chrono::high_resolution_clock::time_point start, end; 
//     std::chrono::duration<double> duration;
//     start = std::chrono::high_resolution_clock::now();

//     double sum = 0.0, sum_square = 0.0;
//     if (likely(exceed_final_slot == false))
//     {
//         for (int i = start_slot_first; i < start_slot_first + block_size; i++)
//         {
//             // if (likely(sample_index_count[i] != 0))
//             // {
//             //     sum += data_vec[i] * sample_index_count[i];
//             //     sum_square += pow(data_vec[i], 2.0) * sample_index_count[i];
//             // }
//             sum += data_vec[i];
//             sum_square += pow(data_vec[i], 2.0);
//         }
//     }
//     else
//     {
//         for (int i = start_slot_first; i < num_element; i++)
//         {
//             sum += data_vec[i];
//             sum_square += pow(data_vec[i], 2.0);
//         }
//         for (int i = 0; i < block_size - (num_element - start_slot_first); i++)
//         {
//             sum += data_vec[i];
//             sum_square += pow(data_vec[i], 2.0);
//         }
//     }
//     MoM_pair.first[block_index] = sum / block_size;
//     MoM_pair.second[block_index] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);

//     end = std::chrono::high_resolution_clock::now();
//     duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//     precise += duration.count();
// }

void ApproxVarianceTopk(const std::string mode, const int topk, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, std::vector<int>& result_idx, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    long count = 0; 
    SampleSizePair sample_size_pair(0, initial_sample_size);

    double pf0 = pf * 1.0 / num_feature / ceil(log2(num_element * 1.0 / initial_sample_size)), sample_index_time = 0.0, access_data_time = 0.0, sort_time = 0.0;
    printf("Element %d Feature %d pf0 %1.15f\n", num_element, num_feature, pf0);
    std::vector<std::pair<double, double>> cumulate_pair(num_feature, std::pair<double, double>(0.0, 0.0));
    std::vector<Feature> feature_vec;
    std::vector<uint32_t> sample_vec;
    
    for (int i = 0; i < num_feature; i++)
    {
        feature_vec.push_back(Feature(i, 0.0));
    }
    result_idx.clear();
    result_value.clear();
    while (!feature_vec.empty() && sample_size_pair.current_ <= num_element)
    {
        std::chrono::high_resolution_clock::time_point start, end; 
        std::chrono::duration<double> duration;
        start = std::chrono::high_resolution_clock::now();

        sample_vec.clear();
        sample_vec.resize(sample_size_pair.current_);
        for (int i = 0; i < sample_vec.size(); i++)
        {
            sample_vec[i] = rand() % num_element;
        }

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        sample_index_time += duration.count();

        for (int i = 0; i < feature_vec.size(); i++)
        {
            int feature_idx = feature_vec[i].idx_;
            // if (mode == "random_colt")
            // {
            //     printf("idx %d ", feature_idx);
            // }
            VarianceEst variance(sample_vec, mode, sample_size_pair, cumulate_pair[feature_idx], data_mat[feature_idx], pf0, access_data_time);
            feature_vec[i].Update(variance.low_, variance.up_, variance.sample_variance_);
        }
        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_low);
        double kth_low = feature_vec[topk - 1].low_;
        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
        sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
        double min_ratio = std::min_element(feature_vec.begin(), feature_vec.begin() + topk, cmp_ratio)->ratio_;
        printf("SampleSize %d ExistedFeature %d\n", sample_size_pair.current_, int(feature_vec.size()));
        for (int i = 0; i < topk; i++)
        {
            printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        } 
        count += sample_size_pair.current_ * feature_vec.size();
        if (min_ratio >= (1 - error_bound))
        {
            printf("ApproximateVarianceTopk has been found.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }
            printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f\n", sample_size_pair.current_, count, sample_index_time, access_data_time);
            break;
        }
        else if (sample_size_pair.current_ >= num_element / 2 && sample_size_pair.current_ < num_element)
        {
            sample_size_pair.previous_ = sample_size_pair.current_;
            sample_size_pair.current_ = num_element; 
            
            std::chrono::high_resolution_clock::time_point start, end; 
            std::chrono::duration<double> duration;
            start = std::chrono::high_resolution_clock::now();
            
            for (int i = 0; i < feature_vec.size(); i++)
            {
                int feature_idx = feature_vec[i].idx_;
                double exact_variance = ExactVarianceScore(data_mat[feature_idx]);
                feature_vec[i].Update(exact_variance, exact_variance, exact_variance);
            }
            std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
            sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
            printf("ApproximateVarianceTopk has been found when elements are used up.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            sort_time += duration.count();

            printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f SortTime %f\n", sample_size_pair.current_, count, sample_index_time, access_data_time, sort_time);
            break;
        }
        else
        {
            sample_size_pair.previous_ = sample_size_pair.current_;
            sample_size_pair.current_ *= 2; 
        }
        
        // if (sample_size_pair.previous_ >= num_element)
        // {
        //     printf("ApproximateVarianceTopk has been found when elements are used up.\n");
        //     for (int i = 0; i < topk; i++)
        //     {
        //         result_idx.push_back(feature_vec[i].idx_);
        //         result_value.push_back(feature_vec[i].avg_);
        //         printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        //     }
        //     break;
        // }
        std::vector<Feature>::const_iterator iter = feature_vec.begin();
        while (iter != feature_vec.end()) 
        {
            if (iter->up_ < kth_low) 
            {
                iter = feature_vec.erase(iter);
            }
            else 
            {
                ++iter;
            }
        }
    } 
    
}

// void ApproxVarianceTopkPred(const std::string mode, const int topk, const PredPair& pred, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, std::vector<int>& result_idx, const double pf, const double error_bound, std::vector<double>& result_value)
// {
//     const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
//     long count = 0; 
//     SampleSizePair sample_size_pair(0, initial_sample_size);

//     double pf0 = pf * 1.0 / num_feature / ceil(log2(num_element * 1.0 / initial_sample_size)), access_data_time = 0.0, sort_time = 0.0, sample_index_time = 0.0;
//     printf("Element %d Feature %d pf0 %1.15f\n", num_element, num_feature, pf0);
//     std::vector<std::pair<double, double>> cumulate_pair(num_feature, std::pair<double, double>(0.0, 0.0));
//     std::vector<Feature> feature_vec;
//     std::vector<uint32_t> sample_vec;
    
//     for (int i = 0; i < num_feature; i++)
//     {
//         if (i != pred.col_)
//         {
//             feature_vec.push_back(Feature(i, 0.0));
//         }
//     }
//     result_idx.clear();
//     result_value.clear();
//     while (!feature_vec.empty() && sample_size_pair.current_ <= num_element)
//     {
//         std::chrono::high_resolution_clock::time_point start, end; 
//         std::chrono::duration<double> duration;
//         start = std::chrono::high_resolution_clock::now();

//         sample_vec.resize(sample_size_pair.current_);
//         for (int i = 0; i < sample_vec.size(); i++)
//         {
//             sample_vec[i] = rand() % num_element;
//         }

//         ApproxSamplePredIndex(sample_size_pair.current_, sample_vec, pred.value_, data_mat[pred.col_], sample_index_time);

//         end = std::chrono::high_resolution_clock::now();
//         duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//         sort_time += duration.count();

//         for (int i = 0; i < feature_vec.size(); i++)
//         {
//             int feature_idx = feature_vec[i].idx_;
//             // if (mode == "random_colt")
//             // {
//             //     printf("idx %d ", feature_idx);
//             // }
//             VarianceEst variance(sample_vec, mode, sample_size_pair, cumulate_pair[feature_idx], data_mat[feature_idx], pf0, access_data_time);
//             feature_vec[i].Update(variance.low_, variance.up_, variance.sample_variance_);
//         }
//         std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_low);
//         double kth_low = feature_vec[topk - 1].low_;
//         std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
//         sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
//         double min_ratio = std::min_element(feature_vec.begin(), feature_vec.begin() + topk, cmp_ratio)->ratio_;
//         printf("SampleSize %d ExistedFeature %d\n", sample_size_pair.current_, int(feature_vec.size()));
//         for (int i = 0; i < topk; i++)
//         {
//             printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
//         } 
//         count += sample_size_pair.current_ * feature_vec.size();
//         if (min_ratio >= (1 - error_bound))
//         {
//             printf("ApproximateVarianceTopk has been found.\n");
//             for (int i = 0; i < topk; i++)
//             {
//                 result_idx.push_back(feature_vec[i].idx_);
//                 result_value.push_back(feature_vec[i].avg_);
//                 printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
//             }
//             printf("SampleSize %d Count %ld SampleDataTime %f AccessDataTime %f\n", sample_size_pair.current_, count, sample_index_time, access_data_time);
//             break;
//         }
//         else if (sample_size_pair.current_ >= num_element / 2 && sample_size_pair.current_ < num_element)
//         {
//             sample_size_pair.previous_ = sample_size_pair.current_;
//             sample_size_pair.current_ = num_element; 
            
//             std::chrono::high_resolution_clock::time_point start, end; 
//             std::chrono::duration<double> duration;
//             start = std::chrono::high_resolution_clock::now();
            
//             for (int i = 0; i < feature_vec.size(); i++)
//             {
//                 int feature_idx = feature_vec[i].idx_;
//                 double exact_variance = ExactVarianceScore(data_mat[feature_idx]);
//                 feature_vec[i].Update(exact_variance, exact_variance, exact_variance);
//             }
//             std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
//             sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
//             printf("ApproximateVarianceTopk has been found when elements are used up.\n");
//             for (int i = 0; i < topk; i++)
//             {
//                 result_idx.push_back(feature_vec[i].idx_);
//                 result_value.push_back(feature_vec[i].avg_);
//                 printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
//             }
//             end = std::chrono::high_resolution_clock::now();
//             duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//             sort_time += duration.count();

//             printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f\n", sample_size_pair.current_, count, sample_index_time, access_data_time);
//             break;
//         }
//         else
//         {
//             sample_size_pair.previous_ = sample_size_pair.current_;
//             sample_size_pair.current_ *= 2; 
//         }
        
//         // if (sample_size_pair.previous_ >= num_element)
//         // {
//         //     printf("ApproximateVarianceTopk has been found when elements are used up.\n");
//         //     for (int i = 0; i < topk; i++)
//         //     {
//         //         result_idx.push_back(feature_vec[i].idx_);
//         //         result_value.push_back(feature_vec[i].avg_);
//         //         printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
//         //     }
//         //     break;
//         // }
//         std::vector<Feature>::const_iterator iter = feature_vec.begin();
//         while (iter != feature_vec.end()) 
//         {
//             if (iter->up_ < kth_low) 
//             {
//                 iter = feature_vec.erase(iter);
//             }
//             else 
//             {
//                 ++iter;
//             }
//         }
//     } 
    
// }

void SelectAbsPredTopkVariance(const std::string mode, const int topk, const PredPair& pred, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, std::vector<int>& result_idx, const double pf, const double abs_error, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    uint64_t count = 0, num_total_scan = 0;
    SampleSizePair sample_size_pair(0, initial_sample_size);

    const double pf0 = pf * 1.0 / num_feature / ceil(log2(num_element * 1.0 / initial_sample_size)), require_diff = abs_error * 2;
    double access_data_time = 0.0, sort_time = 0.0, sample_index_time = 0.0;
    printf("Element %d Feature %d pf0 %1.15f\n", num_element, num_feature, pf0);
    std::vector<std::pair<double, double>> cumulate_pair(num_feature, std::pair<double, double>(0.0, 0.0));
    std::vector<Feature> feature_vec;
    std::vector<uint32_t> sample_vec;
    
    for (int i = 0; i < num_feature; i++)
    {
        if (i != pred.col_)
        {
            feature_vec.push_back(Feature(i, 0.0));
        }
    }
    result_idx.clear();
    result_value.clear();
    while (!feature_vec.empty() && sample_size_pair.current_ <= num_element)
    {
        uint32_t num_scan = 0;
        std::chrono::high_resolution_clock::time_point start, end; 
        std::chrono::duration<double> duration;
        start = std::chrono::high_resolution_clock::now();

        sample_vec.resize(sample_size_pair.current_);
        for (int i = 0; i < sample_vec.size(); i++)
        {
            sample_vec[i] = rand() % num_element;
        }

        ApproxSamplePredIndex(sample_size_pair.current_, sample_vec, pred.value_, data_mat[pred.col_], num_scan, sample_index_time);

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        sort_time += duration.count();

        for (int i = 0; i < feature_vec.size(); i++)
        {
            int feature_idx = feature_vec[i].idx_;
            VarianceEst variance(sample_vec, mode, sample_size_pair, cumulate_pair[feature_idx], data_mat[feature_idx], pf0, access_data_time);
            feature_vec[i].Update(variance.low_, variance.up_, variance.sample_variance_);
        }
        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_low);
        double kth_low = feature_vec[topk - 1].low_;
        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
        sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
        double max_diff = std::max_element(feature_vec.begin(), feature_vec.begin() + topk, cmp_diff)->diff_;
        printf("SampleSize %d ExistedFeature %d\n", sample_size_pair.current_, int(feature_vec.size()));
        for (int i = 0; i < topk; i++)
        {
            printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f requireDiff %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, require_diff);
        } 
        count += sample_size_pair.current_ * feature_vec.size();
        num_total_scan += num_scan;
        if (max_diff <= require_diff)
        {
            printf("ApproximateVarianceTopk has been found.\n");
            sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_low);
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f requireDiff %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, require_diff);
            }
            printf("SampleSize %d Count %ld TotalScan %lu SampleDataTime %f AccessDataTime %f\n", sample_size_pair.current_, count, num_total_scan, sample_index_time, access_data_time);
            break;
        }
        else if (sample_size_pair.current_ >= num_element / 2 && sample_size_pair.current_ < num_element)
        {
            sample_size_pair.previous_ = sample_size_pair.current_;
            sample_size_pair.current_ = num_element; 
            
            std::chrono::high_resolution_clock::time_point start, end; 
            std::chrono::duration<double> duration;
            start = std::chrono::high_resolution_clock::now();
            
            for (int i = 0; i < feature_vec.size(); i++)
            {
                int feature_idx = feature_vec[i].idx_;
                double exact_variance = ExactPredVarianceScore(pred.value_, data_mat[pred.col_], data_mat[feature_idx]);
                feature_vec[i].Update(exact_variance, exact_variance, exact_variance);
            }
            std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
            sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
            printf("ApproximateVarianceTopk has been found when elements are used up.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f requireDiff %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, require_diff);
            }
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            sort_time += duration.count();

            printf("SampleSize %d Count %ld TotalScan %lu SampleIndexTime %f AccessDataTime %f\n", sample_size_pair.current_, count, num_total_scan, sample_index_time, access_data_time);
            break;
        }
        else
        {
            sample_size_pair.previous_ = sample_size_pair.current_;
            sample_size_pair.current_ *= 2; 
        }
        
        // if (sample_size_pair.previous_ >= num_element)
        // {
        //     printf("ApproximateVarianceTopk has been found when elements are used up.\n");
        //     for (int i = 0; i < topk; i++)
        //     {
        //         result_idx.push_back(feature_vec[i].idx_);
        //         result_value.push_back(feature_vec[i].avg_);
        //         printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        //     }
        //     break;
        // }
        std::vector<Feature>::const_iterator iter = feature_vec.begin();
        while (iter != feature_vec.end()) 
        {
            if (iter->up_ < kth_low) 
            {
                iter = feature_vec.erase(iter);
            }
            else 
            {
                ++iter;
            }
        }
    } 
    
}

void EpsilonDeltaApproxVariance(const std::string& mode, const double& empirical_threshold, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size(), num_element = data_mat[0].size(), num_test_feature = test_feature.size();
    long count = 0; 
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    printf("Element %d NumTestFeature %d pf0 %1.15f\n", num_element, num_test_feature, pf0);
     
    for (int j = 0; j < test_feature.size(); j++)
    {
        Feature single_feature(test_feature[j], 0.0);
        SampleSizePair sample_size_pair(0, initial_sample_size);
        std::pair<double, double> cum_pair(0.0, 0.0);
        while (sample_size_pair.current_ < num_element)
        {
            std::vector<int> sample_vec;
            UpdateCumRandom(sample_size_pair.current_, cum_pair, data_mat[test_feature[j]], precise_time);
            VarianceEst variance(mode, sample_size_pair, cum_pair, data_mat[test_feature[j]], pf0, precise_time);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
        
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", sample_size_pair.current_, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
            count += sample_size_pair.current_;
            if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
            {
                printf("ApproximateVariance has been found.\n");
                result_value[j] = single_feature.est_; 
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, precise_time, median_time, sort_time);
                break;
            }
            else if (single_feature.diff_ <= 2.0 * error_bound * empirical_threshold && single_feature.low_ < empirical_threshold)
            {
                printf("ApproximateValue has been found when the lower is smaller than the empirical threshold.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonDelta %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * empirical_threshold);
                printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, precise_time, median_time, sort_time);
                break;
            }
            else if (sample_size_pair.current_ >= num_element / 2 && sample_size_pair.current_ < num_element)
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ = num_element; 
                
                start = std::chrono::high_resolution_clock::now();
                
                double exact_variance = ExactVarianceScore(data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);

                printf("ApproximateVariance has been found when elements are used up.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], sample_size_pair.current_, single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, precise_time, median_time, sort_time);
                break;
            }
            else
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ *= 2; 
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", i, result_value[i]);   
    }
}

void SelectAbsPredApproxVariance(const std::string& mode, const PredPair& pred, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const double pf, const double abs_error, std::vector<double>& result_value)
{
    const uint8_t num_test_feature = test_feature.size();
    const uint32_t num_element = data_mat[0].size(); 
    uint64_t count = 0, num_total_scan = 0;
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
    // select_pr = 0.0, pract_error = 0.0;
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    printf("Element %d NumTestFeature %d pf0 %1.15f\n", num_element, num_test_feature, pf0);
    
     
    for (uint8_t j = 0; j < test_feature.size(); j++)
    {
        Feature single_feature(test_feature[j], 0.0);
        SampleSizePair sample_size_pair(0, initial_sample_size);
        std::pair<double, double> cum_pair(0.0, 0.0);
        while (sample_size_pair.current_ < num_element)
        {
            uint32_t num_scan = 0;
            UpdateCumRandPred(sample_size_pair.current_, cum_pair, pred.value_, data_mat[pred.col_], data_mat[test_feature[j]], num_scan, precise_time);
            // select_pr = 1.0 * sample_size_pair.current_ / num_scan;

            VarianceEst variance(mode, sample_size_pair, cum_pair, data_mat[test_feature[j]], pf0, precise_time);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
            // pract_error = abs_error / sqrt(select_pr);
        
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2AbsError %1.10f\n", sample_size_pair.current_, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * abs_error);
            count += sample_size_pair.current_;
            num_total_scan += num_scan;
            if (single_feature.diff_ <= 2.0 * abs_error)
            {
                printf("PredApproxiVariance has been found.\n");
                result_value[j] = single_feature.est_; 
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2AbsError %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * abs_error);
                printf("SampleSize %d Count %ld TotalScan %lu PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, num_total_scan, precise_time, median_time, sort_time);
                break;
            }
            else if (num_scan > num_element / 2)
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ = num_element; 
                
                start = std::chrono::high_resolution_clock::now();
                
                double exact_variance = ExactPredVarianceScore(pred.value_, data_mat[pred.col_], data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);

                printf("PredApproxVariance has been found when elements are used up.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2AbsError %1.10f\n", test_feature[j], sample_size_pair.current_, single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * abs_error);
                
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %lu TotalScan %lu PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, num_total_scan, precise_time, median_time, sort_time);
                break;
            }
            else
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ *= 2; 
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", test_feature[i], result_value[i]);   
    }
}

void SelectAbsPredRandomVarianceStat(const PredPair& pred, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const double pf, const double abs_error, std::vector<double>& result_value)
{
    const uint8_t num_test_feature = test_feature.size();
    const uint32_t num_element = data_mat[0].size(); 
    uint64_t count = 0, num_total_scan = 0;
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
    // select_pr = 0.0, pract_error = 0.0;
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    printf("Element %d NumTestFeature %d pf0 %1.15f\n", num_element, num_test_feature, pf0);
    
     
    for (uint8_t j = 0; j < test_feature.size(); j++)
    {
        Feature single_feature(test_feature[j], 0.0);
        SampleSizePair sample_size_pair(0, initial_sample_size);
        std::pair<double, double> cum_pair(0.0, 0.0);
        while (sample_size_pair.current_ < num_element)
        {
            uint32_t num_scan = 0, sample_size_count = 0;
            cum_pair = std::pair<double, double>(0.0, 0.0);
    // int sample = rand() % data_vec.size();
    // printf("size %lu sample %d\n", data_vec.size(), sample);
            while(sample_size_count < sample_size_pair.current_)
            {
                int sample_index = rand() % data_mat[test_feature[j]].size();
                num_scan ++;
                if (abs(data_mat[pred.col_][sample_index] - pred.value_) <= 1e-12)
                {
                    cum_pair.first += data_mat[test_feature[j]][sample_index];
                    cum_pair.second += pow(data_mat[test_feature[j]][sample_index], 2.0);
                    sample_size_count ++;
                }
            }
            // select_pr = 1.0 * sample_size_pair.current_ / num_scan;

            VarianceEst variance("random", sample_size_pair, cum_pair, data_mat[test_feature[j]], pf0, precise_time);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
            // pract_error = abs_error / sqrt(select_pr);
        
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2AbsError %1.10f\n", sample_size_pair.current_, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * abs_error);
            count += sample_size_pair.current_;
            num_total_scan += num_scan;
            if (single_feature.diff_ <= 2.0 * abs_error)
            {
                printf("PredApproxiVariance has been found.\n");
                result_value[j] = single_feature.est_; 
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2AbsError %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * abs_error);
                printf("SampleSize %d Count %ld TotalScan %lu PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, num_total_scan, precise_time, median_time, sort_time);
                break;
            }
            else if (num_scan > num_element / 2)
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ = num_element; 
                
                start = std::chrono::high_resolution_clock::now();
                
                double exact_variance = ExactPredVarianceScore(pred.value_, data_mat[pred.col_], data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);

                printf("PredApproxVariance has been found when elements are used up.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2AbsError %1.10f\n", test_feature[j], sample_size_pair.current_, single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * abs_error);
                
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %lu TotalScan %lu PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, num_total_scan, precise_time, median_time, sort_time);
                break;
            }
            else
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ *= 2; 
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", test_feature[i], result_value[i]);   
    }
}

void SelectAbsPredColtVarianceStat(const PredPair& pred, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const double pf, const double abs_error, std::vector<double>& result_value)
{
    const uint8_t num_test_feature = test_feature.size();
    const uint32_t num_element = data_mat[0].size(); 
    uint64_t count = 0, num_total_scan = 0;
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
    // select_pr = 0.0, pract_error = 0.0;
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    printf("Element %d NumTestFeature %d pf0 %1.15f\n", num_element, num_test_feature, pf0);
    
     
    for (uint8_t j = 0; j < test_feature.size(); j++)
    {
        Feature single_feature(test_feature[j], 0.0);
        SampleSizePair sample_size_pair(0, initial_sample_size);
        std::pair<double, double> cum_pair(0.0, 0.0);
        while (sample_size_pair.current_ < num_element)
        {
            uint32_t num_scan = 0, sample_size_count = 0;
            // select_pr = 1.0 * sample_size_pair.current_ / num_scan;
            // uint32_t ;
            cum_pair = std::pair<double, double>(0.0, 0.0);
            // int sample = rand() % data_vec.size();
            // printf("size %lu sample %d\n", data_vec.size(), sample);
            while(sample_size_count < sample_size_pair.current_)
            {
                int sample_index = rand() % data_mat[test_feature[j]].size();
                num_scan ++;
                if (abs(data_mat[pred.col_][sample_index] - pred.value_) <= 1e-12)
                {
                    cum_pair.first += data_mat[test_feature[j]][sample_index];
                    cum_pair.second += pow(data_mat[test_feature[j]][sample_index], 2.0);
                    sample_size_count ++;
                }
            }


            VarianceEst variance("colt", sample_size_pair, cum_pair, data_mat[test_feature[j]], pf0, precise_time);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
            // pract_error = abs_error / sqrt(select_pr);
        
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2AbsError %1.10f\n", sample_size_pair.current_, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * abs_error);
            count += sample_size_pair.current_;
            num_total_scan += num_scan;
            if (single_feature.diff_ <= 2.0 * abs_error)
            {
                printf("PredApproxiVariance has been found.\n");
                result_value[j] = single_feature.est_; 
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2AbsError %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * abs_error);
                printf("SampleSize %d Count %ld TotalScan %lu PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, num_total_scan, precise_time, median_time, sort_time);
                break;
            }
            else if (num_scan > num_element / 2)
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ = num_element; 
                
                start = std::chrono::high_resolution_clock::now();
                
                double exact_variance = ExactPredVarianceScore(pred.value_, data_mat[pred.col_], data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);

                printf("PredApproxVariance has been found when elements are used up.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2AbsError %1.10f\n", test_feature[j], sample_size_pair.current_, single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * abs_error);
                
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %lu TotalScan %lu PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, num_total_scan, precise_time, median_time, sort_time);
                break;
            }
            else
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ *= 2; 
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", test_feature[i], result_value[i]);   
    }
}

void EpsilonDeltaPredApproxVariance(const std::string& mode, const PredPair& pred, const double& empirical_threshold, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size(), num_element = data_mat[0].size(), num_test_feature = test_feature.size();
    long count = 0; 
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    printf("Element %d NumTestFeature %d pf0 %1.15f\n", num_element, num_test_feature, pf0);
    uint num_total_scan = 0;
     
    for (int j = 0; j < test_feature.size(); j++)
    {
        uint num_scan = 0;
        Feature single_feature(test_feature[j], 0.0);
        SampleSizePair sample_size_pair(0, initial_sample_size);
        std::pair<double, double> cum_pair(0.0, 0.0);
        while (sample_size_pair.current_ < num_element)
        {
            // std::vector<uint32_t> sample_vec;
            UpdateCumRandPred(sample_size_pair.current_, cum_pair, pred.value_, data_mat[pred.col_], data_mat[test_feature[j]], num_total_scan, precise_time);
            VarianceEst variance(mode, sample_size_pair, cum_pair, data_mat[test_feature[j]], pf0, precise_time);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
        
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", sample_size_pair.current_, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
            count += sample_size_pair.current_;
            num_total_scan += num_scan;
            if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
            {
                printf("PredApproxiVariance has been found.\n");
                result_value[j] = single_feature.est_; 
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                printf("SampleSize %d Count %ld TotalScan %u PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, num_total_scan, precise_time, median_time, sort_time);
                break;
            }
            else if (single_feature.diff_ <= 2.0 * error_bound * empirical_threshold && single_feature.low_ < empirical_threshold)
            {
                printf("PredApproxValue has been found when the lower is smaller than the empirical threshold.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonDelta %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * empirical_threshold);
                printf("SampleSize %d Count %ld TotalScan %u PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, num_total_scan, precise_time, median_time, sort_time);
                break;
            }
            else if (num_scan >= num_element / 2)
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ = num_element; 
                
                start = std::chrono::high_resolution_clock::now();
                
                double exact_variance = ExactPredVarianceScore(pred.value_, data_mat[pred.col_], data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);

                printf("PredApproxVariance has been found when elements are used up.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], sample_size_pair.current_, single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %ld TotalScan %u PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, num_total_scan, precise_time, median_time, sort_time);
                break;
            }
            else
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ *= 2; 
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", test_feature[i], result_value[i]);   
    }
}

void EpsilonDeltaPredApproxVarianceBlock2(const std::string& mode, const PredPair& pred, const double& empirical_threshold, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size(), num_element = data_mat[0].size(), num_test_feature = test_feature.size();
    long count = 0; 
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), sample_index_time = 0.0, access_data_time = 0.0, sort_time = 0.0;
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    printf("Element %d NumTestFeature %d pf0 %1.15f\n", num_element, num_test_feature, pf0);
     
    for (int j = 0; j < test_feature.size(); j++)
    {
        Feature single_feature(test_feature[j], 0.0);
        SampleSizePair sample_size_pair(0, initial_sample_size);
        std::pair<double, double> cum_pair(0.0, 0.0);
        while (sample_size_pair.current_ < num_element)
        {
            std::vector<int> sample_vec;
            UpdateCumRandPredBlock2(sample_size_pair.current_, sample_vec, cum_pair, pred.value_, data_mat[pred.col_], data_mat[test_feature[j]], sample_index_time, access_data_time);
            VarianceEst variance(mode, sample_size_pair, cum_pair, data_mat[test_feature[j]], pf0, access_data_time);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
        
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", sample_size_pair.current_, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
            count += sample_size_pair.current_;
            if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
            {
                printf("PredApproxiVariance has been found.\n");
                result_value[j] = single_feature.est_; 
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f\n", sample_size_pair.current_, count, sample_index_time, access_data_time);
                break;
            }
            else if (single_feature.diff_ <= 2.0 * error_bound * empirical_threshold && single_feature.low_ < empirical_threshold)
            {
                printf("PredApproxValue has been found when the lower is smaller than the empirical threshold.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonDelta %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * empirical_threshold);
                printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f\n", sample_size_pair.current_, count, sample_index_time, access_data_time);
                break;
            }
            else if (sample_size_pair.current_ >= num_element / 2 && sample_size_pair.current_ < num_element)
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ = num_element; 
                
                start = std::chrono::high_resolution_clock::now();
                
                double exact_variance = ExactVarianceScore(data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);

                printf("PredApproxVariance has been found when elements are used up.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], sample_size_pair.current_, single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f SortTime %f\n", sample_size_pair.current_, count, sample_index_time, access_data_time, sort_time);
                break;
            }
            else
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ *= 2; 
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", test_feature[i], result_value[i]);   
    }
}

// void ApproxMetric(const std::string& metric, const std::string& mode, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const double pf, const double error_bound, std::vector<double>& result_value)
// {
//     const int num_feature = data_mat.size(), num_element = data_mat[0].size(), num_test_feature = test_feature.size();
//     long count = 0; 
//     if (result_value.size() != test_feature.size())
//     {
//         result_value.resize(num_test_feature);
//     }
//     double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
//     std::chrono::high_resolution_clock::time_point start, end; 
//     std::chrono::duration<double> duration;
//     printf("Element %d NumTestFeature %d pf0 %1.15f\n", num_element, num_test_feature, pf0);
     
//     for (int j = 0; j < test_feature.size(); j++)
//     {
//         Feature single_feature(test_feature[j], 0.0);
//         SampleSizePair sample_size_pair(0, initial_sample_size);
//         std::pair<double, double> cumulate_pair(0.0, 0.0);
//         while (sample_size_pair.current_ < num_element)
//         {
//             start = std::chrono::high_resolution_clock::now();
//             std::vector<int> sample_vec(sample_size_pair.current_);
//             for (int i = 0; i < sample_vec.size(); i++)
//             {
//                 sample_vec[i] = rand() % num_element;
//             }
//             end = std::chrono::high_resolution_clock::now();
//             duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//             median_time += duration.count();

//             if (metric == "variance")
//             {
//                 VarianceEst variance(sample_vec, mode, sample_size_pair, cumulate_pair, data_mat[test_feature[j]], pf0, precise_time);
//                 single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
//             }
//             // else if (metric == "average")
//             // {
//             //     AverageEst average(sample_vec, sample_size_pair, cumulate_pair.first, data_mat[test_feature[j]], pf0, precise_time);
//             //     single_feature.Update(average.low_, average.up_);
//             // }
//             else
//             {
//                 printf("Error Metric!\n");
//             }
            
//             printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", sample_size_pair.current_, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
//             count += sample_size_pair.current_;
//             if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
//             {
//                 printf("ApproximateVariance has been found.\n");
//                 result_value[j] = single_feature.est_; 
//                 printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
//                 printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, precise_time, median_time, sort_time);
//                 break;
//             }
//             else if (sample_size_pair.current_ >= num_element / 2 && sample_size_pair.current_ < num_element)
//             {
//                 sample_size_pair.previous_ = sample_size_pair.current_;
//                 sample_size_pair.current_ = num_element; 
                
//                 start = std::chrono::high_resolution_clock::now();
                
//                 if (metric == "variance")
//                 {
//                     double exact_variance = ExactVarianceScore(data_mat[test_feature[j]]);
//                     single_feature.Update(exact_variance, exact_variance, exact_variance);
//                 }
//                 // else if (metric == "average")
//                 // {
//                 //     double exact_average = ExactAverageScore(data_mat[test_feature[j]]);
//                 //     single_feature.Update(exact_average, exact_average);
//                 // }
//                 printf("ApproximateVariance has been found when elements are used up.\n");
//                 result_value[j] = single_feature.est_;
//                 printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], sample_size_pair.current_, single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                
//                 end = std::chrono::high_resolution_clock::now();
//                 duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//                 sort_time += duration.count();

//                 printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, precise_time, median_time, sort_time);
//                 break;
//             }
//             else
//             {
//                 sample_size_pair.previous_ = sample_size_pair.current_;
//                 sample_size_pair.current_ *= 2; 
//             }
//         } 
//     }
//     printf("\n");
//     for (int i = 0; i < result_value.size(); i++)
//     {
//         printf("idx %d %1.12f\n", i, result_value[i]);   
//     }
// }

// void BernsteinMetric(const std::string& metric, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const double pf, const double error_bound, std::vector<double>& result_value, const int block_size)
// {
//     const int num_feature = data_mat.size(), num_element = data_mat[0].size(), num_test_feature = test_feature.size(); 
//     // total_num_block = num_element / block_size;
//     long count = 0; 
//     if (result_value.size() != test_feature.size())
//     {
//         result_value.resize(num_test_feature);
//     }
//     double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
//     std::chrono::high_resolution_clock::time_point start, end; 
//     std::chrono::duration<double> duration;
//     printf("Element %d NumTestFeature %d pf0 %1.15f BernsteinBlockSize %d\n", num_element, num_test_feature, pf0, block_size);
     
//     for (int j = 0; j < test_feature.size(); j++)
//     {
//         Feature single_feature(test_feature[j], 0.0);
//         SampleSizePair sample_size_pair(0, initial_sample_size);
//         double cumulate_value;
//         while (sample_size_pair.current_ < num_element)
//         {
//             start = std::chrono::high_resolution_clock::now();
//             std::vector<int> sample_vec(sample_size_pair.current_ / block_size);
//             for (int i = 0; i < sample_vec.size(); i++)
//             {
//                 sample_vec[i] = rand() % num_element;
//             }
//             end = std::chrono::high_resolution_clock::now();
//             duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//             median_time += duration.count();

//             if (metric == "average")
//             {
//                 AverageEst average(sample_vec, block_size, sample_size_pair, cumulate_value, data_mat[test_feature[j]], pf0, precise_time);
//                 single_feature.Update(average.low_, average.up_);
//             }
//             else
//             {
//                 printf("Error Metric!\n");
//             }
            
//             printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", sample_size_pair.current_, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
//             count += sample_size_pair.current_;
//             if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
//             {
//                 printf("ApproximateVariance has been found.\n");
//                 result_value[j] = single_feature.avg_; 
//                 printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
//                 printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, precise_time, median_time, sort_time);
//                 break;
//             }
//             else if (sample_size_pair.current_ >= num_element / 2 && sample_size_pair.current_ < num_element)
//             {
//                 sample_size_pair.previous_ = sample_size_pair.current_;
//                 sample_size_pair.current_ = num_element; 
                
//                 start = std::chrono::high_resolution_clock::now();
                
//                 if (metric == "variance")
//                 {
//                     double exact_variance = ExactVarianceScore(data_mat[test_feature[j]]);
//                     single_feature.Update(exact_variance, exact_variance);
//                 }
//                 else if (metric == "average")
//                 {
//                     double exact_average = ExactAverageScore(data_mat[test_feature[j]]);
//                     single_feature.Update(exact_average, exact_average);
//                 }
//                 printf("ApproximateVariance has been found when elements are used up.\n");
//                 result_value[j] = single_feature.avg_;
//                 printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], sample_size_pair.current_, single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                
//                 end = std::chrono::high_resolution_clock::now();
//                 duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//                 sort_time += duration.count();

//                 printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", sample_size_pair.current_, count, precise_time, median_time, sort_time);
//                 break;
//             }
//             else
//             {
//                 sample_size_pair.previous_ = sample_size_pair.current_;
//                 sample_size_pair.current_ *= 2; 
//             }
//         } 
//     }
// }

void MomVarianceTopkPred(const double ratio_num_block, const int topk, const PredPair& pred, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, std::vector<int>& result_idx, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    long count = 0; 
    int block_size = initial_sample_size, num_slot = cum_size.size() - 1;
    double pf0 = pf * 1.0 / num_feature / ceil(log2(num_element * 1.0 / initial_sample_size)), sample_index_time = 0.0, access_data_time = 0.0, sort_time = 0.0;
    const int num_var = ceil(ratio_num_block * log(2.0 / pf0));
    std::vector<uint32_t> sample_vec;

    printf("Element %d Feature %d pf0 %1.15f 1_pf0 %f NumVar %d BlockSize %d Block3.5\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var, block_size);

    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    
    std::vector<Feature> feature_vec;
    std::vector<uint32_t> sample_start_pos(num_var, 0);
    std::vector<std::pair<std::vector<double>, std::vector<double>>> MoM_pair(num_feature, std::pair<std::vector<double>, std::vector<double>>(std::vector<double>(num_var), std::vector<double>(num_var)));
    for (uint8_t i = 0; i < num_feature; i++)
    {
        if (i != pred.col_)
        {
            feature_vec.push_back(Feature(i, 0.0));
        }    
    }
    result_idx.clear();
    result_value.clear();
    while (!feature_vec.empty() && block_size * num_var <= num_element)
    {
        int existed_feature = 0;
        const int num_block = num_element / block_size;
        for (uint8_t i = 0; i < sample_start_pos.size(); i++)
        {
            sample_start_pos[i] = rand() % num_element;
        }
        std::sort(sample_start_pos.begin(), sample_start_pos.end());
        // printf("finish sort for block size %d\n", block_size);
        MomSamplePredIndex(block_size, sample_start_pos, sample_vec, pred.value_, data_mat[pred.col_], data_mat[num_feature], cum_size, sample_index_time);

        for (int i = 0; i < feature_vec.size(); i++)
        {
            if (feature_vec[i].ratio_ < (1 - error_bound))
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                // printf("Feature %d\n", feature_idx);
                VarianceEst variance(num_var, block_size, sample_vec, MoM_pair[feature_idx], data_mat[feature_idx], pf0, access_data_time);
                // printf("special idx %3d low %2.6f up %1.6f\n", feature_vec[i].idx_, variance.low_, variance.up_);
                feature_vec[i].Update(variance.low_, variance.up_, variance.sample_variance_);
            }
            // if (feature_vec[i].idx_ == 107 || feature_vec[i].idx_ == 1)
            // {
            //     printf("special idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            // }
        }

        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_low);
        double kth_low = feature_vec[topk - 1].low_;
        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
        sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
        double min_ratio = std::min_element(feature_vec.begin(), feature_vec.begin() + topk, cmp_ratio)->ratio_;
        printf("SampleSize %d ExistedFeature %d\n", block_size * num_var, existed_feature);

        for (int i = 0; i < topk; i++)
        {
            printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        } 
        count += block_size * num_var * existed_feature;
        if (min_ratio >= (1 - error_bound))
        {
            printf("ApproximateVarianceTopk has been found.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }
            printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f\n", block_size * num_var, count, sample_index_time, access_data_time);
            break;
        }
        else if (block_size * num_var > num_element / 2)
        {
            std::chrono::high_resolution_clock::time_point start, end; 
            std::chrono::duration<double> duration;
            start = std::chrono::high_resolution_clock::now();

            for (int i = 0; i < feature_vec.size(); i++)
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                double exact_variance = ExactVarianceScore(data_mat[feature_idx]);
                feature_vec[i].Update(exact_variance, exact_variance, exact_variance);
            }
            std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
            sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
            printf("ApproximateVarianceTopk has been found when elements are used up.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }

            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            sort_time += duration.count();

            printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f SortTime %f\n", num_element, count, sample_index_time, access_data_time, sort_time);
            break;
        }
        else
        {
            block_size *= 2;
        }
        
        std::vector<Feature>::const_iterator iter = feature_vec.begin();
        while (iter != feature_vec.end()) 
        {
            if (iter->up_ < kth_low) 
            {
                iter = feature_vec.erase(iter);
            }
            else 
            {
                ++iter;
            }
        }
    } 
    // printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var / 2, count, precise_time, median_time, sort_time);
}

void SelectAbsPredMomTopkVariance(const double ratio_num_block, const int topk, const PredPair& pred, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, std::vector<int>& result_idx, const double pf, const double abs_error, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    uint64_t count = 0, num_total_scan = 0; 
    int block_size = initial_sample_size, num_slot = cum_size.size() - 1;
    const double pf0 = pf * 1.0 / num_feature / ceil(log2(num_element * 1.0 / initial_sample_size)), require_diff = abs_error * 2;
    double sample_index_time = 0.0, access_data_time = 0.0, sort_time = 0.0;
    const int num_var = ceil(ratio_num_block * log(2.0 / pf0));
    std::vector<uint32_t> sample_vec;

    printf("Element %d Feature %d pf0 %1.15f 1_pf0 %f NumVar %d BlockSize %d Block3.5\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var, block_size);

    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    
    std::vector<Feature> feature_vec;
    std::vector<uint32_t> sample_start_pos(num_var, 0);
    std::vector<std::pair<std::vector<double>, std::vector<double>>> MoM_pair(num_feature, std::pair<std::vector<double>, std::vector<double>>(std::vector<double>(num_var), std::vector<double>(num_var)));
    for (uint8_t i = 0; i < num_feature; i++)
    {
        if (i != pred.col_)
        {
            feature_vec.push_back(Feature(i, 0.0));
        }    
    }
    result_idx.clear();
    result_value.clear();
    while (!feature_vec.empty() && block_size * num_var <= num_element)
    {
        int existed_feature = 0;
        uint32_t num_scan = 0;
        const int num_block = num_element / block_size;
        for (uint8_t i = 0; i < sample_start_pos.size(); i++)
        {
            sample_start_pos[i] = rand() % num_element;
        }
        std::sort(sample_start_pos.begin(), sample_start_pos.end());
        // printf("finish sort for block size %d\n", block_size);
        MomSamplePredIndexBlock2(block_size, sample_start_pos, sample_vec, pred.value_, data_mat[pred.col_], data_mat[num_feature], num_scan, cum_size, sample_index_time);
        printf("Finish sample index!\n");
        for (int i = 0; i < feature_vec.size(); i++)
        {
            if (feature_vec[i].diff_ > require_diff)
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                VarianceEst variance(num_var, block_size, sample_vec, MoM_pair[feature_idx], data_mat[feature_idx], pf0, access_data_time);
                feature_vec[i].Update(variance.low_, variance.up_, variance.sample_variance_);
            }
        }
        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_low);
        double kth_low = feature_vec[topk - 1].low_;
        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
        sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
        double max_diff = std::max_element(feature_vec.begin(), feature_vec.begin() + topk, cmp_diff)->diff_;
        printf("SampleSize %d ExistedFeature %d NumOfScan %d\n", block_size * num_var, existed_feature, num_scan);

        for (int i = 0; i < topk; i++)
        {
            printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f requireDiff %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, require_diff);
        } 
        count += block_size * num_var * existed_feature;
        num_total_scan += num_scan;
        if (max_diff <= require_diff)
        {
            printf("ApproximateVarianceTopk has been found.\n");
            sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_low);
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f requireDiff %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, require_diff);
            }
            printf("SampleSize %d Count %ld TotalScan %lu SampleIndexTime %f AccessDataTime %f\n", block_size * num_var, count, num_total_scan, sample_index_time, access_data_time);
            break;
        }
        else if (block_size * num_var > num_element / 2)
        {
            std::chrono::high_resolution_clock::time_point start, end; 
            std::chrono::duration<double> duration;
            start = std::chrono::high_resolution_clock::now();

            for (int i = 0; i < feature_vec.size(); i++)
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                double exact_variance = ExactPredVarianceScore(pred.value_, data_mat[pred.col_], data_mat[feature_idx]);
                feature_vec[i].Update(exact_variance, exact_variance, exact_variance);
            }
            std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
            sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
            printf("ApproximateVarianceTopk has been found when elements are used up.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f requireDiff %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, require_diff);
            }

            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            sort_time += duration.count();

            printf("SampleSize %d Count %ld TotalScan %lu SampleIndexTime %f AccessDataTime %f SortTime %f\n", num_element, count, num_total_scan, sample_index_time, access_data_time, sort_time);
            break;
        }
        else
        {
            block_size *= 2;
        }
        
        std::vector<Feature>::const_iterator iter = feature_vec.begin();
        while (iter != feature_vec.end()) 
        {
            if (iter->up_ < kth_low) 
            {
                iter = feature_vec.erase(iter);
            }
            else 
            {
                ++iter;
            }
        }
    } 
    // printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var / 2, count, precise_time, median_time, sort_time);
}

void VarianceTopkMoMSimple(const double ratio_num_block, const int topk, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, std::vector<int>& sample_index_count, std::vector<int>& result_idx, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    long count = 0; 
    int block_size = initial_sample_size, num_slot = cum_size.size() - 1;
    double pf0 = pf * 1.0 / num_feature / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
    const int num_var = ceil(ratio_num_block * log(2.0 / pf0));

    printf("Element %d Feature %d pf0 %1.15f 1_pf0 %f NumVar %d BlockSize %d\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var, block_size);

    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    
    std::vector<Feature> feature_vec;
    std::vector<uint32_t> sample_start_pos(num_var, 0);
    std::vector<std::pair<std::vector<double>, std::vector<double>>> MoM_pair(num_feature, std::pair<std::vector<double>, std::vector<double>>(std::vector<double>(num_var), std::vector<double>(num_var)));
    for (uint8_t i = 0; i < num_feature; i++)
    {
        feature_vec.push_back(Feature(i, 0.0));
    }
    result_idx.clear();
    result_value.clear();
    while (!feature_vec.empty() && block_size * num_var <= num_element)
    {
        //sample_vec.clear();
        // bool last_batch = false;
        int existed_feature = 0;
        const int num_block = num_element / block_size;

        for (uint8_t i = 0; i < sample_start_pos.size(); i++)
        {
            sample_start_pos[i] = rand() % num_element;
        }
        std::sort(sample_start_pos.begin(), sample_start_pos.end());
        printf("finish sort for block size %d\n", block_size);
        for (int i = 0; i < sample_start_pos.size(); i++)
        {
            bool exceed_final_slot = false;
            int start_slot = data_mat[num_feature][sample_start_pos[i]], end_slot, sample_stop_pos;
            // printf("initialize start slot\n");
            
            // determine the final sample position
            end_slot = data_mat[num_feature][sample_start_pos[i] + block_size - 1];
            sample_stop_pos = (sample_start_pos[i] + block_size - 1) % num_element;
            if (start_slot > end_slot)
            {
                exceed_final_slot = true;
            }
        
            // printf("initialize end slot\n");

            // first_slot
            int start_slot_first = cum_size[start_slot];
            std::fill(sample_index_count.begin() + cum_size[start_slot], sample_index_count.begin() + cum_size[start_slot + 1], 0);
            // printf("initialize start_slot\n");
            for (int j = 0; j < cum_size[start_slot + 1] - sample_start_pos[i]; j++)
            {
                int random_index = cum_size[start_slot] + Next() % (cum_size[start_slot + 1] - cum_size[start_slot]);
                sample_index_count[random_index] += 1;
            }
            // printf("sample start_slot\n");
            if (exceed_final_slot == false)
            {
                std::fill(sample_index_count.begin() + cum_size[start_slot + 1], sample_index_count.begin() + cum_size[end_slot], 0);
                for (int j = start_slot + 1; j < end_slot; j++)
                {
                    uint8_t slot_size = cum_size[j + 1] - cum_size[j];

                    if (slot_size == 1)
                    {
                        sample_index_count[cum_size[j]] += 1;
                    }
                    else
                    {
                        for (uint8_t k = 0; k < slot_size; k++)
                        {
                            sample_index_count[cum_size[j] + Next() % slot_size] += 1;
                        }
                    }
                }
            }
            else
            {
                std::fill(sample_index_count.begin() + cum_size[start_slot + 1], sample_index_count.end(), 0);
                for (int j = start_slot + 1; j < cum_size.size(); j++)
                {
                    uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                    if (slot_size == 1)
                    {
                        sample_index_count[cum_size[j]] += 1;
                    }
                    else
                    {
                        for (uint8_t k = 0; k < slot_size; k++)
                        {
                            sample_index_count[cum_size[j] + Next() % slot_size] += 1;
                        }
                    }
                }
                std::fill(sample_index_count.begin(), sample_index_count.begin() + cum_size[end_slot], 0);
                for (int j = 0; j < end_slot; j++)
                {
                    uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                    if (slot_size == 1)
                    {
                        sample_index_count[cum_size[j]] += 1;
                    }
                    else
                    {
                        for (uint8_t k = 0; k < slot_size; k++)
                        {
                            sample_index_count[cum_size[j] + Next() % slot_size] += 1;
                        }
                    }
                }
            }

            // last slot
            uint32_t end_slot_last = cum_size[end_slot + 1];
            std::fill(sample_index_count.begin() + cum_size[end_slot], sample_index_count.begin() + cum_size[end_slot + 1], 0);
            uint8_t end_slot_sample_size = sample_stop_pos - cum_size[end_slot] + 1;
            // printf("initialize end_slot\n");
            for (int j = 0; j < end_slot_sample_size; j++)
            {
                sample_index_count[cum_size[end_slot] + Next() % (cum_size[end_slot + 1] - cum_size[end_slot])] += 1;
            }

            // printf("sample end_slot\n");
            for (int j = 0; j < feature_vec.size(); j++)
            {
                if (feature_vec[j].ratio_ < (1 - error_bound))
                {
                    int feature_idx = feature_vec[j].idx_;
                    UpdateSingleBlock(block_size, MoM_pair[feature_idx], i, data_mat[feature_idx], sample_index_count, start_slot_first, end_slot_last, exceed_final_slot, precise_time);
                }
            }
            // printf("sample records\n");
        }

        for (int i = 0; i < feature_vec.size(); i++)
        {
            if (feature_vec[i].ratio_ < (1 - error_bound))
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                // printf("Feature %d\n", feature_idx);
                VarianceEst variance(block_size, MoM_pair[feature_idx], pf0);
                // printf("special idx %3d low %2.6f up %1.6f\n", feature_vec[i].idx_, variance.low_, variance.up_);
                feature_vec[i].Update(variance.low_, variance.up_, variance.sample_variance_);
            }
            // if (feature_vec[i].idx_ == 107 || feature_vec[i].idx_ == 1)
            // {
            //     printf("special idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            // }
        }

        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_low);
        double kth_low = feature_vec[topk - 1].low_;
        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
        sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
        double min_ratio = std::min_element(feature_vec.begin(), feature_vec.begin() + topk, cmp_ratio)->ratio_;
        printf("SampleSize %d ExistedFeature %d\n", block_size * num_var, existed_feature);

        // if (block_size * num_var >= 2000000)
        // {
        //     for (int i = 0; i < feature_vec.size(); ++i)
        //     {
        //         printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        //     }
        // // }
        // printf("\n");

        for (int i = 0; i < topk; i++)
        {
            printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        } 
        count += block_size * num_var * existed_feature;
        if (min_ratio >= (1 - error_bound))
        {
            printf("ApproximateVarianceTopk has been found.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }
            printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var, count, precise_time, median_time, sort_time);
            break;
        }
        else if (block_size * num_var > num_element / 2)
        {
            std::chrono::high_resolution_clock::time_point start, end; 
            std::chrono::duration<double> duration;
            start = std::chrono::high_resolution_clock::now();

            for (int i = 0; i < feature_vec.size(); i++)
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                double exact_variance = ExactVarianceScore(data_mat[feature_idx]);
                feature_vec[i].Update(exact_variance, exact_variance, exact_variance);
            }
            std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
            sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
            printf("ApproximateVarianceTopk has been found when elements are used up.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }

            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            sort_time += duration.count();

            printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", num_element, count, precise_time, median_time, sort_time);
            break;
        }
        else
        {
            block_size *= 2;
        }
        
        std::vector<Feature>::const_iterator iter = feature_vec.begin();
        while (iter != feature_vec.end()) 
        {
            if (iter->up_ < kth_low) 
            {
                iter = feature_vec.erase(iter);
            }
            else 
            {
                ++iter;
            }
        }
    } 
    // printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var / 2, count, precise_time, median_time, sort_time);
}

void VarianceTopkMomBlock2(const double ratio_num_block, const int topk, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, std::vector<int>& result_idx, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    long count = 0; 
    int block_size = initial_sample_size, num_slot = cum_size.size() - 1;
    double pf0 = pf * 1.0 / num_feature / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
    const int num_var = ceil(ratio_num_block * log(2.0 / pf0));

    printf("Element %d Feature %d pf0 %1.15f 1_pf0 %f NumVar %d BlockSize %d NonSimpleMom\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var, block_size);

    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    
    std::vector<Feature> feature_vec;
    std::vector<uint32_t> sample_start_pos(num_var, 0);
    std::vector<std::pair<std::vector<double>, std::vector<double>>> MoM_pair(num_feature, std::pair<std::vector<double>, std::vector<double>>(std::vector<double>(num_var), std::vector<double>(num_var)));
    for (uint8_t i = 0; i < num_feature; i++)
    {
        feature_vec.push_back(Feature(i, 0.0));
    }
    result_idx.clear();
    result_value.clear();
    while (!feature_vec.empty() && block_size * num_var <= num_element)
    {
        //sample_vec.clear();
        // bool last_batch = false;
        int existed_feature = 0;
        const int num_block = num_element / block_size;

        for (uint8_t i = 0; i < sample_start_pos.size(); i++)
        {
            sample_start_pos[i] = rand() % num_element;
        }
        std::sort(sample_start_pos.begin(), sample_start_pos.end());
        // printf("finish sort for block size %d\n", block_size);
        for (int i = 0; i < sample_start_pos.size(); i++)
        {
            std::vector<uint32_t> sample_vec(block_size);
            bool exceed_final_slot = false;
            int start_slot = data_mat[num_feature][sample_start_pos[i]], end_slot, sample_stop_pos, sample_vector_index = 0;
            // printf("initialize start slot\n");
            
            // determine the final sample position
            end_slot = data_mat[num_feature][sample_start_pos[i] + block_size - 1];
            sample_stop_pos = (sample_start_pos[i] + block_size - 1) % num_element;
            if (start_slot > end_slot)
            {
                exceed_final_slot = true;
            }
        
            // printf("initialize end slot\n");

            // first_slot
            int start_slot_first = cum_size[start_slot];
            // printf("initialize start_slot\n");
            for (int j = 0; j < cum_size[start_slot + 1] - sample_start_pos[i]; j++)
            {
                sample_vec[sample_vector_index ++] = cum_size[start_slot] + Next() % (cum_size[start_slot + 1] - cum_size[start_slot]);
            }
            // printf("sample start_slot\n");
            if (exceed_final_slot == false)
            {
                for (int j = start_slot + 1; j < end_slot; j++)
                {
                    uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                    if (slot_size == 1)
                    {
                        sample_vec[sample_vector_index ++] = cum_size[j];
                    }
                    else
                    {
                        for (uint8_t k = 0; k < slot_size; k++)
                        {
                            sample_vec[sample_vector_index ++] = cum_size[j] + Next() % slot_size;
                        }
                    }
                }
            }
            else
            {
                for (int j = start_slot + 1; j < cum_size.size(); j++)
                {
                    uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                    if (slot_size == 1)
                    {
                        sample_vec[sample_vector_index ++] = cum_size[j];
                    }
                    else
                    {
                        for (uint8_t k = 0; k < slot_size; k++)
                        {
                            sample_vec[sample_vector_index ++] = cum_size[j] + Next() % slot_size;
                        }
                    }
                }
                for (int j = 0; j < end_slot; j++)
                {
                    uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                    if (slot_size == 1)
                    {
                        sample_vec[sample_vector_index ++] = cum_size[j];
                    }
                    else
                    {
                        for (uint8_t k = 0; k < slot_size; k++)
                        {
                            sample_vec[sample_vector_index ++] = cum_size[j] + Next() % slot_size;
                        }
                    }
                }
            }

            // last slot
            uint32_t end_slot_last = cum_size[end_slot + 1];
            uint8_t end_slot_sample_size = sample_stop_pos - cum_size[end_slot] + 1;
            // printf("initialize end_slot\n");
            for (int j = 0; j < end_slot_sample_size; j++)
            {
                sample_vec[sample_vector_index ++] = cum_size[end_slot] + Next() % (cum_size[end_slot + 1] - cum_size[end_slot]);
            }
            // printf("exact block size = %d\n", sample_vector_index);
            // printf("sample end_slot\n");
            for (int j = 0; j < feature_vec.size(); j++)
            {
                if (feature_vec[j].ratio_ < (1 - error_bound))
                {
                    int feature_idx = feature_vec[j].idx_;
                    UpdateSingleBlock2(block_size, MoM_pair[feature_idx], i, data_mat[feature_idx], sample_vec, start_slot_first, end_slot_last, exceed_final_slot, precise_time);
                }
            }
            // printf("sample records\n");
        }

        for (int i = 0; i < feature_vec.size(); i++)
        {
            if (feature_vec[i].ratio_ < (1 - error_bound))
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                // printf("Feature %d\n", feature_idx);
                VarianceEst variance(block_size, MoM_pair[feature_idx], pf0);
                // printf("special idx %3d low %2.6f up %1.6f\n", feature_vec[i].idx_, variance.low_, variance.up_);
                feature_vec[i].Update(variance.low_, variance.up_, variance.sample_variance_);
            }
            // if (feature_vec[i].idx_ == 107 || feature_vec[i].idx_ == 1)
            // {
            //     printf("special idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            // }
        }

        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_low);
        double kth_low = feature_vec[topk - 1].low_;
        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
        sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
        double min_ratio = std::min_element(feature_vec.begin(), feature_vec.begin() + topk, cmp_ratio)->ratio_;
        printf("SampleSize %d ExistedFeature %d\n", block_size * num_var, existed_feature);

        // if (block_size * num_var >= 2000000)
        // {
        //     for (int i = 0; i < feature_vec.size(); ++i)
        //     {
        //         printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        //     }
        // // }
        // printf("\n");

        for (int i = 0; i < topk; i++)
        {
            printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        } 
        count += block_size * num_var * existed_feature;
        if (min_ratio >= (1 - error_bound))
        {
            printf("ApproximateVarianceTopk has been found.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }
            printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var, count, precise_time, median_time, sort_time);
            break;
        }
        else if (block_size * num_var > num_element / 2)
        {
            std::chrono::high_resolution_clock::time_point start, end; 
            std::chrono::duration<double> duration;
            start = std::chrono::high_resolution_clock::now();

            for (int i = 0; i < feature_vec.size(); i++)
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                double exact_variance = ExactVarianceScore(data_mat[feature_idx]);
                feature_vec[i].Update(exact_variance, exact_variance, exact_variance);
            }
            std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
            sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
            printf("ApproximateVarianceTopk has been found when elements are used up.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }

            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            sort_time += duration.count();

            printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", num_element, count, precise_time, median_time, sort_time);
            break;
        }
        else
        {
            block_size *= 2;
        }
        
        std::vector<Feature>::const_iterator iter = feature_vec.begin();
        while (iter != feature_vec.end()) 
        {
            if (iter->up_ < kth_low) 
            {
                iter = feature_vec.erase(iter);
            }
            else 
            {
                ++iter;
            }
        }
    } 
    // printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var / 2, count, precise_time, median_time, sort_time);
}

void VarianceTopkMomBlock3(const double ratio_num_block, const int topk, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, std::vector<int>& result_idx, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    long count = 0; 
    int block_size = initial_sample_size, num_slot = cum_size.size() - 1;
    double pf0 = pf * 1.0 / num_feature / ceil(log2(num_element * 1.0 / initial_sample_size)), sample_index_time = 0.0, access_data_time = 0.0, sort_time = 0.0;
    const int num_var = ceil(ratio_num_block * log(2.0 / pf0));
    // std::vector<std::vector<uint32_t>> sample_vec(num_var, std::vector<uint32_t>(0));
    std::vector<uint32_t> sample_vec;

    printf("Element %d Feature %d pf0 %1.15f 1_pf0 %f NumVar %d BlockSize %d Block3.5\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var, block_size);

    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    
    std::vector<Feature> feature_vec;
    std::vector<uint32_t> sample_start_pos(num_var, 0);
    std::vector<std::pair<std::vector<double>, std::vector<double>>> MoM_pair(num_feature, std::pair<std::vector<double>, std::vector<double>>(std::vector<double>(num_var), std::vector<double>(num_var)));
    for (uint8_t i = 0; i < num_feature; i++)
    {
        feature_vec.push_back(Feature(i, 0.0));
    }
    result_idx.clear();
    result_value.clear();
    while (!feature_vec.empty() && block_size * num_var <= num_element)
    {
        //sample_vec.clear();
        // bool last_batch = false;
        int existed_feature = 0;
        const int num_block = num_element / block_size;
        sample_vec.resize(block_size * num_var);

        std::chrono::high_resolution_clock::time_point start, end; 
        std::chrono::duration<double> duration;
        start = std::chrono::high_resolution_clock::now();

        for (uint8_t i = 0; i < sample_start_pos.size(); i++)
        {
            sample_start_pos[i] = rand() % num_element;
        }
        std::sort(sample_start_pos.begin(), sample_start_pos.end());
        // printf("finish sort for block size %d\n", block_size);
        for (int i = 0; i < sample_start_pos.size(); i++)
        {
            bool exceed_final_slot = false;
            int start_slot = data_mat[num_feature][sample_start_pos[i]], end_slot = data_mat[num_feature][sample_start_pos[i] + block_size - 1], sample_stop_pos = (sample_start_pos[i] + block_size - 1) % num_element, sample_vector_index = 0;
            // printf("initialize start slot\n");
            
            if (start_slot > end_slot)
            {
                exceed_final_slot = true;
            }

            // first_slot
            // int start_slot_first = cum_size[start_slot];
            // printf("initialize start_slot\n");
            uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i];
            for (uint8_t j = 0; j < start_slot_sample_size; j++)
            {
                sample_vec[i * block_size + sample_vector_index] = cum_size[start_slot] + Next() % (cum_size[start_slot + 1] - cum_size[start_slot]);
                sample_vector_index ++;
            }
            // printf("sample start_slot\n");
            if (exceed_final_slot == false)
            {
                for (uint32_t j = start_slot + 1; j < end_slot; j++)
                {
                    uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                    if (slot_size == 1)
                    {
                        // if (i * block_size + sample_vector_index >= sample_vec.size())
                        // {
                        //     printf("%d %d %d error!\n", i, block_size, sample_vector_index);
                        // }
                        sample_vec[i * block_size + sample_vector_index] = cum_size[j];
                        sample_vector_index ++;
                    }
                    else
                    {
                        for (uint8_t k = 0; k < slot_size; k++)
                        {
                            // if (i * block_size + sample_vector_index >= sample_vec.size())
                            // {
                            //     printf("%d %d %d error!\n", i, block_size, sample_vector_index);
                            // }
                            sample_vec[i * block_size + sample_vector_index] = cum_size[j] + Next() % slot_size;
                            sample_vector_index ++;
                        }
                    }
                }
            }
            else
            {
                for (int j = start_slot + 1; j < cum_size.size(); j++)
                {
                    uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                    if (slot_size == 1)
                    {
                        // if (i * block_size + sample_vector_index >= sample_vec.size())
                        // {
                        //     printf("%d %d %d error!\n", i, block_size, sample_vector_index);
                        // }
                        sample_vec[i * block_size + sample_vector_index] = cum_size[j];
                        sample_vector_index ++;
                    }
                    else
                    {
                        for (uint8_t k = 0; k < slot_size; k++)
                        {
                            // if (i * block_size + sample_vector_index >= sample_vec.size())
                            // {
                            //     printf("%d %d %d error!\n", i, block_size, sample_vector_index);
                            // }
                            sample_vec[i * block_size + sample_vector_index] = cum_size[j] + Next() % slot_size;
                            sample_vector_index ++;
                        }
                    }
                }
                for (int j = 0; j < end_slot; j++)
                {
                    uint8_t slot_size = cum_size[j + 1] - cum_size[j];
                    if (slot_size == 1)
                    {
                        // if (i * block_size + sample_vector_index >= sample_vec.size())
                        // {
                        //     printf("%d %d %d error!\n", i, block_size, sample_vector_index);
                        // }
                        sample_vec[i * block_size + sample_vector_index] = cum_size[j];
                        sample_vector_index ++;
                    }
                    else
                    {
                        for (uint8_t k = 0; k < slot_size; k++)
                        {
                            // if (i * block_size + sample_vector_index >= sample_vec.size())
                            // {
                            //     printf("%d %d %d error!\n", i, block_size, sample_vector_index);
                            // }
                            sample_vec[i * block_size + sample_vector_index] = cum_size[j] + Next() % slot_size;
                            sample_vector_index ++;
                        }
                    }
                }
            }

            // last slot
            // uint32_t end_slot_last = cum_size[end_slot + 1];
            uint8_t end_slot_sample_size = sample_stop_pos - cum_size[end_slot] + 1;
            // printf("initialize end_slot\n");
            for (int j = 0; j < end_slot_sample_size; j++)
            {
                // if (i * block_size + sample_vector_index >= sample_vec.size())
                // {
                //     printf("%d %d %d error!\n", i, block_size, sample_vector_index);
                // }
                sample_vec[i * block_size + sample_vector_index] = cum_size[end_slot] + Next() % (cum_size[end_slot + 1] - cum_size[end_slot]);
                sample_vector_index ++;
            }
            // printf("exact block size = %d\n", sample_vector_index);
            // printf("sample end_slot\n");
            // for (int j = 0; j < feature_vec.size(); j++)
            // {
            //     if (feature_vec[j].ratio_ < (1 - error_bound))
            //     {
            //         int feature_idx = feature_vec[j].idx_;
            //         UpdateSingleBlock2(block_size, MoM_pair[feature_idx], i, data_mat[feature_idx], sample_vec, start_slot_first, end_slot_last, exceed_final_slot, precise_time);
            //     }
            // }
            // printf("sample records\n");
        }

        end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		sample_index_time += duration.count();

        for (int i = 0; i < feature_vec.size(); i++)
        {
            if (feature_vec[i].ratio_ < (1 - error_bound))
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                // printf("Feature %d\n", feature_idx);
                VarianceEst variance(num_var, block_size, sample_vec, MoM_pair[feature_idx], data_mat[feature_idx], pf0, access_data_time);
                // printf("special idx %3d low %2.6f up %1.6f\n", feature_vec[i].idx_, variance.low_, variance.up_);
                feature_vec[i].Update(variance.low_, variance.up_, variance.sample_variance_);
            }
            // if (feature_vec[i].idx_ == 107 || feature_vec[i].idx_ == 1)
            // {
            //     printf("special idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            // }
        }

        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_low);
        double kth_low = feature_vec[topk - 1].low_;
        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
        sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
        double min_ratio = std::min_element(feature_vec.begin(), feature_vec.begin() + topk, cmp_ratio)->ratio_;
        printf("SampleSize %d ExistedFeature %d\n", block_size * num_var, existed_feature);

        // if (block_size * num_var >= 2000000)
        // {
        //     for (int i = 0; i < feature_vec.size(); ++i)
        //     {
        //         printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        //     }
        // // }
        // printf("\n");

        for (int i = 0; i < topk; i++)
        {
            printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        } 
        count += block_size * num_var * existed_feature;
        if (min_ratio >= (1 - error_bound))
        {
            printf("ApproximateVarianceTopk has been found.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }
            printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f\n", block_size * num_var, count, sample_index_time, access_data_time);
            break;
        }
        else if (block_size * num_var > num_element / 2)
        {
            std::chrono::high_resolution_clock::time_point start, end; 
            std::chrono::duration<double> duration;
            start = std::chrono::high_resolution_clock::now();

            for (int i = 0; i < feature_vec.size(); i++)
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                double exact_variance = ExactVarianceScore(data_mat[feature_idx]);
                feature_vec[i].Update(exact_variance, exact_variance, exact_variance);
            }
            std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
            sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
            printf("ApproximateVarianceTopk has been found when elements are used up.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }

            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            sort_time += duration.count();

            printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f SortTime %f\n", num_element, count, sample_index_time, access_data_time, sort_time);
            break;
        }
        else
        {
            block_size *= 2;
        }
        
        std::vector<Feature>::const_iterator iter = feature_vec.begin();
        while (iter != feature_vec.end()) 
        {
            if (iter->up_ < kth_low) 
            {
                iter = feature_vec.erase(iter);
            }
            else 
            {
                ++iter;
            }
        }
    } 
    // printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var / 2, count, precise_time, median_time, sort_time);
}

void VarianceTopkMomBlock4(const double ratio_num_block, const int topk, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, std::vector<int>& result_idx, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    long count = 0; 
    int block_size = initial_sample_size, num_slot = cum_size.size() - 1;
    double pf0 = pf * 1.0 / num_feature / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
    const int num_var = ceil(ratio_num_block * log(2.0 / pf0));
    // std::vector<std::vector<uint32_t>> sample_vec(num_var, std::vector<uint32_t>(0));
    // std::vector<uint32_t> sample_vec;

    printf("Element %d Feature %d pf0 %1.15f 1_pf0 %f NumVar %d BlockSize %d Block4\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var, block_size);

    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;
    
    std::vector<Feature> feature_vec;
    std::vector<int> sample_start_pos(num_var, 0);
    std::vector<std::pair<std::vector<double>, std::vector<double>>> MoM_pair(num_feature, std::pair<std::vector<double>, std::vector<double>>(std::vector<double>(num_var), std::vector<double>(num_var)));
    for (uint8_t i = 0; i < num_feature; i++)
    {
        feature_vec.push_back(Feature(i, 0.0));
    }
    result_idx.clear();
    result_value.clear();
    while (!feature_vec.empty() && block_size * num_var <= num_element)
    {
        //sample_vec.clear();
        // bool last_batch = false;
        int existed_feature = 0;
        const int num_block = num_element / block_size;
        // sample_vec.resize(block_size * num_var);
        for (uint8_t i = 0; i < sample_start_pos.size(); i++)
        {
            sample_start_pos[i] = rand() % num_element;
        }
        std::sort(sample_start_pos.begin(), sample_start_pos.end());
        // printf("finish sort for block size %d\n", block_size);
        
        for (int i = 0; i < feature_vec.size(); i++)
        {
            if (feature_vec[i].ratio_ < (1 - error_bound))
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                // printf("Feature %d\n", feature_idx);
                UpdateCumSeq(block_size, num_element,  sample_start_pos, MoM_pair[feature_idx], data_mat[feature_idx], data_mat[num_feature], cum_size, precise_time);
                VarianceEst variance(block_size, MoM_pair[feature_idx], pf0);
                // printf("special idx %3d low %2.6f up %1.6f\n", feature_vec[i].idx_, variance.low_, variance.up_);
                feature_vec[i].Update(variance.low_, variance.up_, variance.sample_variance_);
            }
            // if (feature_vec[i].idx_ == 107 || feature_vec[i].idx_ == 1)
            // {
            //     printf("special idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            // }
        }

        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_low);
        double kth_low = feature_vec[topk - 1].low_;
        std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
        sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
        double min_ratio = std::min_element(feature_vec.begin(), feature_vec.begin() + topk, cmp_ratio)->ratio_;
        printf("SampleSize %d ExistedFeature %d\n", block_size * num_var, existed_feature);

        // if (block_size * num_var >= 2000000)
        // {
        //     for (int i = 0; i < feature_vec.size(); ++i)
        //     {
        //         printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        //     }
        // // }
        // printf("\n");

        for (int i = 0; i < topk; i++)
        {
            printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
        } 
        count += block_size * num_var * existed_feature;
        if (min_ratio >= (1 - error_bound))
        {
            printf("ApproximateVarianceTopk has been found.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }
            printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var, count, precise_time, median_time, sort_time);
            break;
        }
        else if (block_size * num_var > num_element / 2)
        {
            std::chrono::high_resolution_clock::time_point start, end; 
            std::chrono::duration<double> duration;
            start = std::chrono::high_resolution_clock::now();

            for (int i = 0; i < feature_vec.size(); i++)
            {
                existed_feature ++;
                int feature_idx = feature_vec[i].idx_;
                double exact_variance = ExactVarianceScore(data_mat[feature_idx]);
                feature_vec[i].Update(exact_variance, exact_variance, exact_variance);
            }
            std::nth_element(feature_vec.begin(), feature_vec.begin() + topk - 1, feature_vec.end(), cmp_up);
            sort(feature_vec.begin(), feature_vec.begin() + topk, cmp_up);
            printf("ApproximateVarianceTopk has been found when elements are used up.\n");
            for (int i = 0; i < topk; i++)
            {
                result_idx.push_back(feature_vec[i].idx_);
                result_value.push_back(feature_vec[i].avg_);
                printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f ratio %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, feature_vec[i].ratio_);
            }

            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            sort_time += duration.count();

            printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", num_element, count, precise_time, median_time, sort_time);
            break;
        }
        else
        {
            block_size *= 2;
        }
        
        std::vector<Feature>::const_iterator iter = feature_vec.begin();
        while (iter != feature_vec.end()) 
        {
            if (iter->up_ < kth_low) 
            {
                iter = feature_vec.erase(iter);
            }
            else 
            {
                ++iter;
            }
        }
    } 
    // printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var / 2, count, precise_time, median_time, sort_time);
}

void EpsilonDeltaMomVariance(const double ratio_num_block, const double empirical_threshold, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size(), num_test_feature = test_feature.size();
    long count = 0; 
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
    int num_var;
    // if (metric == "variance")
    // {
    num_var = ceil(ratio_num_block * log(2.0 / pf0));
    // }
    // else if (metric == "average")
    // {
    //     num_var = ceil(ratio_num_block * log(1.0 / pf0));
    // }
    // else
    // {
    //     printf("Error Metric!\n");
    // }
    std::pair<std::vector<double>, std::vector<double>> MoM_pair{std::vector<double>(num_var), std::vector<double>(num_var)};
    std::vector<int> sample_start_pos(num_var, 0);
     
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;

    printf("Element %d NumTestFeature %d pf0 %1.15f 1_pf0 %f NumVar %d\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var);

    for (int j = 0; j < test_feature.size(); j++)
    {
        Feature single_feature(test_feature[j], 0.0);
        int block_size = initial_sample_size;
        // double max_diff_square = 1.0;
        
        while (block_size * num_var < num_element)
        {
            // bool last_batch = false;
            const int num_block = num_element / block_size;

            start = std::chrono::high_resolution_clock::now();
            // for (int i = 0; i < sample_vec.size(); i++)
            // {
            //     const int random_pos = rand() % num_element;
            //     if (random_pos >= block_size + num_element % block_size)
            //     {
            //         sample_vec[i] = rand() % (num_block - 1);
            //     }
            //     else
            //     {
            //         sample_vec[i] = num_block - 1;
            //         last_batch = true;
            //         printf("RandPos %d BlockSize %d Res %d Total %d\n", random_pos, block_size, int(num_element % block_size), int(block_size + num_element % block_size));
            //     }
            // }
            for (int i = 0; i < sample_start_pos.size(); i++)
            {
                sample_start_pos[i] = Next() % num_element;
            }
            // std::sort(sample_start_pos.begin(), sample_start_pos.end());
            
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            median_time += duration.count();

            UpdateCumSeq(block_size, num_element, sample_start_pos, MoM_pair, data_mat[test_feature[j]], data_mat[num_feature], cum_size, precise_time);
            VarianceEst variance(block_size, MoM_pair, pf0);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
            // else if (metric == "average")
            // {
            //     AverageEst average(sample_vec, block_size, MoM_pair.first, data_mat[test_feature[j]], pf0, precise_time, max_diff_square);
            //     single_feature.Update(average.low_, average.up_);
            // }
            
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", block_size * num_var, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
            count += block_size * num_var;
            if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
            {
                printf("ApproximateValue has been found.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var, count, precise_time, median_time, sort_time);
                break;
            }
            else if (single_feature.diff_ <= 2.0 * error_bound * empirical_threshold && single_feature.low_ < empirical_threshold)
            {
                printf("ApproximateValue has been found when the lower is smaller than the empirical threshold.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonDelta %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * empirical_threshold);
                printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var, count, precise_time, median_time, sort_time);
                break;
            }
            else if (block_size * num_var > num_element / 2)
            {
                start = std::chrono::high_resolution_clock::now();
                // if (metric == "variance")
                // {
                double exact_variance = ExactVarianceScore(data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);
                // }
                // else if (metric == "average")
                // {
                //     double exact_average = ExactAverageScore(data_mat[test_feature[j]]);
                //     single_feature.Update(exact_average, exact_average);
                // }
                result_value[j] = single_feature.est_;
                printf("ApproximateValue has been found when elements are used up.\n");
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], num_element, single_feature.low_, single_feature.up_, single_feature.avg_,  single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);

                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", num_element, count, precise_time, median_time, sort_time);
                break;
            }
            else
            {
                block_size *= 2;
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", i, result_value[i]);   
    }
}

void SelectAbsPredMomVariance(const double ratio_num_block, const PredPair& pred, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const double pf, const double abs_error, std::vector<double>& result_value)
{
    const uint8_t num_feature = data_mat.size() - 1, num_test_feature = test_feature.size();
    const uint32_t num_element = data_mat[0].size();
    uint64_t count = 0, num_total_scan = 0;
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), access_data_time = 0.0, median_time = 0.0, sort_time = 0.0;
    // select_pr = 0.0, pract_error = 0.0
    const uint8_t num_var = ceil(ratio_num_block * log(2.0 / pf0));
    
    std::pair<std::vector<double>, std::vector<double>> MoM_pair{std::vector<double>(num_var), std::vector<double>(num_var)};
    std::vector<uint32_t> sample_start_pos(num_var, 0);
    // std::vector<uint32_t> sample_vec;
     
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;

    printf("Element %d NumTestFeature %d pf0 %1.15f 1_pf0 %f NumVar %d Block2\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var);

    for (uint8_t j = 0; j < test_feature.size(); j++)
    {
        Feature single_feature(test_feature[j], 0.0);
        uint32_t block_size = initial_sample_size;
        while (block_size * num_var < num_element)
        {
            uint32_t num_scan = 0;
            start = std::chrono::high_resolution_clock::now();
            
            for (int i = 0; i < sample_start_pos.size(); i++)
            {
                sample_start_pos[i] = Next() % num_element;
            }
            // std::sort(sample_start_pos.begin(), sample_start_pos.end());
            
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            median_time += duration.count();
            
            UpdateCumSeqPredBlock4(block_size, sample_start_pos, MoM_pair, pred.value_, data_mat[pred.col_], data_mat[test_feature[j]], data_mat[num_feature], num_scan, cum_size, access_data_time);
            // select_pr = 1.0 * block_size * num_var / num_scan;

            VarianceEst variance(block_size, MoM_pair, pf0);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
            // pract_error = abs_error / sqrt(select_pr);
            
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2PractError %1.10f\n", block_size * num_var, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * abs_error);
            count += block_size * num_var;
            num_total_scan += num_scan;
            
            if (single_feature.diff_ <= 2.0 * abs_error)
            {
                printf("ApproximateValue has been found.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2PractError %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * abs_error);
                printf("SampleSize %d Count %lu TotalScan %lu SampleIndexTime %f AccessDataTime %f\n", block_size * num_var, count, num_total_scan, median_time, access_data_time);
                break;
            }
            else if (num_scan > num_element / 2)
            {
                start = std::chrono::high_resolution_clock::now();
                const double exact_variance = ExactPredVarianceScore(pred.value_, data_mat[pred.col_], data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);
                result_value[j] = single_feature.est_;
                printf("ApproximateValue has been found when elements are used up.\n");
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2PractError %1.10f\n", test_feature[j], num_element, single_feature.low_, single_feature.up_, single_feature.avg_,  single_feature.est_, single_feature.diff_, 2.0 * abs_error);

                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %lu TotalScan %lu SampleIndexTime %f AccessDataTime %f SortTime %f\n", num_element, count, num_total_scan, median_time, access_data_time, sort_time);
                break;
            }
            else
            {
                block_size *= 2;
                // printf("Double\n");
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", test_feature[i], result_value[i]);   
    }
}

void SelectAbsPredMomVarianceStat(const double ratio_num_block, const PredPair& pred, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const double pf, const double abs_error, std::vector<double>& result_value)
{
    const uint8_t num_feature = data_mat.size() - 1, num_test_feature = test_feature.size();
    const uint32_t num_element = data_mat[0].size(), num_slot = cum_size.size() - 1;
    uint64_t count = 0, num_total_scan = 0;
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), access_data_time = 0.0, median_time = 0.0, sort_time = 0.0;
    // select_pr = 0.0, pract_error = 0.0
    const uint8_t num_var = ceil(ratio_num_block * log(2.0 / pf0));
    
    std::pair<std::vector<double>, std::vector<double>> MoM_pair{std::vector<double>(num_var), std::vector<double>(num_var)};
    std::vector<uint32_t> sample_start_pos(num_var, 0);
    // std::vector<uint32_t> sample_vec;
     
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;

    printf("Element %d NumTestFeature %d pf0 %1.15f 1_pf0 %f NumVar %d Block2\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var);

    for (uint8_t j = 0; j < test_feature.size(); j++)
    {
        Feature single_feature(test_feature[j], 0.0);
        uint32_t block_size = initial_sample_size;
        while (block_size * num_var < num_element)
        {
            uint32_t num_scan = 0;
            start = std::chrono::high_resolution_clock::now();
            
            for (int i = 0; i < sample_start_pos.size(); i++)
            {
                sample_start_pos[i] = Next() % num_element;
            }
            // std::sort(sample_start_pos.begin(), sample_start_pos.end());
            
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            median_time += duration.count();
            
            // UpdateCumSeqPredBlock4(block_size, sample_start_pos, MoM_pair, pred.value_, data_mat[pred.col_], data_mat[test_feature[j]], data_mat[num_feature], num_scan, cum_size, access_data_time);
            // // select_pr = 1.0 * block_size * num_var / num_scan;
            
            for (uint8_t i = 0; i < sample_start_pos.size(); i++)
            {
                int start_slot = data_mat[num_feature][sample_start_pos[i]];
                double sum = 0.0, sum_square = 0.0;
                uint8_t start_slot_sample_size = cum_size[start_slot + 1] - sample_start_pos[i], slot_size = cum_size[start_slot + 1] - cum_size[start_slot];
                num_scan += start_slot_sample_size;
                uint32_t count = 0, current_slot = start_slot, sign;
                for (uint8_t k = 0; k < start_slot_sample_size; k++)
                {
                    uint32_t sample_index = cum_size[start_slot] + NextShort() % slot_size;
                    sign = my_sign(abs(data_mat[pred.col_][sample_index] - pred.value_) - 1e-12);
                    sum += sign * data_mat[test_feature[j]][sample_index];
                    sum_square += sign * data_mat[test_feature[j]][sample_index] * data_mat[test_feature[j]][sample_index];
                    count += sign;
                }
                // bool flag = true;
                while (count < block_size)
                {
                    current_slot += 1;
                    if (unlikely(current_slot >= num_slot))
                    {
                        current_slot -= num_slot;
                    }
                    slot_size = cum_size[current_slot + 1] - cum_size[current_slot];
                    num_scan += slot_size;
                    if (likely(slot_size == 1))
                    {
                        uint32_t sample_index = cum_size[current_slot];
                        sign = my_sign(abs(data_mat[pred.col_][sample_index] - pred.value_) - 1e-12);
                        sum += sign * data_mat[test_feature[j]][sample_index];
                        sum_square += sign * data_mat[test_feature[j]][sample_index] * data_mat[test_feature[j]][sample_index];
                        count += sign;
                    }
                    else
                    {
                        uint32_t slot_rand = Next();
                        for (uint8_t k = 0; k < slot_size; k++)
                        {
                            uint32_t sample_index = cum_size[current_slot] + (slot_rand & 0xff) % slot_size;
                            slot_rand >>= 8;
                            // printf("sample_index %d current_slot %d start %d\n", sample_index, current_slot, cum_size[current_slot]);
                            sign = my_sign(abs(data_mat[pred.col_][sample_index] - pred.value_) - 1e-12);
                            sum += sign * data_mat[test_feature[j]][sample_index];
                            sum_square += sign * data_mat[test_feature[j]][sample_index] * data_mat[test_feature[j]][sample_index];
                            count += sign;
                            if (unlikely(count >= block_size))
                            {
                                k = slot_size;
                            }
                        }
                    }
                }
                MoM_pair.first[i] = sum / block_size;
                MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
            }

            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            median_time += duration.count();

            VarianceEst variance(block_size, MoM_pair, pf0);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
            // pract_error = abs_error / sqrt(select_pr);
            
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2PractError %1.10f NumofScan %d\n", block_size * num_var, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * abs_error, num_scan);
            count += block_size * num_var;
            num_total_scan += num_scan;
            
            if (single_feature.diff_ <= 2.0 * abs_error)
            {
                printf("ApproximateValue has been found.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2PractError %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * abs_error);
                printf("SampleSize %d Count %lu TotalScan %lu SampleIndexTime %f AccessDataTime %f\n", block_size * num_var, count, num_total_scan, median_time, access_data_time);
                break;
            }
            else if (num_scan > num_element / 2)
            {
                start = std::chrono::high_resolution_clock::now();
                const double exact_variance = ExactPredVarianceScore(pred.value_, data_mat[pred.col_], data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);
                result_value[j] = single_feature.est_;
                printf("ApproximateValue has been found when elements are used up.\n");
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2PractError %1.10f\n", test_feature[j], num_element, single_feature.low_, single_feature.up_, single_feature.avg_,  single_feature.est_, single_feature.diff_, 2.0 * abs_error);

                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %lu TotalScan %lu SampleIndexTime %f AccessDataTime %f SortTime %f\n", num_element, count, num_total_scan, median_time, access_data_time, sort_time);
                break;
            }
            else
            {
                block_size *= 2;
                // printf("Double\n");
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", test_feature[i], result_value[i]);   
    }
}

void EpsilonDeltaPredMomVariance(const double ratio_num_block, const PredPair& pred, const double empirical_threshold, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const uint32_t num_element = data_mat[0].size(), num_test_feature = test_feature.size();
    const uint8_t num_feature = data_mat.size() - 1; 
    long count = 0; 
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), access_data_time = 0.0, median_time = 0.0, sort_time = 0.0;
    const uint8_t num_var = ceil(ratio_num_block * log(2.0 / pf0));
    uint32_t num_total_scan = 0;
    
    std::pair<std::vector<double>, std::vector<double>> MoM_pair{std::vector<double>(num_var), std::vector<double>(num_var)};
    std::vector<uint32_t> sample_start_pos(num_var, 0);
    // std::vector<uint32_t> sample_vec;
     
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;

    printf("Element %d NumTestFeature %d pf0 %1.15f 1_pf0 %f NumVar %d Block2\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var);

    for (int j = 0; j < test_feature.size(); j++)
    {
        Feature single_feature(test_feature[j], 0.0);
        int block_size = initial_sample_size;
        uint32_t num_scan = 0;
        
        // double max_diff_square = 1.0;
        
        while (block_size * num_var < num_element)
        {
            // bool last_batch = false;
            // printf("num_var %d block_size %d mul %d\n", num_var, block_size, num_var * block_size);
            // printf("finish resize\n");
            const int num_block = num_element / block_size;

            start = std::chrono::high_resolution_clock::now();
            
            for (int i = 0; i < sample_start_pos.size(); i++)
            {
                sample_start_pos[i] = Next() % num_element;
            }
            // std::sort(sample_start_pos.begin(), sample_start_pos.end());
            
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            median_time += duration.count();
            
            UpdateCumSeqPred(block_size, sample_start_pos, MoM_pair, pred.value_, data_mat[pred.col_], data_mat[test_feature[j]], data_mat[num_feature], num_scan, cum_size, access_data_time);
            // UpdateCumSeqPredBlock2(block_size, sample_start_pos, sample_vec, MoM_pair, pred.value_, data_mat[pred.col_], data_mat[test_feature[j]], data_mat[num_feature], cum_size, precise_time, median_time);
            VarianceEst variance(block_size, MoM_pair, pf0);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
            // else if (metric == "average")
            // {
            //     AverageEst average(sample_vec, block_size, MoM_pair.first, data_mat[test_feature[j]], pf0, precise_time, max_diff_square);
            //     single_feature.Update(average.low_, average.up_);
            // }
            
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", block_size * num_var, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
            count += block_size * num_var;
            num_total_scan += num_scan;
            if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
            {
                printf("ApproximateValue has been found.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                printf("SampleSize %d Count %ld TotalScan %u SampleIndexTime %f AccessDataTime %f\n", block_size * num_var, count, num_total_scan, median_time, access_data_time);
                break;
            }
            else if (single_feature.diff_ <= 2.0 * error_bound * empirical_threshold && single_feature.low_ < empirical_threshold)
            {
                printf("ApproximateValue has been found when the lower is smaller than the empirical threshold.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonDelta %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * empirical_threshold);
                printf("SampleSize %d Count %ld TotalScan %u SampleIndexTime %f AccessDataTime %f\n", block_size * num_var, count, num_total_scan, median_time, access_data_time);
                break;
            }
            else if (num_scan > num_element / 2)
            {
                start = std::chrono::high_resolution_clock::now();
                // if (metric == "variance")
                // {
                double exact_variance = ExactPredVarianceScore(pred.value_, data_mat[pred.col_], data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);
                // }
                // else if (metric == "average")
                // {
                //     double exact_average = ExactAverageScore(data_mat[test_feature[j]]);
                //     single_feature.Update(exact_average, exact_average);
                // }
                result_value[j] = single_feature.est_;
                printf("ApproximateValue has been found when elements are used up.\n");
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], num_element, single_feature.low_, single_feature.up_, single_feature.avg_,  single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);

                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %ld TotalScan %u SampleIndexTime %f AccessDataTime %f SortTime %f\n", num_element, count, num_total_scan, median_time, access_data_time, sort_time);
                break;
            }
            else
            {
                block_size *= 2;
                // printf("Double\n");
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", test_feature[i], result_value[i]);   
    }
}

void EpsilonDeltaPredMomVarianceBlock2(const double ratio_num_block, const PredPair& pred, const double empirical_threshold, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size(), num_test_feature = test_feature.size();
    long count = 0; 
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), access_data_time = 0.0, median_time = 0.0, sort_time = 0.0;
    uint64_t num_total_scan = 0;
    int num_var = ceil(ratio_num_block * log(2.0 / pf0));
    
    std::pair<std::vector<double>, std::vector<double>> MoM_pair{std::vector<double>(num_var), std::vector<double>(num_var)};
    std::vector<int> sample_start_pos(num_var, 0);
    std::vector<uint32_t> sample_vec;
     
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;

    printf("Element %d NumTestFeature %d pf0 %1.15f 1_pf0 %f NumVar %d Block2\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var);

    for (int j = 0; j < test_feature.size(); j++)
    {
        uint32_t num_scan = 0; 
        Feature single_feature(test_feature[j], 0.0);
        int block_size = initial_sample_size;
        
        // double max_diff_square = 1.0;
        
        while (block_size * num_var < num_element)
        {
            // bool last_batch = false;
            // printf("num_var %d block_size %d mul %d\n", num_var, block_size, num_var * block_size);
            // printf("finish resize\n");
            const int num_block = num_element / block_size;

            // start = std::chrono::high_resolution_clock::now();
            
            for (int i = 0; i < sample_start_pos.size(); i++)
            {
                sample_start_pos[i] = Next() % num_element;
            }
            // std::sort(sample_start_pos.begin(), sample_start_pos.end());
            
            // end = std::chrono::high_resolution_clock::now();
            // duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            // median_time += duration.count();
            
            UpdateCumSeqPredBlock2(block_size, sample_start_pos, sample_vec, MoM_pair, pred.value_, data_mat[pred.col_], data_mat[test_feature[j]], data_mat[num_feature], cum_size, num_scan, median_time, access_data_time);
            // UpdateCumSeqPredBlock2(block_size, sample_start_pos, sample_vec, MoM_pair, pred.value_, data_mat[pred.col_], data_mat[test_feature[j]], data_mat[num_feature], cum_size, precise_time, median_time);
            VarianceEst variance(block_size, MoM_pair, pf0);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
            // else if (metric == "average")
            // {
            //     AverageEst average(sample_vec, block_size, MoM_pair.first, data_mat[test_feature[j]], pf0, precise_time, max_diff_square);
            //     single_feature.Update(average.low_, average.up_);
            // }
            
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", block_size * num_var, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
            count += block_size * num_var;
            num_total_scan += num_scan;
            if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
            {
                printf("ApproximateValue has been found.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                printf("SampleSize %d Count %ld TotalScan %lu SampleIndexTime %f AccessDataTime %f\n", block_size * num_var, count, num_total_scan, median_time, access_data_time);
                break;
            }
            else if (single_feature.diff_ <= 2.0 * error_bound * empirical_threshold && single_feature.low_ < empirical_threshold)
            {
                printf("ApproximateValue has been found when the lower is smaller than the empirical threshold.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonDelta %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * empirical_threshold);
                printf("SampleSize %d Count %ld TotalScan %lu SampleIndexTime %f AccessDataTime %f\n", block_size * num_var, count, num_total_scan, median_time, access_data_time);
                break;
            }
            else if (block_size * num_var > num_element / 2)
            {
                start = std::chrono::high_resolution_clock::now();
                // if (metric == "variance")
                // {
                double exact_variance = ExactVarianceScore(data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);
                // }
                // else if (metric == "average")
                // {
                //     double exact_average = ExactAverageScore(data_mat[test_feature[j]]);
                //     single_feature.Update(exact_average, exact_average);
                // }
                result_value[j] = single_feature.est_;
                printf("ApproximateValue has been found when elements are used up.\n");
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], num_element, single_feature.low_, single_feature.up_, single_feature.avg_,  single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);

                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %ld TotalScan %lu SampleIndexTime %f AccessDataTime %f SortTime %f\n", num_element, count, num_total_scan, median_time, access_data_time, sort_time);
                break;
            }
            else
            {
                block_size *= 2;
                // printf("Double\n");
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", test_feature[i], result_value[i]);   
    }
}

void EpsilonDeltaPredMomVarianceBlock3(const double ratio_num_block, const PredPair& pred, const double empirical_threshold, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size(), num_test_feature = test_feature.size();
    long count = 0; 
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), access_data_time = 0.0, median_time = 0.0, sort_time = 0.0;
    int num_var = ceil(ratio_num_block * log(2.0 / pf0));
    
    std::pair<std::vector<double>, std::vector<double>> MoM_pair{std::vector<double>(num_var), std::vector<double>(num_var)};
    std::vector<int> sample_start_pos(num_var, 0);
    std::vector<uint32_t> sample_vec;
     
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;

    printf("Element %d NumTestFeature %d pf0 %1.15f 1_pf0 %f NumVar %d Block2\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var);

    for (int j = 0; j < test_feature.size(); j++)
    {
        Feature single_feature(test_feature[j], 0.0);
        int block_size = initial_sample_size;
        
        // double max_diff_square = 1.0;
        
        while (block_size * num_var < num_element)
        {
            // bool last_batch = false;
            // printf("num_var %d block_size %d mul %d\n", num_var, block_size, num_var * block_size);
            // printf("finish resize\n");
            const int num_block = num_element / block_size;

            start = std::chrono::high_resolution_clock::now();
            
            for (int i = 0; i < sample_start_pos.size(); i++)
            {
                sample_start_pos[i] = Next() % num_element;
            }
            // std::sort(sample_start_pos.begin(), sample_start_pos.end());
            
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            median_time += duration.count();
            
            UpdateCumSeqPredBlock3(block_size, sample_start_pos, MoM_pair, pred.value_, data_mat[pred.col_], data_mat[test_feature[j]], data_mat[num_feature], cum_size, access_data_time);
            // UpdateCumSeqPredBlock2(block_size, sample_start_pos, sample_vec, MoM_pair, pred.value_, data_mat[pred.col_], data_mat[test_feature[j]], data_mat[num_feature], cum_size, precise_time, median_time);
            VarianceEst variance(block_size, MoM_pair, pf0);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
            // else if (metric == "average")
            // {
            //     AverageEst average(sample_vec, block_size, MoM_pair.first, data_mat[test_feature[j]], pf0, precise_time, max_diff_square);
            //     single_feature.Update(average.low_, average.up_);
            // }
            
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", block_size * num_var, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
            count += block_size * num_var;
            if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
            {
                printf("ApproximateValue has been found.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f\n", block_size * num_var, count, median_time, access_data_time);
                break;
            }
            else if (single_feature.diff_ <= 2.0 * error_bound * empirical_threshold && single_feature.low_ < empirical_threshold)
            {
                printf("ApproximateValue has been found when the lower is smaller than the empirical threshold.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonDelta %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * empirical_threshold);
                printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f\n", block_size * num_var, count, median_time, access_data_time);
                break;
            }
            else if (block_size * num_var > num_element / 2)
            {
                start = std::chrono::high_resolution_clock::now();
                // if (metric == "variance")
                // {
                double exact_variance = ExactVarianceScore(data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);
                // }
                // else if (metric == "average")
                // {
                //     double exact_average = ExactAverageScore(data_mat[test_feature[j]]);
                //     single_feature.Update(exact_average, exact_average);
                // }
                result_value[j] = single_feature.est_;
                printf("ApproximateValue has been found when elements are used up.\n");
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], num_element, single_feature.low_, single_feature.up_, single_feature.avg_,  single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);

                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %ld SampleIndexTime %f AccessDataTime %f SortTime %f\n", num_element, count, median_time, access_data_time, sort_time);
                break;
            }
            else
            {
                block_size *= 2;
                // printf("Double\n");
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", test_feature[i], result_value[i]);   
    }
}

void EpsilonDeltaMomVarianceBlock2(const double ratio_num_block, const double empirical_threshold, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const double pf, const double error_bound, std::vector<double>& result_value)
{
    const int num_feature = data_mat.size() - 1, num_element = data_mat[0].size(), num_test_feature = test_feature.size();
    long count = 0; 
    if (result_value.size() != test_feature.size())
    {
        result_value.resize(num_test_feature);
    }
    double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
    int num_var;
    // if (metric == "variance")
    // {
    num_var = ceil(ratio_num_block * log(2.0 / pf0));
    // }
    // else if (metric == "average")
    // {
    //     num_var = ceil(ratio_num_block * log(1.0 / pf0));
    // }
    // else
    // {
    //     printf("Error Metric!\n");
    // }
    std::pair<std::vector<double>, std::vector<double>> MoM_pair{std::vector<double>(num_var), std::vector<double>(num_var)};
    std::vector<int> sample_start_pos(num_var, 0);
     
    std::chrono::high_resolution_clock::time_point start, end; 
    std::chrono::duration<double> duration;

    printf("Element %d NumTestFeature %d pf0 %1.15f 1_pf0 %f NumVar %d Block2\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var);

    for (int j = 0; j < test_feature.size(); j++)
    {
        Feature single_feature(test_feature[j], 0.0);
        int block_size = initial_sample_size;
        std::vector<uint32_t> sample_vec;
        // double max_diff_square = 1.0;
        
        while (block_size * num_var < num_element)
        {
            // bool last_batch = false;
            const int num_block = num_element / block_size;

            start = std::chrono::high_resolution_clock::now();
            // for (int i = 0; i < sample_vec.size(); i++)
            // {
            //     const int random_pos = rand() % num_element;
            //     if (random_pos >= block_size + num_element % block_size)
            //     {
            //         sample_vec[i] = rand() % (num_block - 1);
            //     }
            //     else
            //     {
            //         sample_vec[i] = num_block - 1;
            //         last_batch = true;
            //         printf("RandPos %d BlockSize %d Res %d Total %d\n", random_pos, block_size, int(num_element % block_size), int(block_size + num_element % block_size));
            //     }
            // }
            for (int i = 0; i < sample_start_pos.size(); i++)
            {
                sample_start_pos[i] = Next() % num_element;
            }
            // std::sort(sample_start_pos.begin(), sample_start_pos.end());
            
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            median_time += duration.count();

            SampleBlockRandomNum(sample_vec, num_element, num_var, block_size, sample_start_pos, data_mat[num_feature], cum_size, precise_time);
            VarianceEst variance(num_var, block_size, sample_vec, MoM_pair, data_mat[test_feature[j]], pf0, precise_time);
            single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
            // else if (metric == "average")
            // {
            //     AverageEst average(sample_vec, block_size, MoM_pair.first, data_mat[test_feature[j]], pf0, precise_time, max_diff_square);
            //     single_feature.Update(average.low_, average.up_);
            // }
            
            printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", block_size * num_var, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
            count += block_size * num_var;
            if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
            {
                printf("ApproximateValue has been found.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
                printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var, count, precise_time, median_time, sort_time);
                break;
            }
            else if (single_feature.diff_ <= 2.0 * error_bound * empirical_threshold && single_feature.low_ < empirical_threshold)
            {
                printf("ApproximateValue has been found when the lower is smaller than the empirical threshold.\n");
                result_value[j] = single_feature.est_;
                printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonDelta %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * empirical_threshold);
                printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var, count, precise_time, median_time, sort_time);
                break;
            }
            else if (block_size * num_var > num_element / 2)
            {
                start = std::chrono::high_resolution_clock::now();
                // if (metric == "variance")
                // {
                double exact_variance = ExactVarianceScore(data_mat[test_feature[j]]);
                single_feature.Update(exact_variance, exact_variance, exact_variance);
                // }
                // else if (metric == "average")
                // {
                //     double exact_average = ExactAverageScore(data_mat[test_feature[j]]);
                //     single_feature.Update(exact_average, exact_average);
                // }
                result_value[j] = single_feature.est_;
                printf("ApproximateValue has been found when elements are used up.\n");
                printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], num_element, single_feature.low_, single_feature.up_, single_feature.avg_,  single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);

                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                sort_time += duration.count();

                printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", num_element, count, precise_time, median_time, sort_time);
                break;
            }
            else
            {
                block_size *= 2;
            }
        } 
    }
    printf("\n");
    for (int i = 0; i < result_value.size(); i++)
    {
        printf("idx %d %1.12f\n", i, result_value[i]);   
    }
}

// void MoMSimpleMetric(const double& ratio_num_block, const std::string& metric, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const double pf, const double error_bound, std::vector<double>& result_value)
// {
//     const int num_feature = data_mat.size(), num_element = data_mat[0].size(), num_test_feature = test_feature.size();
//     long count = 0; 
//     if (result_value.size() != test_feature.size())
//     {
//         result_value.resize(num_test_feature);
//     }
//     double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
//     int num_var;
//     if (metric == "variance")
//     {
//         num_var = ceil(ratio_num_block * log(2.0 / pf0));
//     }
//     else if (metric == "average")
//     {
//         num_var = ceil(ratio_num_block * log(1.0 / pf0));
//     }
//     else
//     {
//         printf("Error Metric!\n");
//     }
     
//     std::chrono::high_resolution_clock::time_point start, end; 
//     std::chrono::duration<double> duration;

//     printf("Element %d NumTestFeature %d pf0 %1.15f 1_pf0 %f NumVar %d\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var);

//     for (int j = 0; j < test_feature.size(); j++)
//     {
//         Feature single_feature(test_feature[j], 0.0);
//         int block_size = initial_sample_size;
//         // double max_diff_square = 1.0;
//         std::pair<std::vector<double>, std::vector<double>> MoM_pair{std::vector<double>(num_var), std::vector<double>(num_var)};
//         std::vector<int> sample_vec(num_var, 0);
//         while (block_size * num_var < num_element)
//         {
//             // bool last_batch = false;
//             const int num_block = num_element / block_size;

//             start = std::chrono::high_resolution_clock::now();
//             // for (int i = 0; i < sample_vec.size(); i++)
//             // {
//             //     const int random_pos = rand() % num_element;
//             //     if (random_pos >= block_size + num_element % block_size)
//             //     {
//             //         sample_vec[i] = rand() % (num_block - 1);
//             //     }
//             //     else
//             //     {
//             //         sample_vec[i] = num_block - 1;
//             //         last_batch = true;
//             //         printf("RandPos %d BlockSize %d Res %d Total %d\n", random_pos, block_size, int(num_element % block_size), int(block_size + num_element % block_size));
//             //     }
//             // }
//             for (int i = 0; i < sample_vec.size(); i++)
//             {
//                 sample_vec[i] = rand() % num_element;
//             }
//             std::sort(sample_vec.begin(), sample_vec.end());
            
//             end = std::chrono::high_resolution_clock::now();
//             duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//             median_time += duration.count();

//             // std::sort(sample_vec.begin(), sample_vec.end());
//             if (metric == "variance")
//             {
//                 VarianceEst variance(sample_vec, block_size, MoM_pair, data_mat[test_feature[j]], pf0, precise_time, median_time);
//                 single_feature.Update(variance.low_, variance.up_, variance.sample_variance_);
//             }
//             // else if (metric == "average")
//             // {
//             //     AverageEst average(sample_vec, block_size, MoM_pair.first, data_mat[test_feature[j]], pf0, precise_time, max_diff_square);
//             //     single_feature.Update(average.low_, average.up_);
//             // }
            
//             printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", block_size * num_var, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
//             count += block_size * num_var;
//             if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
//             {
//                 printf("ApproximateValue has been found.\n");
//                 result_value[j] = single_feature.est_;
//                 printf("idx %3d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
//                 printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var, count, precise_time, median_time, sort_time);
//                 break;
//             }
//             else if (block_size * num_var > num_element / 2)
//             {
//                 start = std::chrono::high_resolution_clock::now();
//                 if (metric == "variance")
//                 {
//                     double exact_variance = ExactVarianceScore(data_mat[test_feature[j]]);
//                     single_feature.Update(exact_variance, exact_variance, exact_variance);
//                 }
//                 // else if (metric == "average")
//                 // {
//                 //     double exact_average = ExactAverageScore(data_mat[test_feature[j]]);
//                 //     single_feature.Update(exact_average, exact_average);
//                 // }
//                 result_value[j] = single_feature.est_;
//                 printf("ApproximateValue has been found when elements are used up.\n");
//                 printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f est %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], num_element, single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.est_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);

//                 end = std::chrono::high_resolution_clock::now();
//                 duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//                 sort_time += duration.count();

//                 printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", num_element, count, precise_time, median_time, sort_time);
//                 break;
//             }
//             else
//             {
//                 block_size *= 2;
//             }
//         } 
//     }
//     // printf("\n");
//     // for (int i = 0; i < result_value.size(); i++)
//     // {
//     //     printf("idx %d %1.12f\n", i, result_value[i]);   
//     // }
// }

// void MoMSimpleMetric(const double& ratio_num_block, const std::string& metric, const std::vector<int>& test_feature, const int initial_sample_size, const std::vector<std::vector<double>>& data_mat, const double pf, const double error_bound, std::vector<double>& result_value)
// {
//     const int num_feature = data_mat.size(), num_element = data_mat[0].size(), num_test_feature = test_feature.size();
//     long count = 0; 
//     if (result_value.size() != test_feature.size())
//     {
//         result_value.resize(num_test_feature);
//     }
//     double pf0 = pf * 1.0 / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0;
//     int num_var;
//     if (metric == "variance")
//     {
//         num_var = ceil(ratio_num_block * log(2.0 / pf0));
//     }
//     else if (metric == "average")
//     {
//         num_var = ceil(ratio_num_block * log(1.0 / pf0));
//     }
//     else
//     {
//         printf("Error Metric!\n");
//     }
     
//     std::chrono::high_resolution_clock::time_point start, end; 
//     std::chrono::duration<double> duration;

//     printf("Element %d NumTestFeature %d pf0 %1.15f 1_pf0 %f NumVar %d\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var);

//     for (int j = 0; j < test_feature.size(); j++)
//     {
//         Feature single_feature(test_feature[j], 0.0, 0.0);
//         int block_size = initial_sample_size;
//         double max_diff_square = 1.0;
//         std::pair<std::vector<double>, std::vector<double>> MoM_pair{std::vector<double>(num_var), std::vector<double>(num_var)};
//         std::vector<int> sample_vec(num_var, 0);
//         while (block_size * num_var < num_element)
//         {
//             // bool last_batch = false;
//             const int num_block = num_element / block_size;

//             start = std::chrono::high_resolution_clock::now();
//             // for (int i = 0; i < sample_vec.size(); i++)
//             // {
//             //     const int random_pos = rand() % num_element;
//             //     if (random_pos >= block_size + num_element % block_size)
//             //     {
//             //         sample_vec[i] = rand() % (num_block - 1);
//             //     }
//             //     else
//             //     {
//             //         sample_vec[i] = num_block - 1;
//             //         last_batch = true;
//             //         printf("RandPos %d BlockSize %d Res %d Total %d\n", random_pos, block_size, int(num_element % block_size), int(block_size + num_element % block_size));
//             //     }
//             // }
//             for (int i = 0; i < sample_vec.size(); i++)
//             {
//                 sample_vec[i] = rand() % num_element;
//             }
//             std::sort(sample_vec.begin(), sample_vec.end());
            
//             end = std::chrono::high_resolution_clock::now();
//             duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//             median_time += duration.count();

//             // std::sort(sample_vec.begin(), sample_vec.end());
//             if (metric == "variance")
//             {
//                 VarianceEst variance(sample_vec, block_size, MoM_pair, data_mat[test_feature[j]], pf0, precise_time, median_time);
//                 single_feature.Update(variance.low_, variance.up_);
//             }
//             else if (metric == "average")
//             {
//                 AverageEst average(sample_vec, block_size, MoM_pair.first, data_mat[test_feature[j]], pf0, precise_time, max_diff_square);
//                 single_feature.Update(average.low_, average.up_);
//             }
            
//             printf("SampleSize %d idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f max_diff_square %f\n", block_size * num_var, test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_, max_diff_square);
//             count += block_size * num_var;
//             if (single_feature.diff_ <= 2.0 * error_bound * single_feature.low_)
//             {
//                 printf("ApproximateValue has been found.\n");
//                 result_value[j] = single_feature.avg_;
//                 printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);
//                 printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", block_size * num_var, count, precise_time, median_time, sort_time);
//                 break;
//             }
//             else if (block_size * num_var > num_element / 2)
//             {
//                 start = std::chrono::high_resolution_clock::now();
//                 if (metric == "variance")
//                 {
//                     double exact_variance = ExactVarianceScore(data_mat[test_feature[j]]);
//                     single_feature.Update(exact_variance, exact_variance);
//                 }
//                 else if (metric == "average")
//                 {
//                     double exact_average = ExactAverageScore(data_mat[test_feature[j]]);
//                     single_feature.Update(exact_average, exact_average);
//                 }
//                 result_value[j] = single_feature.avg_;
//                 printf("ApproximateValue has been found when elements are used up.\n");
//                 printf("idx %3d num_element %d low %1.6f up %1.6f avg %1.6f diff %1.10f 2epsilonLow %1.10f\n", test_feature[j], num_element, single_feature.low_, single_feature.up_, single_feature.avg_, single_feature.diff_, 2.0 * error_bound * single_feature.low_);

//                 end = std::chrono::high_resolution_clock::now();
//                 duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//                 sort_time += duration.count();

//                 printf("SampleSize %d Count %ld PreciseTime %f MedianTime %f SortTime %f\n", num_element, count, precise_time, median_time, sort_time);
//                 break;
//             }
//             else
//             {
//                 block_size *= 2;
//             }
//         } 
//     }
//     // printf("\n");
//     // for (int i = 0; i < result_value.size(); i++)
//     // {
//     //     printf("idx %d %1.12f\n", i, result_value[i]);   
//     // }
// }

// void ApproxVarianceFilter(const std::string mode, const double& threshold, const int& initial_sample_size, const std::vector<std::vector<double>>& data_mat, std::vector<int>& result_idx, const double pf, const double error_bound, std::vector<std::pair<int, double>>& result_pair)
// {
//     const int num_feature = data_mat.size(), num_element = data_mat[0].size();
//     long count = 0;
//     SampleSizePair sample_size_pair(0, initial_sample_size);

//     double pf0 = pf * 1.0 / num_feature / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, sort_time = 0.0;
//     std::vector<std::pair<double, double>> cumulate_pair(num_feature, std::pair<double, double>(0.0, 0.0));
//     std::vector<Feature> feature_vec;
//     std::vector<int> sample_vec;

//     for (int i = 0; i < num_feature; i++)
//     {
//         feature_vec.push_back(Feature(i, 0.0));
//     }
//     result_idx.clear();
//     result_pair.clear();
//     while (!feature_vec.empty() && sample_size_pair.current_ <= num_element)
//     {
//         std::chrono::high_resolution_clock::time_point start, end; 
//         std::chrono::duration<double> duration;
//         start = std::chrono::high_resolution_clock::now();

//         sample_vec.clear();
//         sample_vec.resize(sample_size_pair.current_);
//         for (int i = 0; i < sample_vec.size(); i++)
//         {
//             sample_vec[i] = rand() % num_element;
//         }

//         end = std::chrono::high_resolution_clock::now();
//         duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//         sort_time += duration.count();

//         for (int i = 0; i < feature_vec.size(); ++i)
//         {
//             int feature_idx = feature_vec[i].idx_;
//             VarianceEst variance(sample_vec, mode, sample_size_pair, cumulate_pair[feature_idx], data_mat[feature_idx], pf0, precise_time);
//             feature_vec[i].Update(variance.low_, variance.up_, variance.sample_variance_);
//         }

//         start = std::chrono::high_resolution_clock::now();

//         count += sample_size_pair.current_ * feature_vec.size();
//         std::vector<Feature>::const_iterator iter = feature_vec.begin();
//         printf("SampleSize %d ExistedFeature %d\n", sample_size_pair.current_, int(feature_vec.size()));

//         // if (sample_size_pair.current_ >= 2000000)
//         // {
//         //     for (int i = 0; i < feature_vec.size(); ++i)
//         //     {
//         //         printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f\n", feature_vec[i].idx_, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_);
//         //     }
//         // }

//         while (iter != feature_vec.end())
//         {
//             if (iter->diff_ < 2.0 * error_bound * threshold)
//             {
//                 if (iter->avg_ >= threshold)
//                 {
//                     result_idx.push_back(iter->idx_);
//                     result_pair.push_back(std::pair<int, double>(iter->idx_, iter->avg_));
//                     printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_);
//                 }
//                 // else if (iter->idx_ == 70)
//                 // {
//                 //     printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f lower_th %f\n", iter->idx_, iter->value_low_, iter->value_up_, iter->value_avg_, iter->value_diff_, (1.0 - error_bound) * threshold);
//                 // }
//                 iter = feature_vec.erase(iter);
//             }
//             else if (iter->low_ >= (1.0 - error_bound) * threshold)
//             {
//                 result_idx.push_back(iter->idx_);
//                 result_pair.push_back(std::pair<int, double>(iter->idx_, iter->avg_));
//                 printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f lower_th %f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_, (1.0 - error_bound) * threshold);
//                 iter = feature_vec.erase(iter);
//             }
//             else if (iter->up_ < (1.0 + error_bound) * threshold)
//             {
//                 iter = feature_vec.erase(iter);
//             }

//             // if (iter->low_ >= (1.0 - error_bound) * threshold)
//             // {
//             //     result_idx.push_back(iter->idx_);
//             //     result_pair.push_back(std::pair<int, double>(iter->idx_, iter->avg_));
//             //     printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f lower_th %f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_, (1.0 - error_bound) * threshold);
//             //     iter = feature_vec.erase(iter);
//             // }
//             // else if (iter->up_ < (1.0 + error_bound) * threshold)
//             // {
//             //     iter = feature_vec.erase(iter);
//             // }
//             // else if (iter->diff_ < 2.0 * error_bound * threshold)
//             // {
//             //     if (iter->avg_ >= threshold)
//             //     {
//             //         result_idx.push_back(iter->idx_);
//             //         result_pair.push_back(std::pair<int, double>(iter->idx_, iter->avg_));
//             //         printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_);
//             //     }
//             //     // else if (iter->idx_ == 70)
//             //     // {
//             //     //     printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f lower_th %f\n", iter->idx_, iter->value_low_, iter->value_up_, iter->value_avg_, iter->value_diff_, (1.0 - error_bound) * threshold);
//             //     // }
//             //     iter = feature_vec.erase(iter);
//             // }

//             else
//             {
//                 ++iter;
//             }
//         }
//         if (sample_size_pair.current_ >= num_element / 2 && sample_size_pair.current_ < num_element)
//         {
//             sample_size_pair.previous_ = sample_size_pair.current_;
//             sample_size_pair.current_ = num_element; 
//             for (int i = 0; i < feature_vec.size(); i++)
//             {
//                 int feature_idx = feature_vec[i].idx_;
//                 double exact_variance = ExactVarianceScore(data_mat[feature_idx]);
//                 feature_vec[i].Update(exact_variance, exact_variance, exact_variance);
//                 if (feature_vec[i].low_>= threshold)
//                 {
//                     result_idx.push_back(feature_idx);
//                     result_pair.push_back(std::pair<int, double>(feature_idx, feature_vec[i].avg_));
//                     printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f\n", feature_idx, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_);
//                 }
//             } 
//             printf("ApproximateVarianceFilter has been found when elements are used up.\n");
//             break;
//         }
//         else
//         {
//             sample_size_pair.previous_ = sample_size_pair.current_;
//             sample_size_pair.current_ *= 2; 
//         }

//         end = std::chrono::high_resolution_clock::now();
//         duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//         sort_time += duration.count();

//     }

//     printf("ApproximateVarianceFilter has been found.\n");
//     printf("SampleSize %d Count %ld PreciseTime %f SortTime %f\n", sample_size_pair.previous_, count, precise_time, sort_time);
    
//     // for (int i = 0; i < result_idx.size(); ++i)
//     // {
//     //     printf("%d ", result_idx[i]);
//     // }
//     // printf("\n");
// }

// void VarianceFilterMoMSimple(const double& ratio_num_block, const double& threshold, const int& initial_sample_size, const std::vector<std::vector<double>>& data_mat, std::vector<int>& result_idx, const double& pf, const double& error_bound, std::vector<std::pair<int, double>>& result_pair)
// {
//     const int num_feature = data_mat.size(), num_element = data_mat[0].size();
//     long count = 0;
//     int block_size = initial_sample_size;
//     double pf0 = pf * 1.0 / num_feature / ceil(log2(num_element * 1.0 / initial_sample_size)), precise_time = 0.0, median_time = 0.0, sort_time = 0.0; //
//     const int num_var = ceil(ratio_num_block * log(2.0 / pf0));

//     printf("Element %d Feature %d pf0 %1.15f 1_pf0 %f NumVar %d BlockSize %d\n", num_element, num_feature, pf0, log(1.0 / pf0), num_var, block_size);
    
//     std::vector<Feature> feature_vec;
//     std::vector<int> sample_vec(num_var, 0);
//     std::vector<std::pair<std::vector<double>, std::vector<double>>> MoM_pair(num_feature, std::pair<std::vector<double>, std::vector<double>>(std::vector<double>(num_var), std::vector<double>(num_var)));
//     for (int i = 0; i < num_feature; i++)
//     {
//         feature_vec.push_back(Feature(i, 0.0));
//     }
//     result_idx.clear();
//     result_pair.clear();
//     while (!feature_vec.empty() && block_size * num_var <= num_element)
//     {
//         // bool last_batch = false;
//         // const int num_block = num_element / block_size;

//         // for (int i = 0; i < sample_vec.size(); i++)
//         // {
//         //     const int random_pos = rand() % num_element;
//         //     if (random_pos >= block_size + num_element % block_size)
//         //     {
//         //         sample_vec[i] = rand() % (num_block - 1);
//         //     }
//         //     else
//         //     {
//         //         sample_vec[i] = num_block - 1;
//         //         last_batch = true;
//         //         printf("RandPos %d BlockSize %d Res %d Total %d\n", random_pos, block_size, int(num_element % block_size), int(block_size + num_element % block_size));
//         //     }
//         // }
//         for (int i = 0; i < sample_vec.size(); i++)
//         {
//             sample_vec[i] = rand() % num_element;
//         }
//         std::sort(sample_vec.begin(), sample_vec.end());

//         for (int i = 0; i < feature_vec.size(); ++i)
//         {
//             int feature_idx = feature_vec[i].idx_;
//             VarianceEst variance(sample_vec, block_size, MoM_pair[feature_idx], data_mat[feature_idx], pf0, precise_time, median_time);
//             feature_vec[i].Update(variance.low_, variance.up_, variance.sample_variance_);
//         }
//         count += block_size * num_var * feature_vec.size();
//         std::vector<Feature>::const_iterator iter = feature_vec.begin();
//         printf("SampleSize %d ExistedFeature %d\n", block_size * num_var, int(feature_vec.size()));
//         while (iter != feature_vec.end())
//         {
//             if (iter->diff_ < 2.0 * error_bound * threshold)
//             {
//                 if (iter->avg_ >= threshold)
//                 {
//                     result_idx.push_back(iter->idx_);
//                     result_pair.push_back(std::pair<int, double>(iter->idx_, iter->avg_));
//                     printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_);
//                 }
//                 // else if (iter->idx_ == 33)
//                 // {
//                 //     printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f lower_th %f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_, (1.0 - error_bound) * threshold);
//                 // }
//                 iter = feature_vec.erase(iter);
//             }
//             else if (iter->low_ >= (1.0 - error_bound) * threshold)
//             {
//                 result_idx.push_back(iter->idx_);
//                 result_pair.push_back(std::pair<int, double>(iter->idx_, iter->avg_));
//                 printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f lower_th %f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_, (1.0 - error_bound) * threshold);
//                 iter = feature_vec.erase(iter);
//             }
//             else if (iter->up_ < (1.0 + error_bound) * threshold)
//             {
//                 // if (iter->idx_ == 33)
//                 // {
//                 //     printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f lower_th %f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_, (1.0 + error_bound) * threshold);
//                 // }
//                 iter = feature_vec.erase(iter);
//             }

//             // if (iter->low_ >= (1.0 - error_bound) * threshold)
//             // {
//             //     result_idx.push_back(iter->idx_);
//             //     result_pair.push_back(std::pair<int, double>(iter->idx_, iter->avg_));
//             //     printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f lower_th %f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_, (1.0 - error_bound) * threshold);
//             //     iter = feature_vec.erase(iter);
//             // }
//             // else if (iter->up_ < (1.0 + error_bound) * threshold)
//             // {
//             //     // if (iter->idx_ == 33)
//             //     // {
//             //     //     printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f lower_th %f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_, (1.0 + error_bound) * threshold);
//             //     // }
//             //     iter = feature_vec.erase(iter);
//             // }
//             // else if (iter->diff_ < 2.0 * error_bound * threshold)
//             // {
//             //     if (iter->avg_ >= threshold)
//             //     {
//             //         result_idx.push_back(iter->idx_);
//             //         result_pair.push_back(std::pair<int, double>(iter->idx_, iter->avg_));
//             //         printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_);
//             //     }
//             //     // else if (iter->idx_ == 33)
//             //     // {
//             //     //     printf("idx %3d low %2.6f up %1.6f avg %1.6f diff %1.10f lower_th %f\n", iter->idx_, iter->low_, iter->up_, iter->avg_, iter->diff_, (1.0 - error_bound) * threshold);
//             //     // }
//             //     iter = feature_vec.erase(iter);
//             // }

//             else
//             {
//                 ++iter;
//             }
//         }
//         if (block_size * num_var >= num_element / 2)
//         {
//             std::chrono::high_resolution_clock::time_point start, end; 
//             std::chrono::duration<double> duration;
//             start = std::chrono::high_resolution_clock::now();
//             printf("SampleSize %d ExistedFeature %d\n", num_element, int(feature_vec.size()));
//             for (int i = 0; i < feature_vec.size(); i++)
//             {
//                 int feature_idx = feature_vec[i].idx_;
//                 double exact_variance = ExactVarianceScore(data_mat[feature_idx]);
//                 feature_vec[i].Update(exact_variance, exact_variance, exact_variance);
//                 if (feature_vec[i].low_>= threshold)
//                 {
//                     result_idx.push_back(feature_idx);
//                     result_pair.push_back(std::pair<int, double>(feature_idx, feature_vec[i].avg_));
//                     printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f\n", feature_idx, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_);
//                 }
//                 printf("idx %3d low %1.6f up %1.6f avg %1.6f diff %1.10f threshold %1.10f\n", feature_idx, feature_vec[i].low_, feature_vec[i].up_, feature_vec[i].avg_, feature_vec[i].diff_, threshold);
//             } 
//             end = std::chrono::high_resolution_clock::now();
//             duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//             sort_time += duration.count();

//             printf("ApproximateVarianceFilter has been found when elements are used up.\n");
//             printf("SampleSize %d Count %ld PreciseTime %f SortTime %f\n", num_element, count, precise_time, sort_time);
//             break;
//         }
//         else
//         {
//             block_size *= 2;
//         }
//     }
//     printf("ApproximateVarianceFilter has been found.\n");
    
//     // for (int i = 0; i < result_idx.size(); ++i)
//     // {
//     //     printf("%d ", result_idx[i]);
//     // }
//     // printf("\n");
// }

void CompareVariance(const std::vector<std::vector<double>>& data_mat, const int topk, const int initial_sample_size)
{
    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();

    int num_feature = data_mat.size(), num_element = data_mat[0].size(), max_sample_size = 1048576, num_point = log2(1048576 / 256) + 1, num_iteration = num_element / max_sample_size;
    std::vector<double> exact_feature_vec(num_feature), exact_sort_vec, exact_topk_vec, point_value(num_point);
    

    std::chrono::high_resolution_clock::time_point cal_start, cal_end;
    std::chrono::duration<double> cal_duration;
    cal_start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < exact_feature_vec.size(); i++)
    {
        double exact_variance_score = ExactVarianceScore(data_mat[i]);
        exact_feature_vec[i] = exact_variance_score;
    }
    exact_sort_vec.assign(exact_feature_vec.begin(), exact_feature_vec.end());
    std::sort(exact_sort_vec.begin(), exact_sort_vec.end(), std::greater<double>());
    exact_topk_vec.assign(exact_sort_vec.begin(), exact_sort_vec.begin() + topk);

    cal_end = std::chrono::high_resolution_clock::now();
	cal_duration = std::chrono::duration_cast<std::chrono::duration<double>>(cal_end - cal_start);
    printf("CalDuration %f\n\n", cal_duration.count());

    for (int k = 0; k < num_iteration; k++)
    {
        SampleSizePair sample_size_pair(0, initial_sample_size);
        std::vector<double> approx_feature_vec(num_feature);
        std::vector<std::pair<double, double>> cumulate_pair(num_feature, std::pair<double, double>(0.0, 0.0));
        // while (sample_size_pair.current_ <= max_sample_size)
        for (int l = 0; l < point_value.size(); l++)
        {
            for (int i = 0; i < num_feature; i++)
            {
                double variance_score = 0.0;
                for (int j = k * max_sample_size + sample_size_pair.previous_; j < k * max_sample_size + sample_size_pair.current_; j++)
                {
                    cumulate_pair[i].first += data_mat[i][j];
                    cumulate_pair[i].second += pow(data_mat[i][j], 2.0);
                }
                variance_score = (cumulate_pair[i].second / sample_size_pair.current_ - pow(cumulate_pair[i].first / sample_size_pair.current_, 2.0)) * sample_size_pair.current_ / (sample_size_pair.current_ - 1.0);
                approx_feature_vec[i] = variance_score;
                // printf("idx %d exact %1.10f approx %1.10f\n", i, exact_feature_vec[i], approx_feature_vec[i]);
            }
            std::vector<double> approx_sort_vec, approx_topk_vec;
            approx_sort_vec.assign(approx_feature_vec.begin(), approx_feature_vec.end());
            std::sort(approx_sort_vec.begin(), approx_sort_vec.end(), std::greater<double>());
            approx_topk_vec.assign(approx_sort_vec.begin(), approx_sort_vec.begin() + topk);

            RelativeError total_relative_error =  TopkRelativeError(exact_feature_vec, approx_feature_vec), topk_relative_error =  TopkRelativeError(exact_topk_vec, approx_topk_vec);

            printf("SampleSize %d Total MRE %f ARE %f Top%d MRE %f ARE %f\n", sample_size_pair.current_, total_relative_error.max_, total_relative_error.avg_, topk, topk_relative_error.max_, topk_relative_error.avg_);
            point_value[l] += total_relative_error.avg_;
            if (sample_size_pair.current_ > num_element / 2 && sample_size_pair.current_ < num_element)
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ = num_element;
            }
            else
            {
                sample_size_pair.previous_ = sample_size_pair.current_;
                sample_size_pair.current_ *= 2;
            }
        }
        printf("\n");
    }
    printf("Average Statistics\n");
    for (int i = 0; i < point_value.size(); i++)
    {
        printf("SampleSize %d ARE %f\n", int(initial_sample_size * pow(2, i)), point_value[i] / num_iteration);
    }
}

void CompareAverage(const std::vector<std::vector<double>>& data_mat, const int topk, const int initial_sample_size)
{
    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();

    int num_feature = data_mat.size(), num_element = data_mat[0].size();
    SampleSizePair sample_size_pair(0, initial_sample_size);
    std::vector<double> exact_feature_vec(num_feature), approx_feature_vec(num_feature), exact_sort_vec, exact_topk_vec;
    std::vector<double> cumulate_sum(num_feature,  0.0);

    std::chrono::high_resolution_clock::time_point cal_start, cal_end;
    std::chrono::duration<double> cal_duration;
    cal_start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < exact_feature_vec.size(); i++)
    {
        double exact_average_score = ExactAverageScore(data_mat[i]);
        exact_feature_vec[i] = exact_average_score;
    }
    exact_sort_vec.assign(exact_feature_vec.begin(), exact_feature_vec.end());
    std::sort(exact_sort_vec.begin(), exact_sort_vec.end(), std::greater<double>());
    exact_topk_vec.assign(exact_sort_vec.begin(), exact_sort_vec.begin() + topk);

    cal_end = std::chrono::high_resolution_clock::now();
	cal_duration = std::chrono::duration_cast<std::chrono::duration<double>>(cal_end - cal_start);
    printf("CalDuration %f\n\n", cal_duration.count());

    while (sample_size_pair.current_ <= num_element)
    {
        for (int i = 0; i < num_feature; i++)
        {
            double average_score = 0.0;
            for (int j = sample_size_pair.previous_; j < sample_size_pair.current_; j++)
            {
                cumulate_sum[i] += data_mat[i][j];
            }
            average_score = cumulate_sum[i] / sample_size_pair.current_;
            approx_feature_vec[i] = average_score;
            // printf("idx %d exact %1.10f approx %1.10f\n", i, exact_feature_vec[i], approx_feature_vec[i]);
        }
        std::vector<double> approx_sort_vec, approx_topk_vec;
        approx_sort_vec.assign(approx_feature_vec.begin(), approx_feature_vec.end());
        std::sort(approx_sort_vec.begin(), approx_sort_vec.end(), std::greater<double>());
        approx_topk_vec.assign(approx_sort_vec.begin(), approx_sort_vec.begin() + topk);

        RelativeError total_relative_error =  TopkRelativeError(exact_feature_vec, approx_feature_vec), topk_relative_error =  TopkRelativeError(exact_topk_vec, approx_topk_vec);

        printf("SampleSize %d Total MRE %f ARE %f Top%d MRE %f ARE %f\n", sample_size_pair.current_, total_relative_error.max_, total_relative_error.avg_, topk, topk_relative_error.max_, topk_relative_error.avg_);
        if (sample_size_pair.current_ > num_element / 2 && sample_size_pair.current_ < num_element)
        {
            sample_size_pair.previous_ = sample_size_pair.current_;
            sample_size_pair.current_ = num_element;
        }
        else
        {
            sample_size_pair.previous_ = sample_size_pair.current_;
            sample_size_pair.current_ *= 2;
        }
    }
    printf("\n");
}

void ExactVarianceTopk(const int top_k, const std::vector<std::vector<double>>& data_mat, std::vector<int>& result_idx, std::vector<int>& result_idx_plus, std::vector<double>& variance_score, std::vector<double>& result_value)
{
	std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();

    int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    std::vector<Feature> feature_vec;

    std::chrono::high_resolution_clock::time_point cal_start, cal_end;
    std::chrono::duration<double> cal_duration;
    cal_start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_feature; i++)
    {
        double exact_variance_score = ExactVarianceScore(data_mat[i]);
        feature_vec.push_back(Feature(i, exact_variance_score));
    }

    cal_end = std::chrono::high_resolution_clock::now();
	cal_duration = std::chrono::duration_cast<std::chrono::duration<double>>(cal_end - cal_start);
    printf("CalDuration %f\n\n", cal_duration.count());

    result_idx.clear();
    result_idx_plus.clear();
    result_value.clear();
    std::vector<Feature> temp(feature_vec);
    sort(temp.begin(), temp.end(), cmp_avg);
    for (int i = 0; i < temp.size(); i++)
    {
        if (i < top_k)
        {
            result_idx.push_back(temp[i].idx_);
            result_idx_plus.push_back(temp[i].idx_);  
            result_value.push_back(temp[i].avg_);
        }
        else if (abs(temp[i].avg_ - temp[top_k - 1].avg_) <= 1e-12) // With difference smaller than 1e-12, we regard them as equal values.  
        {
            result_idx_plus.push_back(temp[i].idx_);
        }
        if (i < 20)
        {
            printf("idx %3d avg %1.6f\n", temp[i].idx_, temp[i].avg_);
        }
    }
    variance_score.resize(num_feature);
    for (int i = 0; i < feature_vec.size(); i++)
    {
        variance_score[i] = feature_vec[i].avg_;
        // printf("idx %3d avg %1.6f\n", feature_vec[i].idx_, feature_vec[i].avg_);
    }

    end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    printf("Feature %d Element %d Count %ld ExactDuration %f\n\n", num_feature, num_element, long(num_feature) * long(num_element), duration.count());
}

// void ExactVarianceTopkPred(const int top_k, const PredPair& pred, const std::vector<std::vector<double>>& data_mat, std::vector<int>& result_idx, std::vector<int>& result_idx_plus, std::vector<double>& variance_score, std::vector<double>& result_value)
// {
// 	std::chrono::high_resolution_clock::time_point start, end;
//     std::chrono::duration<double> duration;
//     start = std::chrono::high_resolution_clock::now();

//     int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
//     std::vector<Feature> feature_vec;

//     std::chrono::high_resolution_clock::time_point cal_start, cal_end;
//     std::chrono::duration<double> cal_duration;
//     cal_start = std::chrono::high_resolution_clock::now();

//     for (int i = 0; i < num_feature; i++)
//     {
//         if (i != pred.col_)
//         {
//             double exact_variance_score = ExactPredVarianceScore(pred.value_, data_mat[pred.col_], data_mat[i]);
//             feature_vec.push_back(Feature(i, exact_variance_score));
//         }
//     }

//     cal_end = std::chrono::high_resolution_clock::now();
// 	cal_duration = std::chrono::duration_cast<std::chrono::duration<double>>(cal_end - cal_start);
//     printf("CalDuration %f\n\n", cal_duration.count());

//     result_idx.clear();
//     result_idx_plus.clear();
//     result_value.clear();
//     std::vector<Feature> temp(feature_vec);
//     sort(temp.begin(), temp.end(), cmp_avg);
//     for (int i = 0; i < temp.size(); i++)
//     {
//         if (i < top_k)
//         {
//             result_idx.push_back(temp[i].idx_);
//             result_idx_plus.push_back(temp[i].idx_);  
//             result_value.push_back(temp[i].avg_);
//         }
//         else if (abs(temp[i].avg_ - temp[top_k - 1].avg_) <= 1e-12) // With difference smaller than 1e-12, we regard them as equal values.  
//         {
//             result_idx_plus.push_back(temp[i].idx_);
//         }
//         if (i < 20)
//         {
//             printf("idx %3d avg %1.6f\n", temp[i].idx_, temp[i].avg_);
//         }
//     }
//     variance_score.resize(num_feature);
//     for (int i = 0; i < feature_vec.size(); i++)
//     {
//         variance_score[i] = feature_vec[i].avg_;
//         // printf("idx %3d avg %1.6f\n", feature_vec[i].idx_, feature_vec[i].avg_);
//     }

//     end = std::chrono::high_resolution_clock::now();
// 	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//     printf("Feature %d Element %d Count %ld ExactDuration %f\n\n", num_feature, num_element, long(num_feature) * long(num_element), duration.count());
// }

void ExactVarianceTopkPred(const int top_k, const PredPair& pred, const std::vector<std::vector<double>>& data_mat, std::vector<int>& result_idx, std::vector<int>& result_idx_plus, std::vector<double>& variance_score, std::vector<double>& result_value)
{
	std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double> duration;
    start = std::chrono::high_resolution_clock::now();

    int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    std::vector<Feature> feature_vec;

    std::chrono::high_resolution_clock::time_point cal_start, cal_end;
    std::chrono::duration<double> cal_duration;
    cal_start = std::chrono::high_resolution_clock::now();

    std::vector<int> sample_vec;
    for (int i = 0; i < num_element; i++)
    {
        if (abs(data_mat[pred.col_][i] - pred.value_) <= 1e-14)
        {
            sample_vec.push_back(i);
        }
    }
    for (int i = 0; i < num_feature; i++)
    {
        if (i != pred.col_)
        {
            double sum = 0.0, sum_square = 0.0, variance_score = 0.0;
            for (int j = 0; j < sample_vec.size(); j++)
            {
                int sample_index = sample_vec[j];
                sum += data_mat[i][sample_index];
                sum_square += pow(data_mat[i][sample_index], 2.0);
            }
            variance_score = sum_square / sample_vec.size() - pow(sum / sample_vec.size(), 2.0);
            feature_vec.push_back(Feature(i, variance_score));
        }
    }

    cal_end = std::chrono::high_resolution_clock::now();
	cal_duration = std::chrono::duration_cast<std::chrono::duration<double>>(cal_end - cal_start);
    printf("CalDuration %f\n\n", cal_duration.count());

    result_idx.clear();
    result_idx_plus.clear();
    result_value.clear();
    std::vector<Feature> temp(feature_vec);
    sort(temp.begin(), temp.end(), cmp_avg);
    for (int i = 0; i < temp.size(); i++)
    {
        if (i < top_k)
        {
            result_idx.push_back(temp[i].idx_);
            result_idx_plus.push_back(temp[i].idx_);  
            result_value.push_back(temp[i].avg_);
        }
        else if (abs(temp[i].avg_ - temp[top_k - 1].avg_) <= 1e-12) // With difference smaller than 1e-12, we regard them as equal values.  
        {
            result_idx_plus.push_back(temp[i].idx_);
        }
        if (i < 20)
        {
            printf("idx %3d avg %1.6f\n", temp[i].idx_, temp[i].avg_);
        }
    }
    variance_score.resize(num_feature);
    for (int i = 0; i < feature_vec.size(); i++)
    {
        variance_score[i] = feature_vec[i].avg_;
        // printf("idx %3d avg %1.6f\n", feature_vec[i].idx_, feature_vec[i].avg_);
    }

    end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    printf("Feature %d Element %d Count %ld ExactDuration %f\n\n", num_feature, num_element, long(num_feature) * long(num_element), duration.count());
}

void ExactVarianceFilter(const double& threshold, const std::vector<std::vector<double>>& data_mat, std::vector<int>& result_idx, std::vector<std::pair<int, double>>& result_pair)
{
    result_idx.clear();
    result_pair.clear();
    int num_feature = data_mat.size(), num_element = data_mat[0].size();
    std::vector<Feature> feature_vec;
    for (int i = 0; i < num_feature; i++)
    {
        double exact_variance_score = ExactVarianceScore(data_mat[i]);
        Feature feature(i, exact_variance_score);
        feature_vec.push_back(feature);
        if (exact_variance_score >= threshold)
        {
            result_idx.push_back(i);
            result_pair.push_back(std::pair<int, double>(i, feature.avg_));
            //printf("idx %3d low %1.6f up %1.6f avg %1.6f\n", i, attribute.value_low_, attribute.value_up_, attribute.value_avg_);
        }
        
    }
    printf("Feature %d Element %d Count %ld\n", num_feature, num_element, long(num_feature) * long(num_element));
    // for (int i = 0; i < result_idx.size(); ++i)
    // {
    //     printf("%d ", result_idx[i]);
    // }
    // printf("\n");
}