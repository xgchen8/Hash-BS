// variance_est.h
#ifndef VARIANCE_EST_H
#define VARIANCE_EST_H
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <chrono>
#include "feature.h"
#include "sample_size_pair.h"

#define likely(x) __builtin_expect(!!(x), 1) 
#define unlikely(x) __builtin_expect(!!(x), 0)
class VarianceEst
{
private:
    void UpdateCumRandomIndex(const std::vector<uint32_t>& sample_vec, const SampleSizePair& sample_size_pair, std::pair<double, double>& cum_pair, const std::vector<double>& data_vec, double& precise)
    {
        std::chrono::high_resolution_clock::time_point start, end; 
        std::chrono::duration<double> duration;
        start = std::chrono::high_resolution_clock::now();
        
        cum_pair = std::pair<double, double>(0.0, 0.0);
        for (int i = 0; i < sample_size_pair.current_; i++)
        {
            int sample_index = sample_vec[i];
            cum_pair.first += data_vec[sample_index];
            cum_pair.second += pow(data_vec[sample_index], 2.0);
        }
        end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		precise += duration.count();
    }

    void UpdateCumMomIndex(const int& num_var, const std::vector<uint32_t>& sample_vec, const int& block_size, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const std::vector<double>& data_vec, double& precise)
    {
        std::chrono::high_resolution_clock::time_point start, end; 
        std::chrono::duration<double> duration;
        start = std::chrono::high_resolution_clock::now();
        
        for (int i = 0; i < num_var; i++)
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
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        precise += duration.count();
    }

    void UpdateBound(Feature& value_est, const double& cum_value, const double& fail_pr)
    {
        double a = log(1.0 / fail_pr);
        // printf("low %f up %f\n", pow(sqrt(sample_count + 2.0 / 9.0 * a) - sqrt(a / 2.0), 2.0) - a / 18.0, pow(sqrt(sample_count + a / 2.0) + sqrt(a / 2.0), 2.0));
        value_est.Update((pow(sqrt(cum_value + 2.0 / 9.0 * a) - sqrt(a / 2.0), 2.0) - a / 18.0) * 1.0 / num_sample_, pow(sqrt(cum_value + a / 2.0) + sqrt(a / 2.0), 2.0) * 1.0 / num_sample_, cum_value / num_sample_);
    }

    void UpdateBoundColt(const std::pair<double, double>& cum_pair, const double& fail_pr)
    {
        double a = log(1.0 / fail_pr), sample_stdev = sqrt((cum_pair.second / num_sample_ - pow(cum_pair.first / num_sample_, 2.0)) * num_sample_ / (num_sample_ - 1.0)), bias = sqrt(2 * a / (num_sample_ - 1.0));
        low_ = pow(std::max(sample_stdev - bias, 0.0), 2.0);
        up_ = pow(sample_stdev + bias, 2.0);
        sample_variance_ = pow(sample_stdev, 2.0);
        // printf("sample_variance_ %f\n", sample_variance_);
    }

    double UpdateBoundMoMWR4Wise(std::vector<double>& MoM_vec, const double& fail_pr, const int& block_size)
    {
        int median = MoM_vec.size() / 2;
        double a = log(1.0 / fail_pr);
        std::nth_element(MoM_vec.begin(), MoM_vec.begin() + median, MoM_vec.end());
        double est_mean = MoM_vec[median], 
        p_epsilon = 0.5 - sqrt(a / 2.0 / MoM_vec.size()), 
        indicator = sqrt(3.0 / block_size / block_size / p_epsilon + 4.0 * sqrt(3.0 / p_epsilon) / block_size * est_mean * (1.0 - est_mean)+ 2.0 / block_size / block_size / sqrt(3.0 * p_epsilon) + 2.0 / p_epsilon / block_size / block_size / block_size), 
        denominator = 2.0 * (1.0 + sqrt(3.0 / p_epsilon) / block_size), 
        base = 2.0 * est_mean + sqrt(3.0 / p_epsilon) / block_size,
        // printf("Bias %f Tempa %f Tempr %d L %f U %f max_diff %f\n", bias, a, int(MoM_vec.size()), min_mean, max_mean, max_diff);
        min_mean = std::max((base - indicator) / denominator, 0.0),
        max_mean = std::min((base + indicator) / denominator, 1.0), 
        max_diff = std::max(max_mean - 0.0, 1.0 - min_mean);
        // printf("Tempa %f Tempr %d L %f U %f max_diff %f\n", a, int(MoM_vec.size()), min_mean, max_mean, max_diff);
        // printf("Tempa %f Tempr %d estMean %f L %f U %f Base %f Denominator %f Indicator %f\n", a, int(MoM_vec.size()), est_mean, low_, up_, base, denominator, indicator);
        // max_diff_square = std::max(up_ - 0.0, 1.0 - low_) * std::max(up_ - 0.0, 1.0 - low_);
        return max_diff * max_diff;
    }

    void UpdateBoundMoVWRTest(std::vector<double>& MoV_vec, const double& fail_pr, const int& block_size, const double& max_diff_square)
    {
        int median = MoV_vec.size() / 2;
        double a = log(1.0 / fail_pr);
        std::nth_element(MoV_vec.begin(), MoV_vec.begin() + median, MoV_vec.end());
        sample_variance_ = MoV_vec[median]; 
        double constant = 0.5 - sqrt(a / 2.0 / MoV_vec.size()),
        indicator = sqrt(max_diff_square * max_diff_square + 4.0 * block_size * constant * sample_variance_ * (max_diff_square - (block_size - 3.0) / (block_size - 1.0) * sample_variance_)),
        denominator = 2.0 * (block_size * constant + (block_size - 3.0) / (block_size - 1.0)), 
        base = 2.0 * block_size * constant * sample_variance_ + max_diff_square;
        low_ = (base - indicator) / denominator;
        up_ = (base + indicator) / denominator;
        static bool onceMoM = true;
        if (onceMoM)
        {
            // printf("Hoeffding Constant %f Tempa %f Tempk %d L %f U %f\n", constant, a, int(MoV.size()), (base - indicator) / denominator, (base + indicator) / denominator);
            printf("Hoeffding Constant %f Tempa %f Tempk %d L %f U %f\n", constant, a, int(MoV_vec.size()), low_, up_);
        }
        onceMoM = false;
        // printf("\n");
    }
    
	int num_sample_, num_element_;
	double fail_pr_;
    std::string mode_;
    Feature variable_, variable_square_;
public:
	VarianceEst(): num_sample_(0), num_element_(0), fail_pr_(0.25), low_(0.0), up_(0.0), sample_variance_(0.0)
    {}
	VarianceEst(const std::string& mode, const SampleSizePair& sample_size_pair, std::pair<double, double>& cum_pair, const std::vector<double>& data_vec, const double& fail_pr, double& precise): num_sample_(sample_size_pair.current_), num_element_(data_vec.size()), fail_pr_(fail_pr), mode_(mode), low_(0.0), up_(0.0)
    {
        if (mode_ == "random")
        {
            double single_fail_pr = fail_pr_ / 4.0;
            UpdateBound(variable_, cum_pair.first, single_fail_pr);
            UpdateBound(variable_square_, cum_pair.second, single_fail_pr);
            low_ = std::max(variable_square_.low_ - pow(variable_.up_, 2.0), 0.0);
            up_ = std::min(variable_square_.up_ - pow(std::max(variable_.low_, 0.0), 2.0), 1.0);
            sample_variance_ = (low_ + up_) / 2.0;
        }
        else if (mode_ == "colt")
        {
            double single_fail_pr = fail_pr_ / 2.0;
            UpdateBoundColt(cum_pair, single_fail_pr);
        }
        else
        {
            printf("The mode is incorrect!\n");
        }
    }
    VarianceEst(const std::vector<uint32_t>& sample_vec, const std::string& mode, const SampleSizePair& sample_size_pair, std::pair<double, double>& cum_pair, const std::vector<double>& data_vec, const double& fail_pr, double& precise): num_sample_(sample_size_pair.current_), num_element_(data_vec.size()), fail_pr_(fail_pr), mode_(mode), low_(0.0), up_(0.0)
    {
        UpdateCumRandomIndex(sample_vec, sample_size_pair, cum_pair, data_vec, precise);
        if (mode_ == "random")
        {
            double single_fail_pr = fail_pr_ / 4.0;
            UpdateBound(variable_, cum_pair.first, single_fail_pr);
            UpdateBound(variable_square_, cum_pair.second, single_fail_pr);
            low_ = std::max(variable_square_.low_ - pow(variable_.up_, 2.0), 0.0);
            up_ = std::min(variable_square_.up_ - pow(std::max(variable_.low_, 0.0), 2.0), 1.0);
            sample_variance_ = (low_ + up_) / 2.0;
        }
        else if (mode_ == "colt")
        {
            double single_fail_pr = fail_pr_ / 2.0;
            UpdateBoundColt(cum_pair, single_fail_pr);
        }
        else
        {
            printf("The mode is incorrect!\n");
        }
    }
    VarianceEst(const int& block_size, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const double& fail_pr): fail_pr_(fail_pr), low_(0.0), up_(0.0)
    {
        double single_fail_pr = fail_pr_ / 2.0;
        // UpdateMoVTestAnyPosition(sample_vec, block_size, MoM_pair, data_vec, precise, median_time);
        double max_diff_square = UpdateBoundMoMWR4Wise(MoM_pair.first, single_fail_pr, block_size);
        UpdateBoundMoVWRTest(MoM_pair.second, single_fail_pr, block_size, max_diff_square);

        // UpdateMoV(sample_vec, block_size, MoM_pair.first, data_vec, precise, median_time, last_batch);
        // UpdateBoundMoVWR(variable_, MoM_pair.first, single_fail_pr, block_size);
        
        // UpdateBoundMoMPairwise(variable_square_, MoM_pair.second, single_fail_pr, block_size);
    }
    VarianceEst(const int& num_var, const int& block_size, std::vector<uint32_t>& sample_vec, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const std::vector<double>& data_vec, const double& fail_pr, double& precise): num_element_(data_vec.size()), fail_pr_(fail_pr), low_(0.0), up_(0.0)
    {
        double single_fail_pr = fail_pr_ / 2.0;
        // UpdateMoVTestAnyPosition(sample_vec, block_size, MoM_pair, data_vec, precise, median_time);
        UpdateCumMomIndex(num_var, sample_vec, block_size, MoM_pair, data_vec, precise);
        double max_diff_square = UpdateBoundMoMWR4Wise(MoM_pair.first, single_fail_pr, block_size);
        UpdateBoundMoVWRTest(MoM_pair.second, single_fail_pr, block_size, max_diff_square);

        // UpdateMoV(sample_vec, block_size, MoM_pair.first, data_vec, precise, median_time, last_batch);
        // UpdateBoundMoVWR(variable_, MoM_pair.first, single_fail_pr, block_size);
        
        // UpdateBoundMoMPairwise(variable_square_, MoM_pair.second, single_fail_pr, block_size);
    }
	~VarianceEst()
    {}
    double low_, up_, sample_variance_;

};
#endif

    // void UpdateMoVTest(const std::vector<int>& sample_vec, const int& block_size, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const std::vector<double>& data_vec, double& precise, double& median_time, const bool& last_batch)
    // {
    //     std::chrono::high_resolution_clock::time_point start, end; 
    //     std::chrono::duration<double> duration;
    //     start = std::chrono::high_resolution_clock::now();
    //     if (likely(last_batch == false)) 
    //     {
    //         for (int i = 0; i < sample_vec.size(); i++)
    //         {
    //             double sum = 0.0, sum_square = 0.0;
    //             for (int j = sample_vec[i] * block_size; j < (sample_vec[i] + 1) * block_size; j++)
    //             {
    //                 sum += data_vec[j];
    //                 sum_square += pow(data_vec[j], 2.0);/* code */
    //             }
    //             MoM_pair.first[i] = sum / block_size;
    //             MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
    //         }
    //     }
    //     else
    //     {
    //         for (int i = 0; i < sample_vec.size(); i++)
    //         {
    //             int res = 0;
    //             double sum = 0.0, sum_square = 0.0;
    //             for (int j = sample_vec[i] * block_size; j < (sample_vec[i] + 1) * block_size; j++)
    //             {
    //                 sum += data_vec[j];
    //                 sum_square += pow(data_vec[j], 2.0);/* code */
    //             }        
    //             if (sample_vec[i] == num_element_ / block_size - 1) // last block
    //             {
    //                 printf("StartIdx %d BlockSize %d\n", sample_vec[i], block_size);
    //                 for (int j = (sample_vec[i] + 1) * block_size; j < num_element_; j++)
    //                 {
    //                     sum += data_vec[j];
    //                     sum_square += pow(data_vec[j], 2.0);
    //                 }
    //                 res = num_element_ % block_size;
    //             }
    //             MoM_pair.first[i] = sum / (block_size + res);
    //             MoM_pair.second[i] = (sum_square / (block_size + res) - pow(sum / (block_size + res), 2.0)) * (block_size + res) / (block_size + res - 1.0);
    //         }
    //     }
    //     end = std::chrono::high_resolution_clock::now();
    //     duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    //     precise += duration.count();
    // }

    // void UpdateMoVTestAnyPosition(const std::vector<int>& sample_vec, const int& block_size, std::pair<std::vector<double>, std::vector<double>>& MoM_pair, const std::vector<double>& data_vec, double& precise, double& median_time)
    // {
    //     std::chrono::high_resolution_clock::time_point start, end; 
    //     std::chrono::duration<double> duration;
    //     start = std::chrono::high_resolution_clock::now();
    //     for (int i = 0; i < sample_vec.size(); i++)
    //     {
    //         double sum = 0.0, sum_square = 0.0;
    //         if (likely(sample_vec[i] <= num_element_ - block_size))
    //         {
    //             for (int j = sample_vec[i]; j < sample_vec[i] + block_size; j++)
    //             {
    //                 sum += data_vec[j];
    //                 sum_square += pow(data_vec[j], 2.0);/* code */
    //             }
    //         }
    //         else
    //         {   
    //             int res = sample_vec[i] + block_size - num_element_;
    //             printf("initialPos %d blockSize %d numElement %d res %d\n", sample_vec[i], block_size, num_element_, res);
    //             for (int j = sample_vec[i]; j < num_element_; j++)
    //             {
    //                 sum += data_vec[j];
    //                 sum_square += pow(data_vec[j], 2.0);/* code */
    //             }
    //             for (int j = 0; j < res; j++)
    //             {
    //                 sum += data_vec[j];
    //                 sum_square += pow(data_vec[j], 2.0);/* code */
    //             }
    //         }     
    //         MoM_pair.first[i] = sum / block_size;
    //         MoM_pair.second[i] = (sum_square / block_size - pow(sum / block_size, 2.0)) * block_size / (block_size - 1.0);
    //     }
    //     end = std::chrono::high_resolution_clock::now();
    //     duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    //     precise += duration.count();
    // }

    
    // void UpdateBoundMoVWR(Feature& value_est, std::vector<double>& MoV, const double& fail_pr, const int& block_size)
    // {
    //     int median = MoV.size() / 2;
    //     double a = log(1.0 / fail_pr);
    //     std::nth_element(MoV.begin(), MoV.begin() + median, MoV.end());
    //     double est_mean = MoV[median], 
    //     bias = 1.0 / sqrt(block_size * (0.5 + 2.0 * a / 3.0 / MoV.size() - sqrt(4.0 * a * a / 9.0 / MoV.size() / MoV.size() +  a / MoV.size()))), 
    //     indicator = sqrt(bias * bias * bias * bias + 4.0 * est_mean * bias * bias - 4 * (block_size - 3.0) / (block_size - 1.0) * bias * bias * est_mean * est_mean), denominator = 2.0 * (1.0 + bias * bias * (block_size - 3.0) / (block_size - 1.0)), 
    //     base = 2.0 * est_mean + bias * bias;
    //     value_est.Update((base - indicator) / denominator, (base + indicator) / denominator);
        
    //     static bool onceMoM = true;
    //     if (onceMoM)
    //     {
    //         printf("Bias %f Tempa %f Tempk %d L %f U %f\n", bias, a, int(MoV.size()), 
    //         (base - indicator) / denominator, (base + indicator) / denominator);
    //     }
    //     onceMoM = false;
    //     // printf("\n");
    // }

    // double UpdateBoundMoMWR(std::vector<double>& MoM, const double& fail_pr, const int& block_size)
    // {
    //     int median = MoM.size() / 2;
    //     double a = log(1.0 / fail_pr);
    //     std::nth_element(MoM.begin(), MoM.begin() + median, MoM.end());
    //     double est_mean = MoM[median], 
    //     bias = 1.0 / sqrt(block_size * (0.5 + 2.0 * a / 3.0 / MoM.size() - sqrt(4.0 * a * a / 9.0 / MoM.size() / MoM.size() +  a / MoM.size()))), 
    //     indicator = sqrt(bias * bias * bias * bias + 4.0 * est_mean * bias * bias), 
    //     denominator = 2.0, 
    //     base = 2.0 * est_mean + bias * bias,
    //     min_mean = std::max((base - indicator) / denominator, 0.0), 
    //     max_mean = std::min((base + indicator) / denominator, 1.0),
    //     max_diff = std::max(max_mean - 0.0, 1.0 - min_mean);
    //     printf("Bias %f Tempa %f Tempr %d L %f U %f max_diff %f\n", bias, a, int(MoM.size()), min_mean, max_mean, max_diff);
    //     //     (base - indicator) / denominator, (base + indicator) / denominator);
    //     // static bool onceMoM = true;
    //     // if (onceMoM)
    //     // {
    //     //     
    //     // }
    //     // onceMoM = false;
    //     // printf("\n");
    //     return max_diff * max_diff;
    // }