#ifndef QUERY_H
#define QUERY_H
#include "functions.h"
#include "config.h"
#include "measure.h"
#include "file_ctrl.h"
#include <random>
#include <string.h>
void VarianceTopkQuery(const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const Config& config)
{
	int num_feature = data_mat.size() - 1, num_element = data_mat[0].size(), initial_sample_size = config.initial_size_;
    printf("NumofFeature %d NumofRecord %d\n", num_feature, num_element);
	double pf = 1.0 / double(num_element), exact_time = 0.0, random_time = 0.0, colt_time = 0.0, mom_time = 0.0, mom3_time = 0.0, mom4_time = 0.0;
	RelativeError random_relative_error(0.0, 0.0), colt_relative_error(0.0, 0.0), mom_relative_error(0.0, 0.0), mom3_relative_error(0.0, 0.0), mom4_relative_error(0.0, 0.0), temp_relative_error;  
	Measure random_measure(0.0, 0.0, 0.0, 0.0, 0.0), colt_measure(0.0, 0.0, 0.0, 0.0, 0.0), mom_measure(0.0, 0.0, 0.0, 0.0, 0.0), mom3_measure(0.0, 0.0, 0.0, 0.0, 0.0), mom4_measure(0.0, 0.0, 0.0, 0.0, 0.0), temp_measure;
	std::vector<int> sample_index_count(num_element);

	printf("Initial Sample Size %d\nTop_%d\n", initial_sample_size, config.topk_);
	std::chrono::high_resolution_clock::time_point start, end;
	std::chrono::duration<double> duration;
    for (int t = 0; t < config.query_num_; t++)
    {
		// int rand_cnt = 0;
		// std::vector<uint8_t> rand_seq(num_element);
		// for (int i = 0; i < rand_seq.size(); i++)
		// {
		// 	rand_seq[i] = rand();
		// 	if (i % 10000000 == 0)
		// 	{
		// 		printf("i = %d: %d\n", i, rand_seq[i]);
		// 	}
		// }
		std::vector<double> exact_variance_score, exact_value, approx_value;
		std::vector<int> exact_variance_topk, exact_variance_topk_plus, approx_variance_topk;
		start = std::chrono::high_resolution_clock::now();
		ExactVarianceTopk(config.topk_, data_mat, exact_variance_topk, exact_variance_topk_plus, exact_variance_score, exact_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		printf("ExactVarianceTopkTime: %f\n", duration.count());
		exact_time += duration.count();

		start = std::chrono::high_resolution_clock::now();
		ApproxVarianceTopk("random", config.topk_, initial_sample_size, data_mat, approx_variance_topk, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		temp_measure = Measure(Precision(exact_variance_topk_plus, approx_variance_topk), Recall(exact_variance_topk_plus, approx_variance_topk), F1Measure(exact_variance_topk_plus, approx_variance_topk), NDCG(exact_variance_score, exact_variance_topk, approx_variance_topk), 0.0);
		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
		printf("VarianceTopkRandomTime: %f\n", duration.count());
		random_time += duration.count();
		temp_relative_error = TopkRelativeError(exact_value, approx_value);
        random_relative_error += temp_relative_error;
		random_measure += temp_measure; 

		start = std::chrono::high_resolution_clock::now();
		ApproxVarianceTopk("colt", config.topk_, initial_sample_size, data_mat, approx_variance_topk, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		temp_measure = Measure(Precision(exact_variance_topk_plus, approx_variance_topk), Recall(exact_variance_topk_plus, approx_variance_topk), F1Measure(exact_variance_topk_plus, approx_variance_topk), NDCG(exact_variance_score, exact_variance_topk, approx_variance_topk), 0.0);
		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
		printf("VarianceTopkColtTime: %f\n", duration.count());
		colt_time += duration.count();
		temp_relative_error = TopkRelativeError(exact_value, approx_value);
        colt_relative_error += temp_relative_error;
		colt_measure += temp_measure; 

		// start = std::chrono::high_resolution_clock::now();
		// VarianceTopkMoMSimple(config.ratio_num_block_, config.topk_, initial_sample_size, data_mat, cum_size, sample_index_count, approx_variance_topk, pf, config.epsilon_, approx_value);
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// temp_measure = Measure(Precision(exact_variance_topk_plus, approx_variance_topk), Recall(exact_variance_topk_plus, approx_variance_topk), F1Measure(exact_variance_topk_plus, approx_variance_topk), NDCG(exact_variance_score, exact_variance_topk, approx_variance_topk), 0.0);
		// printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
		// printf("ApproxVarianceTopkMoMUpdate3Time: %f\n", duration.count());
		// mom_time += duration.count();
		// temp_relative_error = TopkRelativeError(exact_value, approx_value);
        // mom_relative_error += temp_relative_error;
		// mom_measure += temp_measure; 

		start = std::chrono::high_resolution_clock::now();
		VarianceTopkMomBlock3(config.ratio_num_block_, config.topk_, initial_sample_size, data_mat, cum_size, approx_variance_topk, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		temp_measure = Measure(Precision(exact_variance_topk_plus, approx_variance_topk), Recall(exact_variance_topk_plus, approx_variance_topk), F1Measure(exact_variance_topk_plus, approx_variance_topk), NDCG(exact_variance_score, exact_variance_topk, approx_variance_topk), 0.0);
		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
		printf("ApproxVarianceTopkMonNonSimpleTime: %f\n", duration.count());
		mom3_time += duration.count();
		temp_relative_error = TopkRelativeError(exact_value, approx_value);
        mom3_relative_error += temp_relative_error;
		mom3_measure += temp_measure; 

		// start = std::chrono::high_resolution_clock::now();
		// VarianceTopkMomBlock4(config.ratio_num_block_, config.topk_, initial_sample_size, data_mat, cum_size, approx_variance_topk, pf, config.epsilon_, approx_value);
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// temp_measure = Measure(Precision(exact_variance_topk_plus, approx_variance_topk), Recall(exact_variance_topk_plus, approx_variance_topk), F1Measure(exact_variance_topk_plus, approx_variance_topk), NDCG(exact_variance_score, exact_variance_topk, approx_variance_topk), 0.0);
		// printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
		// printf("ApproxVarianceTopkMonNonSimpleTime: %f\n", duration.count());
		// mom4_time += duration.count();
		// temp_relative_error = TopkRelativeError(exact_value, approx_value);
        // mom4_relative_error += temp_relative_error;
		// mom4_measure += temp_measure; 
	}
	exact_time /= config.query_num_;
	printf("Exact %f\n", exact_time); 

	random_time /= config.query_num_;
	random_relative_error /= config.query_num_;
	random_measure /= config.query_num_;
	printf("Random %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", random_time, random_relative_error.max_, random_relative_error.avg_, random_measure.precision_, random_measure.recall_, random_measure.f1_, random_measure.NDCG_);

	colt_time /= config.query_num_;
	colt_relative_error /= config.query_num_;
	colt_measure /= config.query_num_;
	printf("Colt   %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", colt_time, colt_relative_error.max_, colt_relative_error.avg_, colt_measure.precision_, colt_measure.recall_, colt_measure.f1_, colt_measure.NDCG_);

	// mom_time /= config.query_num_;
	// mom_relative_error /= config.query_num_;
	// mom_measure /= config.query_num_;
	// printf("MoM    %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", mom_time, mom_relative_error.max_, mom_relative_error.avg_, mom_measure.precision_, mom_measure.recall_, mom_measure.f1_, mom_measure.NDCG_);

	mom3_time /= config.query_num_;
	mom3_relative_error /= config.query_num_;
	mom3_measure /= config.query_num_;
	printf("MomBlock3   %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", mom3_time, mom3_relative_error.max_, mom3_relative_error.avg_, mom3_measure.precision_, mom3_measure.recall_, mom3_measure.f1_, mom3_measure.NDCG_);

	// mom4_time /= config.query_num_;
	// mom4_relative_error /= config.query_num_;
	// mom4_measure /= config.query_num_;
	// printf("MomBlock4  %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", mom4_time, mom4_relative_error.max_, mom4_relative_error.avg_, mom4_measure.precision_, mom4_measure.recall_, mom4_measure.f1_, mom4_measure.NDCG_);
}

// void PredVarianceTopkQuery(const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const Config& config)
// {
// 	int num_feature = data_mat.size() - 1, num_element = data_mat[0].size(), initial_sample_size = config.initial_size_;
//     printf("NumofFeature %d NumofRecord %d\n", num_feature, num_element);
// 	double pf = 1.0 / double(num_element), exact_time = 0.0, random_time = 0.0, colt_time = 0.0, mom_time = 0.0, mom3_time = 0.0;
// 	RelativeError random_relative_error(0.0, 0.0), colt_relative_error(0.0, 0.0), mom_relative_error(0.0, 0.0), mom3_relative_error(0.0, 0.0), temp_relative_error;  
// 	Measure random_measure(0.0, 0.0, 0.0, 0.0, 0.0), colt_measure(0.0, 0.0, 0.0, 0.0, 0.0), mom_measure(0.0, 0.0, 0.0, 0.0, 0.0), mom3_measure(0.0, 0.0, 0.0, 0.0, 0.0), temp_measure;

// 	printf("Initial Sample Size %d\nTop_%d\n", initial_sample_size, config.topk_);
// 	std::chrono::high_resolution_clock::time_point start, end;
// 	std::chrono::duration<double> duration;
//     for (int t = 0; t < config.query_num_; t++)
//     {
// 		std::vector<double> exact_variance_score, exact_value, approx_value;
// 		std::vector<int> exact_variance_topk, exact_variance_topk_plus, approx_variance_topk;
// 		start = std::chrono::high_resolution_clock::now();
// 		ExactVarianceTopk(config.topk_, data_mat, exact_variance_topk, exact_variance_topk_plus, exact_variance_score, exact_value);
// 		end = std::chrono::high_resolution_clock::now();
// 		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
// 		printf("ExactVarianceTopkTime: %f\n", duration.count());
// 		exact_time += duration.count();

// 		start = std::chrono::high_resolution_clock::now();
// 		ApproxVarianceTopkPred("random", config.topk_, config.pred_, initial_sample_size, data_mat, approx_variance_topk, pf, config.epsilon_, approx_value);
// 		end = std::chrono::high_resolution_clock::now();
// 		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
// 		temp_measure = Measure(Precision(exact_variance_topk_plus, approx_variance_topk), Recall(exact_variance_topk_plus, approx_variance_topk), F1Measure(exact_variance_topk_plus, approx_variance_topk), NDCG(exact_variance_score, exact_variance_topk, approx_variance_topk), 0.0);
// 		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
// 		printf("PredVarianceTopkRandomTime: %f\n", duration.count());
// 		random_time += duration.count();
// 		temp_relative_error = TopkRelativeError(exact_value, approx_value);
//         random_relative_error += temp_relative_error;
// 		random_measure += temp_measure; 

// 		start = std::chrono::high_resolution_clock::now();
// 		ApproxVarianceTopkPred("colt", config.topk_, config.pred_, initial_sample_size, data_mat, approx_variance_topk, pf, config.epsilon_, approx_value);
// 		end = std::chrono::high_resolution_clock::now();
// 		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
// 		temp_measure = Measure(Precision(exact_variance_topk_plus, approx_variance_topk), Recall(exact_variance_topk_plus, approx_variance_topk), F1Measure(exact_variance_topk_plus, approx_variance_topk), NDCG(exact_variance_score, exact_variance_topk, approx_variance_topk), 0.0);
// 		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
// 		printf("PredVarianceTopkColtTime: %f\n", duration.count());
// 		colt_time += duration.count();
// 		temp_relative_error = TopkRelativeError(exact_value, approx_value);
//         colt_relative_error += temp_relative_error;
// 		colt_measure += temp_measure; 

// 		start = std::chrono::high_resolution_clock::now();
// 		MomVarianceTopkPred(config.ratio_num_block_, config.topk_, config.pred_, initial_sample_size, data_mat, cum_size, approx_variance_topk, pf, config.epsilon_, approx_value);
// 		end = std::chrono::high_resolution_clock::now();
// 		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
// 		temp_measure = Measure(Precision(exact_variance_topk_plus, approx_variance_topk), Recall(exact_variance_topk_plus, approx_variance_topk), F1Measure(exact_variance_topk_plus, approx_variance_topk), NDCG(exact_variance_score, exact_variance_topk, approx_variance_topk), 0.0);
// 		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
// 		printf("PredApproxVarianceTopkMonNonSimpleTime: %f\n", duration.count());
// 		mom3_time += duration.count();
// 		temp_relative_error = TopkRelativeError(exact_value, approx_value);
//         mom3_relative_error += temp_relative_error;
// 		mom3_measure += temp_measure; 
// 	}
// 	exact_time /= config.query_num_;
// 	printf("Exact %f\n", exact_time); 

// 	random_time /= config.query_num_;
// 	random_relative_error /= config.query_num_;
// 	random_measure /= config.query_num_;
// 	printf("Random %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", random_time, random_relative_error.max_, random_relative_error.avg_, random_measure.precision_, random_measure.recall_, random_measure.f1_, random_measure.NDCG_);

// 	colt_time /= config.query_num_;
// 	colt_relative_error /= config.query_num_;
// 	colt_measure /= config.query_num_;
// 	printf("Colt   %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", colt_time, colt_relative_error.max_, colt_relative_error.avg_, colt_measure.precision_, colt_measure.recall_, colt_measure.f1_, colt_measure.NDCG_);

// 	// mom_time /= config.query_num_;
// 	// mom_relative_error /= config.query_num_;
// 	// mom_measure /= config.query_num_;
// 	// printf("MoM    %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", mom_time, mom_relative_error.max_, mom_relative_error.avg_, mom_measure.precision_, mom_measure.recall_, mom_measure.f1_, mom_measure.NDCG_);

// 	mom3_time /= config.query_num_;
// 	mom3_relative_error /= config.query_num_;
// 	mom3_measure /= config.query_num_;
// 	printf("MomBlock3   %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", mom3_time, mom3_relative_error.max_, mom3_relative_error.avg_, mom3_measure.precision_, mom3_measure.recall_, mom3_measure.f1_, mom3_measure.NDCG_);
// }

void SelectAbsPredTopkQuery(const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const Config& config)
{
	int num_feature = data_mat.size() - 1, num_element = data_mat[0].size(), initial_sample_size = config.initial_size_;
    printf("NumofFeature %d NumofRecord %d\n", num_feature, num_element);
	double pf = 1.0 / double(num_element), exact_time = 0.0, random_time = 0.0, colt_time = 0.0, mom_time = 0.0, mom3_time = 0.0;
	AbsoluteError random_abs_error(0.0, 0.0), colt_abs_error(0.0, 0.0), mom_abs_error(0.0, 0.0), mom3_abs_error(0.0, 0.0), temp_abs_error;  
	Measure random_measure(0.0, 0.0, 0.0, 0.0, 0.0), colt_measure(0.0, 0.0, 0.0, 0.0, 0.0), mom_measure(0.0, 0.0, 0.0, 0.0, 0.0), mom3_measure(0.0, 0.0, 0.0, 0.0, 0.0), temp_measure;

	printf("Initial Sample Size %d\nTop_%d\n", initial_sample_size, config.topk_);
	std::chrono::high_resolution_clock::time_point start, end;
	std::chrono::duration<double> duration;
    for (int t = 0; t < config.query_num_; t++)
    {
		std::vector<double> exact_variance_score, exact_value, approx_value;
		std::vector<int> exact_variance_topk, exact_variance_topk_plus, approx_variance_topk;
		start = std::chrono::high_resolution_clock::now();
		ExactVarianceTopkPred(config.topk_, config.pred_, data_mat, exact_variance_topk, exact_variance_topk_plus, exact_variance_score, exact_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		printf("ExactVarianceTopkTime: %f\n", duration.count());
		exact_time += duration.count();

		start = std::chrono::high_resolution_clock::now();
		SelectAbsPredTopkVariance("random", config.topk_, config.pred_, initial_sample_size, data_mat, approx_variance_topk, pf, config.abs_error_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		temp_measure = Measure(Precision(exact_variance_topk_plus, approx_variance_topk), Recall(exact_variance_topk_plus, approx_variance_topk), F1Measure(exact_variance_topk_plus, approx_variance_topk), NDCG(exact_variance_score, exact_variance_topk, approx_variance_topk), 0.0);
		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
		printf("PredVarianceTopkRandomTime: %f\n", duration.count());
		random_time += duration.count();
		temp_abs_error = TopkAbsoluteError(exact_value, approx_value);
        random_abs_error += temp_abs_error;
		random_measure += temp_measure; 

		start = std::chrono::high_resolution_clock::now();
		SelectAbsPredTopkVariance("colt", config.topk_, config.pred_, initial_sample_size, data_mat, approx_variance_topk, pf, config.abs_error_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		temp_measure = Measure(Precision(exact_variance_topk_plus, approx_variance_topk), Recall(exact_variance_topk_plus, approx_variance_topk), F1Measure(exact_variance_topk_plus, approx_variance_topk), NDCG(exact_variance_score, exact_variance_topk, approx_variance_topk), 0.0);
		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
		printf("PredVarianceTopkColtTime: %f\n", duration.count());
		colt_time += duration.count();
		temp_abs_error = TopkAbsoluteError(exact_value, approx_value);
        colt_abs_error += temp_abs_error;
		colt_measure += temp_measure; 

		start = std::chrono::high_resolution_clock::now();
		SelectAbsPredMomTopkVariance(config.ratio_num_block_, config.topk_, config.pred_, initial_sample_size, data_mat, cum_size, approx_variance_topk, pf, config.abs_error_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		temp_measure = Measure(Precision(exact_variance_topk_plus, approx_variance_topk), Recall(exact_variance_topk_plus, approx_variance_topk), F1Measure(exact_variance_topk_plus, approx_variance_topk), NDCG(exact_variance_score, exact_variance_topk, approx_variance_topk), 0.0);
		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.NDCG_);
		printf("PredApproxVarianceTopkMonNonSimpleTime: %f\n", duration.count());
		mom3_time += duration.count();
		temp_abs_error = TopkAbsoluteError(exact_value, approx_value);
        mom3_abs_error += temp_abs_error;
		mom3_measure += temp_measure; 
	}
	exact_time /= config.query_num_;
	printf("Exact %f\n", exact_time); 

	random_time /= config.query_num_;
	random_abs_error /= config.query_num_;
	random_measure /= config.query_num_;
	printf("Random %f MAE %f AAE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", random_time, random_abs_error.max_, random_abs_error.avg_, random_measure.precision_, random_measure.recall_, random_measure.f1_, random_measure.NDCG_);

	colt_time /= config.query_num_;
	colt_abs_error /= config.query_num_;
	colt_measure /= config.query_num_;
	printf("Colt   %f MAE %f AAE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", colt_time, colt_abs_error.max_, colt_abs_error.avg_, colt_measure.precision_, colt_measure.recall_, colt_measure.f1_, colt_measure.NDCG_);

	// mom_time /= config.query_num_;
	// mom_abs_error /= config.query_num_;
	// mom_measure /= config.query_num_;
	// printf("MoM    %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", mom_time, mom_abs_error.max_, mom_abs_error.avg_, mom_measure.precision_, mom_measure.recall_, mom_measure.f1_, mom_measure.NDCG_);

	mom3_time /= config.query_num_;
	mom3_abs_error /= config.query_num_;
	mom3_measure /= config.query_num_;
	printf("Mom    %f MAE %f AAE %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", mom3_time, mom3_abs_error.max_, mom3_abs_error.avg_, mom3_measure.precision_, mom3_measure.recall_, mom3_measure.f1_, mom3_measure.NDCG_);
}

void EpsilonMetric(const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const Config& config)
{
	int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    printf("NumofFeature %d NumofRecord %d\n", num_feature, num_element);
	double pf = 1.0 / double(num_element), exact_time = 0.0, random_time = 0.0, colt_time = 0.0, mom_time = 0.0, mom2_time = 0.0, bernstein_time = 0.0;
	RelativeError random_relative_error(0.0, 0.0), colt_relative_error(0.0, 0.0), mom_relative_error(0.0, 0.0), mom2_relative_error(0.0, 0.0), bernstein_relative_error(0.0, 0.0), temp_relative_error;  
	printf("Initial Sample Size %d\nTop_%d\n", config.initial_size_, config.topk_);
	std::chrono::high_resolution_clock::time_point start, end;
	std::chrono::duration<double> duration;
	// std::vector<uint32_t> sample_index_count(num_element);
	// std::vector<int> test_feature = {72};
	std::vector<int> test_feature(num_feature);
	//std::vector<int> test_feature(config.num_test_feature_); // = {6, 2, 106, 3, 18, 114, 58, 59, 36, 37}; //{6, 2, 106, 3, 18, 114, 58, 59, 36, 37};
	printf("Idx ");
	// for (int i = 0; i < test_feature.size(); i++)
	// {
	// 	test_feature[i] = rand() % num_feature;
	// 	printf(" %d", test_feature[i]);
	// }

	for (int i = 0; i < test_feature.size(); i++)
	{
		test_feature[i] = i;
		printf(" %d", test_feature[i]);
	}
	printf("\n");
	
    for (int t = 0; t < config.query_num_; t++)
    {
		std::vector<double> exact_value(test_feature.size()), approx_value(test_feature.size());
		start = std::chrono::high_resolution_clock::now();
		// if (config.metric_ == "variance")
		// {	
		for (int i = 0; i < test_feature.size(); i++)
		{
			exact_value[i] = ExactVarianceScore(data_mat[test_feature[i]]);
			printf("idx %3d variance %1.12f\n", i, exact_value[i]);
		}
		// }
		// else if (config.metric_ == "average")
		// {
		// 	for (int i = 0; i < test_feature.size(); i++)
		// 	{
		// 		exact_value[i] = ExactAverageScore(data_mat[test_feature[i]]);
		// 		printf("idx %3d average %f\n", i, exact_value[i]);
		// 	}
		// }
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		printf("ExactVarianceTopkTime: %f\n", duration.count());
		exact_time += duration.count();

		start = std::chrono::high_resolution_clock::now();
		EpsilonDeltaApproxVariance("random", config.empirical_threshold_, test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value);
		// ApproxMetric(config.metric_, "random", test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		random_time += duration.count();
		temp_relative_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
		random_relative_error += temp_relative_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("EpsilonMetricRandomTime: %f MRE %f ARE %f\n\n", duration.count(), temp_relative_error.max_, temp_relative_error.avg_);
        
		// if (config.metric_ == "variance")
		// {
		start = std::chrono::high_resolution_clock::now();
		EpsilonDeltaApproxVariance("colt", config.empirical_threshold_, test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value);
		// ApproxMetric(config.metric_, "colt", test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		colt_time += duration.count();
		temp_relative_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
		colt_relative_error += temp_relative_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("EpsilonVarianceColtTime: %f MRE %f ARE %f\n\n", duration.count(), temp_relative_error.max_, temp_relative_error.avg_);
		// }

		start = std::chrono::high_resolution_clock::now();
		EpsilonDeltaMomVariance(config.ratio_num_block_, config.empirical_threshold_, test_feature, config.initial_size_, data_mat, cum_size, pf, config.epsilon_, approx_value);
		// MoMSimpleMetric(config.ratio_num_block_, config.metric_, test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		mom_time += duration.count();
		temp_relative_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
        mom_relative_error += temp_relative_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("EpsilonMetricMoMUpdate3Time: %f MRE %f ARE %f\n\n", duration.count(), temp_relative_error.max_, temp_relative_error.avg_);

		start = std::chrono::high_resolution_clock::now();
		EpsilonDeltaMomVarianceBlock2(config.ratio_num_block_, config.empirical_threshold_, test_feature, config.initial_size_, data_mat, cum_size, pf, config.epsilon_, approx_value);
		// MoMSimpleMetric(config.ratio_num_block_, config.metric_, test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		mom2_time += duration.count();
		temp_relative_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
        mom2_relative_error += temp_relative_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("EpsilonMetricMomBlock2Time: %f MRE %f ARE %f\n\n", duration.count(), temp_relative_error.max_, temp_relative_error.avg_);

		// if (config.metric_ == "average")
		// {
		// 	start = std::chrono::high_resolution_clock::now();
		// 	BernsteinMetric(config.metric_, test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value, config.block_size_);
		// 	end = std::chrono::high_resolution_clock::now();
		// 	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// 	bernstein_time += duration.count();
		// 	temp_relative_error = TopkRelativeError(exact_value, approx_value);
		// 	bernstein_relative_error += temp_relative_error;
		// 	printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		// 	printf("EpsilonAverageBernsteinTime: %f MRE %f ARE %f\n\n", duration.count(), temp_relative_error.max_, temp_relative_error.avg_);
		// }
		printf("\n");
	}
	exact_time /= (config.query_num_ * test_feature.size());
	printf("Exact %f\n", exact_time); 

	random_time /= (config.query_num_ * test_feature.size());
	random_relative_error /= config.query_num_;
	printf("Random %f MRE %f ARE %f\n", random_time, random_relative_error.max_, random_relative_error.avg_);

	// if (config.metric_ == "variance")
	// {
	colt_time /= (config.query_num_ * test_feature.size());
	colt_relative_error /= config.query_num_;
	printf("Colt %f MRE %f ARE %f\n", colt_time, colt_relative_error.max_, colt_relative_error.avg_);
	// }

	mom_time /= (config.query_num_ * test_feature.size());
	mom_relative_error /= config.query_num_;
	printf("MoM %f MRE %f ARE %f\n", mom_time, mom_relative_error.max_, mom_relative_error.avg_);

	mom2_time /= (config.query_num_ * test_feature.size());
	mom2_relative_error /= config.query_num_;
	printf("Mom2 %f MRE %f ARE %f\n", mom2_time, mom2_relative_error.max_, mom2_relative_error.avg_);

	// if (config.metric_ == "average")
	// {
	// 	bernstein_time /= (config.query_num_ * test_feature.size());
	// 	bernstein_relative_error /= config.query_num_;
	// 	printf("Bernstein %f MRE %f ARE %f\n", bernstein_time, bernstein_relative_error.max_, bernstein_relative_error.avg_);
	// }
}

void EpsilonDeltaPredMetric(const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const Config& config)
{
	int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    printf("NumofFeature %d NumofRecord %d\n", num_feature, num_element);
	double pf = 1.0 / double(num_element), exact_time = 0.0, rand_time = 0.0, colt_time = 0.0, colt2_time = 0.0, mom_time = 0.0, mom2_time = 0.0, mom3_time = 0.0;
	RelativeError rand_rel_error(0.0, 0.0), colt_rel_error(0.0, 0.0), colt2_rel_error(0.0, 0.0), mom_rel_error(0.0, 0.0), mom2_rel_error(0.0, 0.0), mom3_rel_error(0.0, 0.0), temp_rel_error;  
	printf("Initial Sample Size %d\n\n", config.initial_size_);
	std::chrono::high_resolution_clock::time_point start, end;
	std::chrono::duration<double> duration;
	// std::vector<uint32_t> sample_index_count(num_element);
	// std::vector<int> test_feature = {72};
	std::vector<int> test_feature;
	//std::vector<int> test_feature(config.num_test_feature_); // = {6, 2, 106, 3, 18, 114, 58, 59, 36, 37}; //{6, 2, 106, 3, 18, 114, 58, 59, 36, 37};
	printf("Idx ");
	// for (int i = 0; i < test_feature.size(); i++)
	// {
	// 	test_feature[i] = rand() % num_feature;
	// 	printf(" %d", test_feature[i]);
	// }

	for (int i = 0; i < num_feature; i++)
	{
		if (i != config.pred_.col_)
		{
			test_feature.push_back(i);
			printf(" %d", i);
		}
	}
	printf("\n");
	
    for (int t = 0; t < config.query_num_; t++)
    {
		std::vector<double> exact_value(test_feature.size()), approx_value(test_feature.size());

		// start = std::chrono::high_resolution_clock::now();
		// for (int i = 0; i < test_feature.size(); i++)
		// {
		// 	exact_value[i] = TestPredVarianceScore(config.pred_.value_, data_mat[config.pred_.col_], data_mat[test_feature[i]], 333312);
		// 	printf("Test idx %3d variance %1.12f\n", test_feature[i], exact_value[i]);
		// }
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// printf("TestPredVarianceTime: %f\n\n", duration.count());

		// start = std::chrono::high_resolution_clock::now();
		// for (int i = 0; i < test_feature.size(); i++)
		// {
		// 	exact_value[i] = TestPredVarianceScore2(config.pred_.value_, data_mat[config.pred_.col_], data_mat[test_feature[i]], 333312);
		// 	printf("Test2 idx %3d variance %1.12f\n", test_feature[i], exact_value[i]);
		// }
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// printf("Test2PredVarianceTime: %f\n\n", duration.count());

		// start = std::chrono::high_resolution_clock::now();
		// for (int i = 0; i < test_feature.size(); i++)
		// {
		// 	exact_value[i] = TestPredVarianceScore3(config.pred_.value_, data_mat[config.pred_.col_], data_mat[test_feature[i]], 333312);
		// 	printf("Test2 idx %3d variance %1.12f\n", test_feature[i], exact_value[i]);
		// }
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// printf("Test3PredVarianceTime: %f\n\n", duration.count());

		start = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < test_feature.size(); i++)
		{
			exact_value[i] = ExactPredVarianceScore(config.pred_.value_, data_mat[config.pred_.col_], data_mat[test_feature[i]]);
			printf("idx %3d variance %1.12f\n", test_feature[i], exact_value[i]);
		}
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		printf("ExactPredVarianceTime: %f\n\n", duration.count());
		exact_time += duration.count();

		start = std::chrono::high_resolution_clock::now();
		EpsilonDeltaPredApproxVariance("random", config.pred_, config.empirical_threshold_, test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value);
		// ApproxMetric(config.metric_, "random", test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		rand_time += duration.count();
		temp_rel_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
		rand_rel_error += temp_rel_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("EpsilonPredMetricRandomTime: %f MRE %f ARE %f\n\n", duration.count(), temp_rel_error.max_, temp_rel_error.avg_);

		start = std::chrono::high_resolution_clock::now();
		EpsilonDeltaPredApproxVariance("colt", config.pred_, config.empirical_threshold_, test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		colt_time += duration.count();
		temp_rel_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
		colt_rel_error += temp_rel_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("EpsilonPredVarianceColtTime: %f MRE %f ARE %f\n\n", duration.count(), temp_rel_error.max_, temp_rel_error.avg_);

		start = std::chrono::high_resolution_clock::now();
		EpsilonDeltaPredApproxVarianceBlock2("colt", config.pred_, config.empirical_threshold_, test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		colt2_time += duration.count();
		temp_rel_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
		colt2_rel_error += temp_rel_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("EpsilonPredVarianceColtBlock2Time: %f MRE %f ARE %f\n\n", duration.count(), temp_rel_error.max_, temp_rel_error.avg_);

		start = std::chrono::high_resolution_clock::now();
		EpsilonDeltaPredMomVariance(config.ratio_num_block_, config.pred_, config.empirical_threshold_, test_feature, config.initial_size_, data_mat, cum_size, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		mom_time += duration.count();
		temp_rel_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
        mom_rel_error += temp_rel_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("EpsilonPredMetricMoMTime: %f MRE %f ARE %f\n\n", duration.count(), temp_rel_error.max_, temp_rel_error.avg_);

		start = std::chrono::high_resolution_clock::now();
		EpsilonDeltaPredMomVarianceBlock2(config.ratio_num_block_, config.pred_, config.empirical_threshold_, test_feature, config.initial_size_, data_mat, cum_size, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		mom2_time += duration.count();
		temp_rel_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
        mom2_rel_error += temp_rel_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("EpsilonPredMetricMoMBlock2Time: %f MRE %f ARE %f\n\n", duration.count(), temp_rel_error.max_, temp_rel_error.avg_);

		start = std::chrono::high_resolution_clock::now();
		EpsilonDeltaPredMomVarianceBlock3(config.ratio_num_block_, config.pred_, config.empirical_threshold_, test_feature, config.initial_size_, data_mat, cum_size, pf, config.epsilon_, approx_value);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		mom3_time += duration.count();
		temp_rel_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
        mom3_rel_error += temp_rel_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("EpsilonPredMetricMoMBlock3Time: %f MRE %f ARE %f\n\n", duration.count(), temp_rel_error.max_, temp_rel_error.avg_);

		printf("\n");
	}
	exact_time /= (config.query_num_ * test_feature.size());
	printf("Exact %f\n", exact_time); 

	rand_time /= (config.query_num_ * test_feature.size());
	rand_rel_error /= config.query_num_;
	printf("Random %f MRE %f ARE %f\n", rand_time, rand_rel_error.max_, rand_rel_error.avg_);

	colt_time /= (config.query_num_ * test_feature.size());
	colt_rel_error /= config.query_num_;
	printf("Colt   %f MRE %f ARE %f\n", colt_time, colt_rel_error.max_, colt_rel_error.avg_);

	colt2_time /= (config.query_num_ * test_feature.size());
	colt2_rel_error /= config.query_num_;
	printf("Colt2  %f MRE %f ARE %f\n", colt2_time, colt2_rel_error.max_, colt2_rel_error.avg_);

	mom_time /= (config.query_num_ * test_feature.size());
	mom_rel_error /= config.query_num_;
	printf("MoM    %f MRE %f ARE %f\n", mom_time, mom_rel_error.max_, mom_rel_error.avg_);

	mom2_time /= (config.query_num_ * test_feature.size());
	mom2_rel_error /= config.query_num_;
	printf("MoM2   %f MRE %f ARE %f\n", mom2_time, mom2_rel_error.max_, mom2_rel_error.avg_);

	mom3_time /= (config.query_num_ * test_feature.size());
	mom3_rel_error /= config.query_num_;
	printf("MoM3   %f MRE %f ARE %f\n", mom3_time, mom3_rel_error.max_, mom3_rel_error.avg_);
}

void SelectAbsPredQuery(const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const Config& config)
{
	int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    printf("NumofFeature %d NumofRecord %d\n", num_feature, num_element);
	double pf = 1.0 / double(num_element), exact_time = 0.0, rand_time = 0.0, colt_time = 0.0, mom_time = 0.0;
	AbsoluteError rand_abs_error(0.0, 0.0), colt_abs_error(0.0, 0.0), mom_abs_error(0.0, 0.0), temp_abs_error;  
	printf("Initial Sample Size %d\n\n", config.initial_size_);
	std::chrono::high_resolution_clock::time_point start, end;
	std::chrono::duration<double> duration;
	std::vector<int> test_feature;
	printf("Idx ");

	for (int i = 0; i < num_feature; i++)
	{
		if (i != config.pred_.col_)
		{
			test_feature.push_back(i);
			printf(" %d", i);
		}
	}
	printf("\n");
	
    for (int t = 0; t < config.query_num_; t++)
    {
		std::vector<double> exact_variance(test_feature.size()), approx_variance(test_feature.size());
		// std::vector<std::pair<double, double>> ;

		start = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < test_feature.size(); i++)
		{
			exact_variance[i] = ExactPredVarianceScore(config.pred_.value_, data_mat[config.pred_.col_], data_mat[test_feature[i]]);
			printf("idx %3d variance %1.12f\n", test_feature[i], exact_variance[i]);
		}
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		printf("ExactPredVarianceTime: %f\n\n", duration.count());
		exact_time += duration.count();

		start = std::chrono::high_resolution_clock::now();
		SelectAbsPredApproxVariance("random", config.pred_, test_feature, config.initial_size_, data_mat, pf, config.abs_error_, approx_variance);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		rand_time += duration.count();
		temp_abs_error = SelectAbsoluteError(exact_variance, approx_variance);
		rand_abs_error += temp_abs_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("SelectAbsPredRandomTime: %f MAE %f AAE %f\n\n", duration.count(), temp_abs_error.max_, temp_abs_error.avg_);

		start = std::chrono::high_resolution_clock::now();
		SelectAbsPredApproxVariance("colt", config.pred_, test_feature, config.initial_size_, data_mat, pf, config.abs_error_, approx_variance);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		colt_time += duration.count();
		temp_abs_error = SelectAbsoluteError(exact_variance, approx_variance);
		colt_abs_error += temp_abs_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("SelectAbsPredColtTime: %f MAE %f AAE %f\n\n", duration.count(), temp_abs_error.max_, temp_abs_error.avg_);

		// start = std::chrono::high_resolution_clock::now();
		// EpsilonDeltaPredApproxVarianceBlock2("colt", config.pred_, config.empirical_threshold_, test_feature, config.initial_size_, data_mat, pf, config.epsilon_, approx_value);
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// colt2_time += duration.count();
		// temp_abs_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
		// colt2_abs_error += temp_abs_error;
		// printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		// printf("EpsilonPredVarianceColtBlock2Time: %f MRE %f ARE %f\n\n", duration.count(), temp_abs_error.max_, temp_abs_error.avg_);

		start = std::chrono::high_resolution_clock::now();
		SelectAbsPredMomVariance(config.ratio_num_block_, config.pred_, test_feature, config.initial_size_, data_mat, cum_size, pf, config.abs_error_, approx_variance);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		mom_time += duration.count();
		temp_abs_error = SelectAbsoluteError(exact_variance, approx_variance);
        mom_abs_error += temp_abs_error;
		printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		printf("SelectAbsPredMoMTime: %f MAE %f AAE %f\n\n", duration.count(), temp_abs_error.max_, temp_abs_error.avg_);

		// start = std::chrono::high_resolution_clock::now();
		// EpsilonDeltaPredMomVarianceBlock2(config.ratio_num_block_, config.pred_, config.empirical_threshold_, test_feature, config.initial_size_, data_mat, cum_size, pf, config.epsilon_, approx_value);
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// mom2_time += duration.count();
		// temp_rel_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
        // mom2_rel_error += temp_rel_error;
		// printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		// printf("EpsilonPredMetricMoMBlock2Time: %f MRE %f ARE %f\n\n", duration.count(), temp_rel_error.max_, temp_rel_error.avg_);

		// start = std::chrono::high_resolution_clock::now();
		// EpsilonDeltaPredMomVarianceBlock3(config.ratio_num_block_, config.pred_, config.empirical_threshold_, test_feature, config.initial_size_, data_mat, cum_size, pf, config.epsilon_, approx_value);
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// mom3_time += duration.count();
		// temp_rel_error = EpsilonDeltaRelativeError(exact_value, approx_value, config.empirical_threshold_);
        // mom3_rel_error += temp_rel_error;
		// printf("Approx%1.4f %f\n", config.epsilon_, duration.count());
		// printf("EpsilonPredMetricMoMBlock3Time: %f MRE %f ARE %f\n\n", duration.count(), temp_rel_error.max_, temp_rel_error.avg_);

		printf("\n");
	}
	exact_time /= (config.query_num_ * test_feature.size());
	printf("Exact %f\n", exact_time); 

	rand_time /= (config.query_num_ * test_feature.size());
	rand_abs_error /= config.query_num_;
	printf("Random %f MAE %f AAE %f\n", rand_time, rand_abs_error.max_, rand_abs_error.avg_);

	colt_time /= (config.query_num_ * test_feature.size());
	colt_abs_error /= config.query_num_;
	printf("Colt   %f MAE %f AAE %f\n", colt_time, colt_abs_error.max_, colt_abs_error.avg_);

	// colt2_time /= (config.query_num_ * test_feature.size());
	// colt2_abs_error /= config.query_num_;
	// printf("Colt2  %f MRE %f ARE %f\n", colt2_time, colt2_abs_error.max_, colt2_abs_error.avg_);

	mom_time /= (config.query_num_ * test_feature.size());
	mom_abs_error /= config.query_num_;
	printf("MoM    %f MAE %f AAE %f\n", mom_time, mom_abs_error.max_, mom_abs_error.avg_);

	// mom2_time /= (config.query_num_ * test_feature.size());
	// mom2_rel_error /= config.query_num_;
	// printf("MoM2   %f MRE %f ARE %f\n", mom2_time, mom2_rel_error.max_, mom2_rel_error.avg_);

	// mom3_time /= (config.query_num_ * test_feature.size());
	// mom3_rel_error /= config.query_num_;
	// printf("MoM3   %f MRE %f ARE %f\n", mom3_time, mom3_rel_error.max_, mom3_rel_error.avg_);
}

void CacheMissStat(const std::vector<std::vector<double>>& data_mat, const std::vector<uint32_t>& cum_size, const Config& config)
{
	int num_feature = data_mat.size() - 1, num_element = data_mat[0].size();
    printf("NumofFeature %d NumofRecord %d\n", num_feature, num_element);
	double pf = 1.0 / double(num_element), exact_time = 0.0, rand_time = 0.0, colt_time = 0.0, mom_time = 0.0;
	printf("Initial Sample Size %d\n\n", config.initial_size_);
	std::vector<int> test_feature;
	printf("Idx ");
	for (int i = 0; i < num_feature; i++)
	{
		if (i != config.pred_.col_)
		{
			test_feature.push_back(i);
			printf(" %d", i);
		}
	}
	printf("%s %s\n", config.action_.c_str(), config.cache_miss_mode_.c_str());
	std::vector<double> exact_variance(test_feature.size()), approx_variance(test_feature.size());
	std::vector<double> exact_variance_score, exact_value, approx_value;
	std::vector<int> exact_variance_topk, exact_variance_topk_plus, approx_variance_topk;
	std::chrono::high_resolution_clock::time_point start, end;
	std::chrono::duration<double> duration;
	start = std::chrono::high_resolution_clock::now();

	if (config.cache_miss_mode_ == "single_exact")
	{
		std::vector<double> exact_variance(test_feature.size());
		for (int i = 0; i < test_feature.size(); i++)
		{
			exact_variance[i] = ExactPredVarianceScore(config.pred_.value_, data_mat[config.pred_.col_], data_mat[test_feature[i]]);
			printf("idx %3d variance %1.12f\n", test_feature[i], exact_variance[i]);
		}
	}
	else if (config.cache_miss_mode_ == "topk_exact")
	{	
		ExactVarianceTopkPred(config.topk_, config.pred_, data_mat, exact_variance_topk, exact_variance_topk_plus, exact_variance_score, exact_value);
	}
	
	else if (config.cache_miss_mode_ == "single_random")
	{
		SelectAbsPredRandomVarianceStat(config.pred_, test_feature, config.initial_size_, data_mat, pf, config.abs_error_, approx_variance);
	}
	else if (config.cache_miss_mode_ == "topk_random")
	{
		SelectAbsPredTopkVariance("random", config.topk_, config.pred_, config.initial_size_, data_mat, approx_variance_topk, pf, config.abs_error_, approx_value);
	}
	else if (config.cache_miss_mode_ == "single_colt")
	{
		SelectAbsPredColtVarianceStat(config.pred_, test_feature, config.initial_size_, data_mat, pf, config.abs_error_, approx_variance);
	}
	else if (config.cache_miss_mode_ == "topk_colt")
	{
		SelectAbsPredTopkVariance("colt", config.topk_, config.pred_, config.initial_size_, data_mat, approx_variance_topk, pf, config.abs_error_, approx_value);
	}
	
    else if (config.cache_miss_mode_ == "single_mom")
	{
		SelectAbsPredMomVarianceStat(config.ratio_num_block_, config.pred_, test_feature, config.initial_size_, data_mat, cum_size, pf, config.abs_error_, approx_variance);
	}
	else if (config.cache_miss_mode_ == "topk_mom")
	{
		SelectAbsPredMomTopkVariance(config.ratio_num_block_, config.topk_, config.pred_, config.initial_size_, data_mat, cum_size, approx_variance_topk, pf, config.abs_error_, approx_value);
	}
	
	end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    printf("Duration %f\n", duration);
}

void HashFile(const std::string in_file_name, const Config& config)
{
    std::vector<std::string> temp_file;
	std::string out_file_name = config.prefix_ + config.dataset_ + "_" + config.hash_model_  + ".csv", first_line, line_str, slot_size_file_name = config.prefix_ + config.dataset_ + "_" + config.hash_model_  + "_cumulate_size.csv";
    std::ifstream in_file(in_file_name, std::ios::in);
	std::ofstream out_file(out_file_name, std::ios::out), slot_size_file(slot_size_file_name, std::ios::out);
    getline(in_file, first_line);
	out_file << first_line << '\n';
    int count = 1;
    while (getline(in_file, line_str))
    {
        //printf("%d\n", count);
        count++;
        //cout << line_str << endl; 
        temp_file.push_back(line_str);
    }
    in_file.close();
    int num_record = temp_file.size();
    std::vector<int> permutation_vec(num_record);
    for (int i = 0; i < permutation_vec.size(); i++)
    {
        permutation_vec[i] = i;
    }	
	if (config.hash_model_ == "tabulation")
	{
		HashChain hash_chain(4, 25);
		for (int i = 0; i < permutation_vec.size(); i++)
		{
			hash_chain.Insert(permutation_vec[i]);
		}
		for (int i = 0; i < hash_chain.table_.size(); i++)
		{
			for (std::list<int>::const_iterator iter = hash_chain.table_[i].begin(); iter != hash_chain.table_[i].end(); iter++)
			{
				out_file << temp_file[*iter] << ',' << i << '\n';
			}
		}
	}
	// else if (config.hash_model_ == "pairwise")
	// {
	// 	printf("enter pairwise\n");
	// 	HashChainPairwise hash_chain(config.hash_prime_, 25);
	// 	for (int i = 0; i < permutation_vec.size(); i++)
	// 	{
	// 		hash_chain.Insert(permutation_vec[i]);
	// 	}
	// 	for (int i = 0; i < hash_chain.table_.size(); i++)
	// 	{
	// 		for (std::list<int>::const_iterator iter = hash_chain.table_[i].begin(); iter != hash_chain.table_[i].end(); iter++)
	// 		{
	// 			out_file << temp_file[*iter] << ',' << i << '\n';
	// 		}
	// 	}
	// }
	else if (config.hash_model_ == "4wise")
	{
		HashChainKwise hash_chain(4, config.hash_prime_);
		int cum_size = 0;
		for (int i = 0; i < permutation_vec.size(); i++)
		{
			hash_chain.Insert(permutation_vec[i]);
		}
		int bucket_idx = 0;
		for (int i = 0; i < hash_chain.table_.size(); i++)
		{
			if (!hash_chain.table_[i].empty())
			{
				for (std::list<int>::const_iterator iter = hash_chain.table_[i].begin(); iter != hash_chain.table_[i].end(); iter++)
				{
					out_file << temp_file[*iter] << ',' << bucket_idx << '\n';
				}
				slot_size_file << cum_size << '\n';
				cum_size += hash_chain.table_[i].size();
				bucket_idx ++;
			}
		}
		slot_size_file << cum_size << '\n';
		hash_chain.DisplayTable();
	}
	else
	{
		printf("Hash model is incorrect!\n");
	}
	out_file.close();
	slot_size_file.close();
}

void ReconstructFile(const std::string in_file_name, const Config& config)
{
    std::vector<std::string> temp_file;
	std::string out_file_name = config.prefix_ + config.dataset_ + "_reconstruct.csv", first_line, line_str;
    std::ifstream in_file(in_file_name, std::ios::in);
	std::ofstream out_file(out_file_name, std::ios::out);
    getline(in_file, first_line);
	out_file << first_line << '\n';
    int count = 1;
    while (getline(in_file, line_str))
    {
        //printf("%d\n", count);
        count++;
        //cout << line_str << endl; 
        temp_file.push_back(line_str);
    }
    in_file.close();
    for (int i = 0; i < temp_file.size(); i++)
    {
        out_file << temp_file[i] << '\n';
    }
	out_file.close();
}

void SaveSerializeFile(const std::string in_file_name, std::vector<std::vector<double>>& data_mat, const Config& config)
{
    ReadFile(in_file_name, data_mat);
	std::string file_int = config.prefix_ + config.dataset_ + "_serialize.file";
    save_serialized_graph(file_int, data_mat);
}

void LoadSerializeFile(const std::string in_file_name, std::vector<std::vector<double>>& data_mat, const Config& config)
{
    std::string file_int = config.prefix_ + config.dataset_ + "_serialize.file";
	load_serialized_graph(file_int, data_mat);
}

// void VarianceFilterQuery(const std::vector<std::vector<double>>& data_mat, const Config& config)
// {
// 	int num_feature = data_mat.size(), num_element = data_mat[0].size(), initial_sample_size = config.initial_size_;
//     printf("NumofFeature %d NumofRecord %d\n", num_feature, num_element);
// 	double pf = 1.0 / double(num_element), exact_time = 0.0, random_time = 0.0, colt_time = 0.0, mom_time = 0.0;
// 	RelativeError random_relative_error(0.0, 0.0), colt_relative_error(0.0, 0.0), mom_relative_error(0.0, 0.0), temp_relative_error;  
// 	Measure random_measure(0.0, 0.0, 0.0, 0.0, 0.0), colt_measure(0.0, 0.0, 0.0, 0.0, 0.0), mom_measure(0.0, 0.0, 0.0, 0.0, 0.0), temp_measure;
// 	printf("Initial Sample Size %d\nThreshold_%f\n", initial_sample_size, config.threshold_);
// 	std::chrono::high_resolution_clock::time_point start, end;
// 	std::chrono::duration<double> duration;
//     for (int t = 0; t < config.query_num_; t++)
//     {	
// 		std::vector<int> exact_variance_filter, approx_variance_filter;
// 		std::vector<std::pair<int, double>> exact_filter_pair, approx_filter_pair;
// 		start = std::chrono::high_resolution_clock::now();
// 		ExactVarianceFilter(config.threshold_, data_mat, exact_variance_filter, exact_filter_pair);
// 		end = std::chrono::high_resolution_clock::now();
// 		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
// 		printf("ExactVarianceFilterTime: %f\n\n", duration.count());
// 		exact_time += duration.count();

// 		start = std::chrono::high_resolution_clock::now();
// 		ApproxVarianceFilter("random", config.threshold_, initial_sample_size, data_mat, approx_variance_filter, pf, config.epsilon_, approx_filter_pair);
// 		end = std::chrono::high_resolution_clock::now();
// 		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
// 		temp_measure = Measure(Precision(exact_variance_filter, approx_variance_filter), Recall(exact_variance_filter, approx_variance_filter), F1Measure(exact_variance_filter, approx_variance_filter), 0.0, AccuracyFilter(num_feature, exact_variance_filter, approx_variance_filter));
// 		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f Acc %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.accuracy_);
// 		printf("ApproxVarianceFilterRandomTime: %f\n\n", duration.count());
// 		random_time += duration.count();
// 		temp_relative_error = FilterRelativeError(exact_filter_pair, approx_filter_pair);
//         random_relative_error += temp_relative_error;
// 		random_measure += temp_measure;

// 		start = std::chrono::high_resolution_clock::now();
// 		ApproxVarianceFilter("colt", config.threshold_, initial_sample_size, data_mat, approx_variance_filter, pf, config.epsilon_, approx_filter_pair);
// 		end = std::chrono::high_resolution_clock::now();
// 		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
// 		temp_measure = Measure(Precision(exact_variance_filter, approx_variance_filter), Recall(exact_variance_filter, approx_variance_filter), F1Measure(exact_variance_filter, approx_variance_filter), 0.0, AccuracyFilter(num_feature, exact_variance_filter, approx_variance_filter));
// 		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f Acc %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.accuracy_);
// 		printf("ApproxVarianceFilterRandomColtTime: %f\n\n", duration.count());
// 		colt_time += duration.count();
// 		temp_relative_error = FilterRelativeError(exact_filter_pair, approx_filter_pair);
//         colt_relative_error += temp_relative_error;
// 		colt_measure += temp_measure;

// 		start = std::chrono::high_resolution_clock::now();
// 		VarianceFilterMoMSimple(config.ratio_num_block_, config.threshold_, initial_sample_size, data_mat, approx_variance_filter, pf, config.epsilon_, approx_filter_pair);
// 		end = std::chrono::high_resolution_clock::now();
// 		duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
// 		temp_measure = Measure(Precision(exact_variance_filter, approx_variance_filter), Recall(exact_variance_filter, approx_variance_filter), F1Measure(exact_variance_filter, approx_variance_filter), 0.0, AccuracyFilter(num_feature, exact_variance_filter, approx_variance_filter));
// 		printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f Acc %1.4f\n", config.epsilon_, duration.count(), temp_measure.precision_, temp_measure.recall_, temp_measure.f1_, temp_measure.accuracy_);
// 		printf("ApproxVarianceFilterMoMUpdate3Time: %f\n\n", duration.count());
// 		mom_time += duration.count();
// 		temp_relative_error = FilterRelativeError(exact_filter_pair, approx_filter_pair);
//         mom_relative_error += temp_relative_error;
// 		mom_measure += temp_measure;
// 	}
// 	exact_time /= config.query_num_;
// 	printf("Exact %f\n", exact_time); 

// 	random_time /= config.query_num_;
// 	random_relative_error /= config.query_num_;
// 	random_measure /= config.query_num_;
// 	printf("Random %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f Acc %1.4f\n", random_time, random_relative_error.max_, random_relative_error.avg_, random_measure.precision_, random_measure.recall_, random_measure.f1_, random_measure.accuracy_);

// 	colt_time /= config.query_num_;
// 	colt_relative_error /= config.query_num_;
// 	colt_measure /= config.query_num_;
// 	printf("Colt   %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f Acc %1.4f\n", colt_time, colt_relative_error.max_, colt_relative_error.avg_, colt_measure.precision_, colt_measure.recall_, colt_measure.f1_, colt_measure.accuracy_);

// 	mom_time /= config.query_num_;
// 	mom_relative_error /= config.query_num_;
// 	mom_measure /= config.query_num_;
// 	printf("MoM    %f MRE %f ARE %f P %1.4f R %1.4f F1 %1.4f Acc %1.4f\n", mom_time, mom_relative_error.max_, mom_relative_error.avg_, mom_measure.precision_, mom_measure.recall_, mom_measure.f1_, mom_measure.accuracy_);
// }

// void NewTry(const std::vector<std::vector<double>>& data_mat, const Config& config)
// {
// 	int num_feature = data_mat.size(), num_element = data_mat[0].size(), sample_size = config.initial_size_;
//     printf("NumofFeature %d NumofRecord %d\n", num_feature, num_element);
// 	double pf = 1.0 / double(num_element) / num_feature, exact_time = 0.0, random_time = 0.0, colt_time = 0.0, mom_time = 0.0;
    
// 	exact_time /= config.query_num_;
// 	printf("Exact %f\n", exact_time); 
// 	while (sample_size < 1000000)
// 	{
// 		std::vector<double> tmp;
// 		double avg = 0.0, gamma4, c, se, sigma_hat, sum4 = 0.0, sum2 = 0.0, z=6.0;
// 		for (int i = 0; i < sample_size; i++)
// 		{
// 			double value = data_mat[61][rand() % num_element];
// 			avg += value;
// 			tmp.push_back(value);
// 		}
// 		avg /= sample_size;
// 		for (int i = 0; i < tmp.size(); i++)
// 		{
// 			sum4 += pow(tmp[i] - avg, 4.0);
// 			sum2 += pow(tmp[i] - avg, 2.0);
// 		}
// 		gamma4 = sample_size * sum4 / pow(sum2, 2.0);
// 		sigma_hat = sum2 / (sample_size - 1.0);
// 		c = sample_size / (sample_size - z);
// 		se = c * sqrt(gamma4 * (sample_size - 3.0) / sample_size / (sample_size - 1.0));
// 		printf("sampleSize %d avg %f sigma_hat %f gamma4 %f c %f se %f low %f up %f\n", sample_size, avg, sigma_hat, gamma4, c, se, exp(log(c*sigma_hat)-z*se), exp(log(c*sigma_hat)+z*se));
// 		sample_size *= 2;
// 	}
// }
#endif
// start = std::chrono::high_resolution_clock::now();
		// ApproxVarianceFilter("random_update1", config.threshold_, initial_sample_size, data_mat, approx_variance_filter, pf, config.epsilon_);
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f Acc %1.4f\n", config.epsilon_, duration.count(), Precision(exact_variance_filter, approx_variance_filter), Recall(exact_variance_filter, approx_variance_filter), F1Measure(exact_variance_filter, approx_variance_filter), AccuracyFilter(num_feature, exact_variance_filter, approx_variance_filter));
		// printf("ApproxVarianceFilterRandomUpdate1Time: %f\n", duration.count());
		// random_update1_time += duration.count();
		
		// start = std::chrono::high_resolution_clock::now();
		// ApproxVarianceFilter("stratified", config.threshold_, initial_sample_size, data_mat, approx_variance_filter, pf, config.epsilon_);
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f NDCG %1.4f\n", config.epsilon_, duration.count(), Precision(exact_variance_filter, approx_variance_filter), Recall(exact_variance_filter, approx_variance_filter), F1Measure(exact_variance_filter, approx_variance_filter), AccuracyFilter(num_feature, exact_variance_filter, approx_variance_filter));
		// printf("ApproxVarianceFilterStratifiedTime: %f\n", duration.count());
		// stratified_time += duration.count();

		// start = std::chrono::high_resolution_clock::now();
		// ApproxVarianceFilter("batch", config.threshold_, initial_sample_size, data_mat, approx_variance_filter, pf, config.epsilon_);
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f Acc %1.4f\n", config.epsilon_, duration.count(), Precision(exact_variance_filter, approx_variance_filter), Recall(exact_variance_filter, approx_variance_filter), F1Measure(exact_variance_filter, approx_variance_filter), AccuracyFilter(num_feature, exact_variance_filter, approx_variance_filter));
		// printf("ApproxVarianceFilterBatchTime: %f\n", duration.count());
		// batch_time += duration.count();

		// start = std::chrono::high_resolution_clock::now();
		// ApproxVarianceFilter("batch_update1", config.threshold_, initial_sample_size, data_mat, approx_variance_filter, pf, config.epsilon_);
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f Acc %1.4f\n", config.epsilon_, duration.count(), Precision(exact_variance_filter, approx_variance_filter), Recall(exact_variance_filter, approx_variance_filter), F1Measure(exact_variance_filter, approx_variance_filter), AccuracyFilter(num_feature, exact_variance_filter, approx_variance_filter));
		// printf("ApproxVarianceFilterBatchUpdate1Time: %f\n", duration.count());
		// batch_update1_time += duration.count();

		// start = std::chrono::high_resolution_clock::now();
		// ApproxVarianceFilter("batch_random", config.threshold_, initial_sample_size, data_mat, approx_variance_filter, pf, config.epsilon_);
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f Acc %1.4f\n", config.epsilon_, duration.count(), Precision(exact_variance_filter, approx_variance_filter), Recall(exact_variance_filter, approx_variance_filter), F1Measure(exact_variance_filter, approx_variance_filter), AccuracyFilter(num_feature, exact_variance_filter, approx_variance_filter));
		// printf("ApproxVarianceFilterBatchRandomTime: %f\n", duration.count());
		// batch_random_time += duration.count();

		// start = std::chrono::high_resolution_clock::now();
		// ApproxVarianceFilter("batch_random_update1", config.threshold_, initial_sample_size, data_mat, approx_variance_filter, pf, config.epsilon_);
		// end = std::chrono::high_resolution_clock::now();
		// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		// printf("Approx%1.4f %f P %1.4f R %1.4f F1 %1.4f Acc %1.4f\n", config.epsilon_, duration.count(), Precision(exact_variance_filter, approx_variance_filter), Recall(exact_variance_filter, approx_variance_filter), F1Measure(exact_variance_filter, approx_variance_filter), AccuracyFilter(num_feature, exact_variance_filter, approx_variance_filter));
		// printf("ApproxVarianceFilterBatchRandomUpdate1Time: %f\n", duration.count());
		// batch_random_update1_time += duration.count();