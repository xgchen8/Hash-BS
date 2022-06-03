// test.cpp
#include "config.h"
#include "query.h"
#include <chrono>
#include <random>
#include <string>

int main(int argc, char **argv)
{
	Config config;
	for (int i = 0; i < argc; i++) {
        std::string help_str = ""
			"main action [options]\n"
			"\n"
			"options: \n"
            "  --cache_miss_mode <cache miss mode>\n"
			"  --prefix <prefix>\n"
			"  --dataset <dataset>\n"
			"  --support_file <support_file>\n"
			"  --epsilon <epsilon>\n"
			"  --topk <topk>\n"
            "  --threshold <threshold>\n"
            "  --hashing <bool>\n"
            "  --initial_size <intial_size>\n"
            "  --query_num <query_num>\n"
            "  --hash_model <hash_model>\n"
            "  --hash_prime <prime number>\n"
            "  --metric <metric>\n"
            "  --block_size <block_size>";
        if (std::string(argv[i]) == "--help") {
            std::cout << help_str << std::endl;
            exit(0);
        }
	}
    config.action_ = argv[1];
    std::cout << "action: " << config.action_ << std::endl;
    for (int i = 0; i < argc; i++) {
        std::string arg = argv[i];
        if (std::string(argv[i]) == "--prefix") {
            config.prefix_ = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--cache_miss_mode") {
            config.cache_miss_mode_ = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--dataset") {
            config.dataset_ = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--epsilon") {
            config.epsilon_ = std::stod(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--abs_error") {
            config.abs_error_ = std::stod(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--topk"){
            config.topk_ = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--threshold"){
            config.threshold_ = std::stod(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--hashing"){
            if (std::string(argv[i + 1]) == "1")
            {
                config.hashing_ = true;
            }
            else if (std::string(argv[i + 1]) == "0")
            {
                config.hashing_ = false;
            }
            else
            {
                std::cerr << "command not recognize " << arg << std::endl;
            }
        }
        else if (std::string(argv[i]) == "--initial_size"){
            config.initial_size_ = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--query_num"){
            config.query_num_ = std::stoi(argv[i + 1]);
        }       
        else if (std::string(argv[i]) == "--hash_model"){
            config.hash_model_ = argv[i + 1];
            std::cout << config.hash_model_ << std::endl;
        }    
        else if (std::string(argv[i]) == "--hash_prime"){
            config.hash_prime_ = std::stoi(argv[i + 1]);
            std::cout << config.hash_prime_ << std::endl;
        }    
        // else if (std::string(argv[i]) == "--metric"){
        //     config.metric_ = argv[i + 1];
        //     std::cout << config.metric_ << std::endl;
        // }   
        else if (std::string(argv[i]) == "--block_size"){
            config.block_size_ = std::stoi(argv[i + 1]);
            std::cout << config.block_size_ << std::endl;
        }
        else if (std::string(argv[i]) == "--ratio_num_block"){
            config.ratio_num_block_ = std::stod(argv[i + 1]);
            std::cout << config.ratio_num_block_ << std::endl;
        }
        else if (std::string(argv[i]) == "--empirical_threshold"){
            config.empirical_threshold_ = std::stod(argv[i + 1]);
            std::cout << config.empirical_threshold_ << std::endl;
        }
        else if (std::string(argv[i]) == "--predicate"){
            config.pred_.col_ = std::stoi(argv[i + 1]);
            config.pred_.value_ = std::stod(argv[i + 2]);
        }
        else if (arg.substr(0, 2) == "--") {
            std::cerr << "command not recognize " << arg << std::endl;
            exit(1);
        }
    }
    // printf("col %d pred_value %f\n", config.pred_.col_, config.pred_.value_);
    srand((unsigned)time(NULL));
    seed[0] = (unsigned)time(NULL);
    seed[1] = (unsigned)time(NULL);
    // s[0] = (unsigned)time(NULL);
    // s[1] = (unsigned)time(NULL);
    // random_init();
    // Xoshiro256plus::setSeed((unsigned)time(NULL));
	std::ostringstream oss_data_mat, oss_cumulate_size;
	oss_data_mat << config.prefix_ << config.dataset_ << ".csv";
    oss_cumulate_size << config.prefix_ << config.dataset_ << "_cumulate_size.csv";
	std::vector<std::vector<double>> data_mat;
    std::vector<uint32_t> cum_size;
	
    if (config.action_ == "variance_topk")
    {
        ReadCumSize(oss_cumulate_size.str(), cum_size);
        LoadSerializeFile(oss_data_mat.str(), data_mat, config);
        
        if (config.pred_.value_ < 0)
        {
            VarianceTopkQuery(data_mat, cum_size, config);
        }
        else
        {
            printf("Predicate Mode\n");
            SelectAbsPredTopkQuery(data_mat, cum_size, config);
            // PredVarianceTopkQuery(data_mat, cum_size, config);
        }
    }
    else if (config.action_ == "hash_file")
    {
        printf("enter\n");
        HashFile(oss_data_mat.str(), config);
    }
    else if (config.action_ == "reconstruct_file")
    {
        ReconstructFile(oss_data_mat.str(), config);
    }
    // else if (config.action_ == "variance_filter")
    // {
    //     LoadSerializeFile(oss_data_mat.str(), data_mat, config);
    //     variance_filter_query(data_mat, config);
    // }
    else if (config.action_ == "save_serialize_file")
    {
        SaveSerializeFile(oss_data_mat.str(), data_mat, config);
    }
    else if (config.action_ == "compare_variance")
    {
        printf("enter compare variance!\n");
        LoadSerializeFile(oss_data_mat.str(), data_mat, config);
        CompareVariance(data_mat, config.topk_, config.initial_size_);
    }
    else if (config.action_ == "compare_average")
    {
        printf("enter compare average!\n");
        LoadSerializeFile(oss_data_mat.str(), data_mat, config);
        CompareAverage(data_mat, config.topk_, config.initial_size_);
    }
    else if (config.action_ == "epsilon_metric")
    {
        printf("epsilon approximate variance!\n");
        ReadCumSize(oss_cumulate_size.str(), cum_size);
        LoadSerializeFile(oss_data_mat.str(), data_mat, config);
        if (config.pred_.value_ < 0)
        {
            EpsilonMetric(data_mat, cum_size, config);
        }
        else
        {
            printf("col %d pred_value %f\n", config.pred_.col_, config.pred_.value_);
            EpsilonDeltaPredMetric(data_mat, cum_size, config);
        }
    }
    else if (config.action_ == "absolute_metric")
    {
        printf("absolute approximate variance!\n");
        ReadCumSize(oss_cumulate_size.str(), cum_size);
        LoadSerializeFile(oss_data_mat.str(), data_mat, config);
        printf("col %d pred_value %f\n", config.pred_.col_, config.pred_.value_);
        SelectAbsPredQuery(data_mat, cum_size, config);
    }
    else if (config.action_ == "cache_miss")
    {
        printf("cache miss!\n");
        ReadCumSize(oss_cumulate_size.str(), cum_size);
        LoadSerializeFile(oss_data_mat.str(), data_mat, config);
        printf("col %d pred_value %f\n", config.pred_.col_, config.pred_.value_);
        CacheMissStat(data_mat, cum_size, config);
    }
    else
    {
        printf("Command Error!\n");
    }
    return 0;
}