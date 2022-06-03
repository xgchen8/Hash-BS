// config.h
#ifndef	CONFIG_H
#define CONFIG_H
#include <string>
#include "pred_pair.h"
class Config 
{
public:
    std::string action_ = "", cache_miss_mode_ = "", prefix_ = "", dataset_ = "cdc82_01", hash_model_ = "tabulation";
    PredPair pred_;
    double epsilon_ = 0.1, threshold_ = 0.2, ratio_num_block_ = 4.5, empirical_threshold_ = 0.0, abs_error_ = 0.0;
    int topk_ = 4, initial_size_ = 512, query_num_ = 1, num_test_feature_ = 10, hash_prime_ = 60000011, block_size_ = 32;
    bool hashing_ = true;
    Config(): pred_(0, -1.0)
    {}
};
#endif