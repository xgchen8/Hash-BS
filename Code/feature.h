// feature.h
#ifndef FEATURE_H
#define FEATURE_H
#include <stdlib.h> 
class Feature
{
public:
    Feature(): idx_(0), low_(0.0), up_(0.0), avg_(0.0), est_(0.0), diff_(1.0), ratio_(0.0)// , is_bound_(true)
    {}
    Feature(const int& idx, const double& init_value): idx_(idx), low_(init_value), up_(init_value), avg_(init_value), est_(init_value), diff_(1.0)//, is_bound_(true)
    {
        if (up_ != 0.0)
        {   
            ratio_ = low_ / up_;
        }
        else
        {
            ratio_ = 0.0;
        }
    }
    // Feature(const int& idx, const double& low, const double& up): idx_(idx), low_(low), up_(up)//, is_bound_(true)
    // {
    //     avg_ = (low_ + up_) / 2.0;
    //     diff_ = abs(up_ - low_);
    //     if (up > 0)
    //     {   
    //         ratio_ = low_ / up_;
    //     }
    //     else
    //     {
    //         ratio_ = 0.0;
    //     }
    //     est_ = diff_;
    // }
    int idx_;
    double low_, up_, avg_, est_, diff_, ratio_;
    void Update(const double& low, const double& up, const double& sample_est)
    {
        low_ = low;
        up_ = up;
        avg_ = (low_ + up_) / 2.0;
        diff_ = abs(up_ - low_);
        if (up_ != 0)
        {   
            ratio_ = low_ / up_;
        }
        else
        {
            ratio_ = 0.0;
        }
        if (sample_est != low_)
        {
            est_ = avg_;
        }
        else
        {
            est_ = sample_est;
        }
    }
    // void Update(const double& low, const double& up)
    // {
    //     low_ = low;
    //     up_ = up;
    //     avg_ = (low_ + up_) / 2.0;
    //     diff_ = abs(up_ - low_);
    //     if (up > 0)
    //     {   
    //         ratio_ = low_ / up_;
    //     }
    //     else
    //     {
    //         ratio_ = 0.0;
    //     }
    //     est_ = diff_;
    // }
};
#endif