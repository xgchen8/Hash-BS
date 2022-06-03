#ifndef	KWISE_HASH_H
#define KWISE_HASH_H
#include <stdlib.h>
#include <bits/stdc++.h>
#include <vector>
class KwiseHash
{
private:
    std::vector<int> a_vec_;
    int prime_ = 60000011;
public:
    KwiseHash()
    {}
    KwiseHash(const int k, const int prime): prime_(prime)
    {
        printf("%d-wise hashing! prime %d\n", k, prime_);
        a_vec_.resize(k);
        for (int i = 0; i < a_vec_.size(); i++)
        {
            if (i == a_vec_.size() - 1)
            {
                a_vec_[i] = rand() % (prime_ - 1) + 1;
            }
            else
            {
                a_vec_[i] = rand() % prime_;
            }
            printf("a_%d=%d ", i, a_vec_[i]);
        }
        printf("\n");
    }
    ~KwiseHash()
    {}
    int HashValue(const int key) 
    {
        long long hash_value = 0, temp_key = 1;
        for (int i = 0; i < a_vec_.size(); i++)
        {
            hash_value = (hash_value + (a_vec_[i] * temp_key)) % prime_;
            temp_key = (temp_key * key) % prime_;
        }
        // hash_value = hash_value % (1 << hash_bit_length_);
        if (hash_value < 0)
        {
            printf("key %d p %d hash_value %lld\n", key, prime_, hash_value);
        }
        return int(hash_value);
    }
};
#endif