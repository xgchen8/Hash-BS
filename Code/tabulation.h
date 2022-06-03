#ifndef	TABULATION_H
#define TABULATION_H
#include <stdlib.h>
#include <bits/stdc++.h>
#include <vector>
static int get32rand() {
    return rand();
    // (((unsigned int) rand() <<  0) & 0x0000FFFFull) |
    // (((unsigned int) rand() << 16) & 0x7FFF0000ull);
}
class Tabulation
{
private:
    int key_length_ = 4, hash_bit_length_ = 25;
    std::vector<std::vector<int>> hash_table_;
public:
    Tabulation()
    {}
    Tabulation(int key_length, int hash_bit_length): key_length_(key_length), hash_bit_length_(hash_bit_length)
    {
        hash_table_.resize(key_length_, std::vector<int>(1 << CHAR_BIT));
        printf("Init %d %d\n", int(hash_table_.size()), int(hash_table_[0].size()));
        for (int i = 0 ; i < key_length_ ; i++) 
        {
		    for (int j = 0 ; j < (1 << CHAR_BIT); j++) 
            {
			    hash_table_[i][j] = get32rand() % (1 << hash_bit_length_);
                if (get32rand() % (1 << hash_bit_length_) < 0)
                {
                    printf("Error %d %d\n", get32rand(), (1 << hash_bit_length_));
                }
		    }
        }
        printf("numTable %d tableSize %d\n", int(hash_table_.size()), int(hash_table_[0].size()));
    }
    ~Tabulation()
    {}
    int HashValue(const int key) 
    {
        int hash_value = 0, temp_key = key;
        for (int i = 0; i < key_length_; i++)
        {
            hash_value ^= hash_table_[i][temp_key % (1 << CHAR_BIT)];
            if (key == 55601 || key == 55602)
            {
                printf("%d %d %d\n", key, temp_key % (1 << CHAR_BIT), hash_table_[i][temp_key % (1 << CHAR_BIT)]);
            }
            temp_key = temp_key >> CHAR_BIT;
        }
        if (key == 55601 || key == 55602)
        {
            printf("key %d\n", hash_value);
        }
        if (hash_value < 0)
        {
            printf("NegativeValue %d key %d\n", hash_value, key);
        }
        return hash_value;
    }
};
#endif