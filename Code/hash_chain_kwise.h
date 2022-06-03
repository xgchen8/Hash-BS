#ifndef	HASH_CHAIN_KWISE_H
#define HASH_CHAIN_KWISE_H
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include "kwise.h"
class HashChainKwise: public KwiseHash
{
private:
    int size_, count_;               // count: number of data
public:
    std::vector<std::list<int>> table_;            // hash table with linked list
    HashChainKwise() {};                       
    HashChainKwise(const int k, const int prime): KwiseHash(k, prime), size_(prime), count_(0)
    {   
        table_.resize(size_);              // allocate memory for each slot
    }
    void Insert(const int key)
    {// two steps: 1. get index from hash function
    //            2. insert data at the front of linked list
        int index = HashValue(key);
        if (index < 0)
        {
            printf("%d\n", index);
        }
        
        if (index > table_.size())
        {
            printf("tablesize %d index %d\n", int(table_.size()), index);
        }
        table_[index].push_front(key);  
    };               
    void Delete(const int key)
    {
        int index = HashValue(key);
        for (std::list<int>::const_iterator iter = table_[index].begin(); iter != table_[index].end(); iter++) 
        {
            if (*iter == key) 
            {
                table_[index].erase(iter);
            }
        }
    }              
    int Search(const int key)
    {
    // two steps: 1. get index from hash function
    //            2. traversal in linked list
        int index = HashValue(key);
        for (std::list<int>::const_iterator iter = table_[index].begin(); iter != table_[index].end(); iter++) 
        {
            if (*iter == key) 
            {
                return *iter;
            }
        }
        return -1;
    }     
    void DisplayTable()
    {
        // for (int i = 0; i < table_.size(); i++) 
        for (int i = 0; i < 100; i++) 
        {
            printf("Slot#%d: ", i);
            for (std::list<int>::const_iterator iter = table_[i].begin(); iter != table_[i].end(); iter++) 
            {
                printf("%d ", *iter);
            }
            printf("\n");
        }
        printf("table size %d\n", int(table_.size()));
    }         
    ~HashChainKwise() {};         
};
#endif