#ifndef	HASH_CHAIN_H
#define HASH_CHAIN_H
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include "tabulation.h"
class HashChain: public Tabulation
{
private:
    int size_, count_;               // count: number of data
public:
    std::vector<std::list<int>> table_;            // hash table with linked list
    HashChain() {};                       
    HashChain(const int key_length, const int hash_bit_length): Tabulation(key_length, hash_bit_length), size_(1 << hash_bit_length), count_(0)
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
        for (int i = 0; i < table_.size(); i++) 
        {
            printf("Slot#%d: ", i);
            for (std::list<int>::const_iterator iter = table_[i].begin(); iter != table_[i].end(); iter++) 
            {
                printf("%d ", *iter);
            }
            printf("\n");
        }
        printf("\n");
    }         
    ~HashChain() {};         
};
#endif
// #include <iostream>
// #include <vector>
// #include <list>
// #include <string>

// using std::vector;
// using std::list;
// using std::string;
// using std::cout;
// using std::endl;

// struct dict{                        // self-defined dictionary
//     string key;                     //  key  for Name (eg:Jordan)
//     string value;                   // value for Team (eg:Bulls)
//     dict():key(""),value(""){};     
//     dict(string Key, string Value):key(Key),value(Value){};  
// };

// class HashChain_std{
// private:
//     int size,                // size of table
//         count;               // count: number of data

//     vector<list<dict> > table;            // hash table with linked list

//     int PreHashing(string key_str);       // turn string_type_key to int_type_key
//     int HashFunction(string key_str);     // using Division method

// public:
//     HashChain_std(){};                       
//     HashChain_std(int m):size(m),count(0){   
//         table.resize(size);               // allocate memory for each slot
//     }

//     void Insert(dict data);               
//     void Delete(string key);              
//     string Search(string key);            
//     void DisplayTable();                  
// };

// string HashChain_std::Search(string key_str){
//     // two steps: 1. get index from hash function
//     //            2. traversal in linked list
//     int index = HashFunction(key_str);
//     for (list<dict>::iterator iter = table[index].begin(); iter != table[index].end(); iter++) {
//         if ((*iter).key == key_str) {
//             return (*iter).value;
//         }
//     }
//     return "...\nno such data";
// }

// void HashChain_std::Delete(string key_str){
//     // two steps: 1. get index from hash function
//     //            2. traversal in linked list
//     int index = HashFunction(key_str);
//     for (list<dict>::iterator iter = table[index].begin(); iter != table[index].end(); iter++) {
//         if ((*iter).key == key_str) {
//             table[index].erase(iter);
//         }
//     }
// }

// void HashChain_std::Insert(dict data){
//     // two steps: 1. get index from hash function
//     //            2. insert data at the front of linked list
//     int index = HashFunction(data.key);
//     table[index].push_front(data);             
// }

// int HashChain_std::PreHashing(string key_str){
//     // if   key_str = Jordan, exp = 9
//     // then key_int = ASCII(J)*9^5+ASCII(o)*9^4+ASCII(r)*9^3
//     //               +ASCII(d)*9^2+ASCII(a)*9^1+ASCII(n)*9^0

//     int exp = 9,        // choose randomly 
//         key_int = 0,
//         p = 1;

//     for (int i = (int)key_str.size()-1; i >= 0; i--) {
//         key_int += key_str[i]*p;
//         p *= exp;
//     }
//     return key_int;
// }

// int HashChain_std::HashFunction(string key_str){

//     return (PreHashing(key_str) % this->size);     // Division method
// }

// void HashChain_std::DisplayTable(){

//     for (int i = 0; i < table.size(); i++) {
//         cout << "slot#" << i << ": ";
//         for (list<dict>::iterator iter = table[i].begin(); iter != table[i].end(); iter++) {
//             cout << "(" << (*iter).key << "," << (*iter).value << ") ";
//         }
//         cout << endl;
//     }
//     cout << endl;
// }
