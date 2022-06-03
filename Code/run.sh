#!/bin/sh

timeEnem72="/home/MoMRandom/dataFiles/output021022"

# 59_0.13953488372093025838 101_0.33333333333333331483 88_0.00 28_0.0

for col in '101' # '59' '88' '28' # '59' '65' '88' '28'
do
    for i in '01' # '005' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
    do  
        for j in '512' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
        do
            for k in '04' # '16' '08' '02' '01' # '01' '02' '08' '16' '04' 
            do
                if ((${col} == '59')); then
                    stdbuf -oL ./test variance_topk --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 59 0.13953488372093025838 > ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}.txt
                elif ((${col} == '101')); then
                    stdbuf -oL ./test variance_topk --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 101 0.33333333333333331483 > ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}.txt
                elif ((${col} == '88')); then
                    stdbuf -oL ./test variance_topk --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 88 0.00 > ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}.txt
                else
                    stdbuf -oL ./test variance_topk --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 28 0.0 > ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}.txt
                fi
            done
        done
    done
done

# for col in '101' #'59' '88' '28' # '59' '65' '88' '28'
# do
#     for i in '01' # '005' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             for k in '04' # '16' '08' '04' '02' # '01' '02' '08' '16' '04' 
#             do
#                 if ((${col} == '101')); then
#                     perf record -o ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}mom.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u,page-faults:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_mom --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 101 0.33333333333333331483 > ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}mom.txt

#                     perf record -o ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}exact.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u,page-faults:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_exact --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 101 0.33333333333333331483 > ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}exact.txt
                    
#                     # perf record -o ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}random.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_random --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 101 0.33333333333333331483 > ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}random.txt

#                     # perf record -o ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}colt.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_colt --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 101 0.33333333333333331483 > ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}colt.txt

#                     perf record -o ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}mom.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_mom --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 101 0.33333333333333331483 > ${timeEnem72}${k}Top${k}enem72SelectCol${col}AbsError0${i}_initSize${j}mom.txt
#                 fi
#             done
#         done
#     done
# done

# 59_0.13953488372093025838 65_0.40000000000000002220 88_0.00 28_0.0
# 59_0.13953488372093025838 101_0.33333333333333331483 88_0.00 28_0.0

# for col in '101' # '59' # '65' '88' '28'
# do
#     for i in '025' '05' # '005' '01' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '128' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             if ((${col} == '59')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 59 0.13953488372093025838 > ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}.txt
#             elif ((${col} == '101')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 101 0.33333333333333331483 > ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}Test2.txt
#             elif ((${col} == '88')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 88 0.00 > ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}.txt
#             else
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 28 0.0 > ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}.txt
#             fi
#         done
#     done
# done

# for col in '101' # '65' '88' '28'
# do
#     for i in '01' # '005' '006' '008' '01' '025' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '128' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             if ((${col} == '101')); then
#                 perf record -o ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}Mom.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u,page-faults:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_mom --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 101 0.33333333333333331483 > ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}Mom.txt
                
#                 perf record -o ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}Exact.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u,page-faults:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_exact --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 101 0.33333333333333331483 > ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}Exact.txt

#                 perf record -o ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}Random.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u,page-faults:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_random --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 101 0.33333333333333331483 > ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}Random.txt

#                 perf record -o ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}Colt.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u,page-faults:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_colt --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 101 0.33333333333333331483 > ${timeEnem72}30enem72SelectCol${col}AbsError0${i}_initSize${j}Colt.txt
#             fi
#         done
#     done
# done

# stdbuf -oL ./test hash_file --prefix /home/Downloads/enem7/ --dataset enem72ori_01 --hash_model 4wise --hash_prime 49610017 > ${timeEnem72}01enem72_4wise_hashing.txt
# ./test save_serialize_file --prefix /home/Downloads/enem7/ --dataset enem72ori_01_4wise

timeOntime="/home/MoMRandom/dataFiles/output020911"
#1_0.00000000000000000000 32_0.0_0.5 53_1.0_1.0 34_0.07142857142857139685_0.5
# for col in '6' '43' '49' # '6' '55' '43' '49'
# do
#     for i in '01' # '0025' '005' '01' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             for k in '04' # '01' '02' '04' '08' '16' # '' 
#             do
#             # if ((${col} == '14')); then
#             #     stdbuf -oL ./test variance_topk --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 1 --initial_size ${j} --predicate 14 0.01639344262295079971 > ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 if ((${col} == '6')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 6 0.2727272727272727 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '55')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 55 0.10000000000000000555 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '43')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 43 0.0 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '32')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 32 0.0 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '48')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 48 0.0 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '53')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 53 1.0 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '45')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 45 0.07142857142857139685 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '49')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 49 0.0 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
                
#                 fi
#             done
#         done
#     done
# done

# for col in '55' # '6' '43' '49' # '6' '55' '43' '49'
# do
#     for i in '01' # '0025' '005' '01' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             for k in '04' #'01' '02' '08' # '' 
#             do
#                 if ((${col} == '55')); then
#                     perf record -o ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}exact.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_exact --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 55 0.10000000000000000555 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}exact.txt

#                     perf record -o ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}random.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_random --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 55 0.10000000000000000555 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}random.txt

#                     perf record -o ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}colt.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_colt --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 55 0.10000000000000000555 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}colt.txt

#                     perf record -o ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}mom.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_mom --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 55 0.10000000000000000555 > ${timeOntime}${k}Top${k}ontimeSelectCol${col}AbsError0${i}_initSize${j}mom.txt
#                 fi
#             done
#         done
#     done
# done

# 0.04232592357363166


# 14_0.01639344262295079971 6_0.2727272727272727 1_0.6666666666666666 43_0.0 49_0.0

# for col in '55' # '1' '43' '49'
# do
#     for i in '025' '05' # '0025' '005' '01' '025' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             if ((${col} == '14')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 1 --initial_size ${j} --predicate 14 0.01639344262295079971 > ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#             elif ((${col} == '55')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 55 0.10000000000000000555 > ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}Test2.txt
#             elif ((${col} == '1')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 1 0.6666666666666666 > ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#             elif ((${col} == '43')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 43 0.0 > ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#             else
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 49 0.0 > ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}.txt
#             fi
#         done
#     done
# done

# for col in '55' # '1' '43' '49'
# do
#     for i in '01' # '005' '01' '025' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             if ((${col} == '55')); then
#                 perf record -o ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}Exact.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_exact --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 55 0.10000000000000000555 > ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}Exact.txt

#                 perf record -o ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}Random.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_random --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 55 0.10000000000000000555 > ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}Random.txt

#                 perf record -o ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}Colt.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_colt --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 55 0.10000000000000000555 > ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}Colt.txt

#                 perf record -o ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}Mom.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_mom --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 55 0.10000000000000000555 > ${timeOntime}30ontimeSelectCol${col}AbsError0${i}_initSize${j}Mom.txt
#             fi
#         done
#     done
# done

# stdbuf -oL ./test hash_file --prefix /home/Downloads/ontime/ --dataset ontimeori_01 --epsilon 0.${i} --topk ${j} --hashing 0 --hash_model 4wise --hash_prime 199874833 > ${timeOntime}ontime_4wise_hashing.txt
# ./test save_serialize_file --prefix /home/Downloads/ontime/ --dataset ontimeori_01_4wise --epsilon 0.${i} --threshold 0.${j} --hashing 0 --initial_size 64 --query_num 10

Green="/home/MoMRandom/dataFiles/output020224"
# 7_0.24528301886792455710 12_0.40000000000000002220 11_0.04762362855855657717 0_0.5
# for col in '12' # '7' '12' '11' '0'
# do
#     for i in '01' # '0025' '005' '01' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '128' '512' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             for k in '01' '02' '08' '16' # '01' '02' '08' '16' '04' 
#             do
#                 if ((${col} == '7')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/local/green_tripdata/ --dataset greenori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 7 0.24528301886792455710 > ${Green}${k}Top${k}greenSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '12')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/local/green_tripdata/ --dataset greenori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 12 0.40000000000000002220 > ${Green}${k}Top${k}greenSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '11')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/local/green_tripdata/ --dataset greenori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 11 0.04762362855855657717 > ${Green}${k}Top${k}greenSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 else
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/local/green_tripdata/ --dataset greenori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 0 0.5 > ${Green}${k}Top${k}greenSelectCol${col}AbsError0${i}_initSize${j}.txt
#                 fi
#             done
#         done
#     done
# done

# for col in '12' # '11' '0'
# do
#     for i in '0025' # '005' '01' '025' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' '128' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             if ((${col} == '7')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/local/green_tripdata/ --dataset greenori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 7 0.24528301886792455710 > ${Green}30greenSelectCol${col}AbsError0${i}_initSize${j}.txt
#             elif ((${col} == '12')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/local/green_tripdata/ --dataset greenori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 12 0.40000000000000002220 > ${Green}30greenSelectCol${col}AbsError0${i}_initSize${j}.txt
#             elif ((${col} == '11')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/local/green_tripdata/ --dataset greenori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 11 0.04762362855855657717 > ${Green}30greenSelectCol${col}AbsError0${i}_initSize${j}.txt
#             else
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/local/green_tripdata/ --dataset greenori_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 0 0.5 > ${Green}30greenSelectCol${col}AbsError0${i}_initSize${j}.txt
#             fi
#         done
#     done
# done

# stdbuf -oL ./test hash_file --prefix /home/Downloads/local/green_tripdata/ --dataset greenori_01 --hash_model 4wise --hash_prime 81563653 > ${Green}green_4wise_hashing.txt
# ./test save_serialize_file --prefix /home/Downloads/local/green_tripdata/ --dataset greenori_01_4wise

# 7_0.24528301886792455710 12_0.40000000000000002220 11_0.04762362855855657717 0_0.5

TimeYellow="/home/MoMRandom/dataFiles/output012615"
# 16_0.17054869000000000279 12_0.07368001600000000095 12_0.07321661999999999615 14_0.99974348000000001768

# for col in '14' # '16' '12' '14'
# do
#     for i in '005' # '006' '008' '01' '025' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' #'128' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             if ((${col} == '16')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/yellow/ --dataset yellow_taxi_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 16 0.17054869000000000279 > ${TimeYellow}30yellowSelectCol${col}AbsError0${i}_initSize${j}.txt
#             elif ((${col} == '12')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/yellow/ --dataset yellow_taxi_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 12 0.07368001600000000095 > ${TimeYellow}30yellowSelectCol${col}025AbsError0${i}_initSize${j}.txt
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/yellow/ --dataset yellow_taxi_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 12 0.07321661999999999615 > ${TimeYellow}30yellowSelectCol${col}05AbsError0${i}_initSize${j}.txt
#             else
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/yellow/ --dataset yellow_taxi_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 14 0.99974348000000001768 > ${TimeYellow}30yellowSelectCol${col}AbsError0${i}_initSize${j}Swapoff.txt
#             fi
#         done
#     done
# done

# stdbuf -oL ./test hash_file --prefix /home/Downloads/yellow/ --dataset yellow_taxi_01 --hash_model 4wise --hash_prime 1173057943 > ${TimeYellow}07yellow_4wise_hashing.txt
# ./test save_serialize_file --prefix /home/Downloads/yellow/ --dataset yellow_taxi_01_4wise

TimePus7="/home/MoMRandom/dataFiles/output020912"
# 33_1.0 130_1.00000000000000000000 138_1.0 38_0.0 
# 2_0.14285714285714284921
# for col in '33' '138' '38' # '33' '130' '138' '38'
# do
#     for i in '01' # '0025' '005' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             for k in '04' # '01' '02' '08' '16' '04' 
#             do
#                 if ((${col} == '33')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 33 1.0 > ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '130')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 130 1.00000000000000000000 > ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '138')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 138 1.0 > ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}.txt
#                 else
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 38 0.0 > ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}.txt
#                 fi
#             done
#         done
#     done
# done

# for col in '130' #'33' '138' '38' # '33' '130' '138' '38'
# do
#     for i in  '01' # '005' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             for k in '04' # '01' '02' '08' '16' '04' 
#             do
#                 if ((${col} == '130')); then
#                     perf record -o ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}exact.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_exact --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 130 1.00000000000000000000 > ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}exact.txt

#                     perf record -o ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}random.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_random --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 130 1.00000000000000000000 > ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}random.txt

#                     perf record -o ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}colt.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_colt --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 130 1.00000000000000000000 > ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}colt.txt

#                     perf record -o ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}mom.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_mom --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 130 1.00000000000000000000 > ${TimePus7}${k}Top${k}pus7SelectCol${col}AbsError0${i}_initSize${j}mom.txt
#                 fi
#             done
#         done
#     done
# done

# 33_1.0 130_1.00000000000000000000 138_1.0 38_0.0 
# for col in '130' # '33' '2' '138' '38' # '16' '12' '14'
# do
#     for i in '0025' '005' '01' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '128' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             if ((${col} == '33')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 33 1.0 > ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}.txt
#             elif ((${col} == '130')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 130 1.00000000000000000000 > ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}.txt
#             elif ((${col} == '138')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 138 1.0 > ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}.txt
#             else
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 38 0.0 > ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}.txt
#             fi
#         done
#     done
# done

# for col in '130' # '33' '2' '138' '38' # '16' '12' '14'
# do
#     for i in '01' # '0025' '005' '01' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '128' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             if ((${col} == '130')); then
#                 perf record -o ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}Exact.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_exact --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 130 1.00000000000000000000 > ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}Exact.txt

#                 perf record -o ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}Random.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_random --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 130 1.00000000000000000000 > ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}Random.txt

#                 perf record -o ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}Colt.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_colt --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 130 1.00000000000000000000 > ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}Colt.txt

#                 perf record -o ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}Mom.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_mom --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 130 1.00000000000000000000 > ${TimePus7}39pus7SelectCol${col}AbsError0${i}_initSize${j}Mom.txt
#             fi
#         done
#     done
# done

# stdbuf -oL ./test hash_file --prefix /home/Downloads/5year/pus/ --dataset pus714_01 --hash_model 4wise --hash_prime 109869037 > ${TimePus7}39pus7_4wise_hashing.txt
# ./test save_serialize_file --prefix /home/Downloads/5year/pus/ --dataset pus714_01_4wise

TimeHus11="/home/MoMRandom/dataFiles/output020913"
# 27_0.5 42_0.33333333333333331483 46_0.0 34_0.0
# stdbuf -oL ./test hash_file --prefix /home/Downloads/5year/hus/ --dataset hus11_01 --hash_model 4wise --hash_prime 79728347 > ${TimeHus11}39hus11_4wise_hashing.txt
# ./test save_serialize_file --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise
# for col in '27' '46' '34' # '27' '42' '46' '34'
# do
#     for i in '01' # '0025' '005' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             for k in '04' # '04' '01' '02' '08' '16'
#             do
#                 if ((${col} == '27')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 27 0.5 > ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '42')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 42 0.33333333333333331483 > ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}.txt
#                 elif ((${col} == '46')); then
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 46 0.0 > ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}.txt
#                 else
#                     stdbuf -oL ./test variance_topk --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 34 0.0 > ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}.txt
#                 fi
#             done
#         done
#     done
# done

# for col in '42' # '27' '46' '34' # '27' '42' '46' '34'
# do
#     for i in '01' # '04' # '0025' '005' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             for k in '04' # '02' '08' 
#             do
#                 if ((${col} == '42')); then
#                     perf record -o ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}exact.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_exact --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 42 0.33333333333333331483 > ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}exact.txt

#                     perf record -o ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}random.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_random --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 42 0.33333333333333331483 > ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}random.txt

#                     perf record -o ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}colt.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_colt --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 42 0.33333333333333331483 > ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}colt.txt

#                     perf record -o ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}mom.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u stdbuf -oL ./test cache_miss --cache_miss_mode topk_mom --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --topk ${k} --hashing 0 --query_num 10 --initial_size ${j} --predicate 42 0.33333333333333331483 > ${TimeHus11}${k}Top${k}hus11SelectCol${col}AbsError0${i}_initSize${j}mom.txt
#                 fi
#             done
#         done
#     done
# done

# 27_0.5 42_0.33333333333333331483 46_0.0 34_0.0
# for col in '27' '46' '34' # '27' '42' '46' '34'
# do
#     for i in '01' # '0025' '005' '01' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '128' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             if ((${col} == '27')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 27 0.5 > ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}.txt
#             elif ((${col} == '42')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 42 0.33333333333333331483 > ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}.txt
#             elif ((${col} == '46')); then
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 46 0.0 > ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}.txt
#             else
#                 stdbuf -oL ./test absolute_metric --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 34 0.0 > ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}.txt
#             fi
#         done
#     done
# done

# for col in '42' # '27' '46' '34' # '27' '42' '46' '34'
# do
#     for i in '01' # '0025' '005' '01' '02' '025' '04' '05' # '2' '4' '8'  # '01' '025' '05' '1' '25' '5'
#     do  
#         for j in '512' # '128' # '16' '32' '64' '128' '256' '512' '1024' '2048' '4096' '8192' 
#         do
#             if ((${col} == '42')); then
#                 perf record -o ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}Exact.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_exact --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 42 0.33333333333333331483 > ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}Exact.txt

#                 perf record -o ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}Random.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_random --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 42 0.33333333333333331483 > ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}Random.txt

#                 perf record -o ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}Colt.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_colt --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 42 0.33333333333333331483 > ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}Colt.txt

#                 perf record -o ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}Mom.data -e cache-misses:u,branch-misses:u,LLC-load-misses:u,LLC-store-misses:u,dTLB-load-misses:u,dTLB-store-misses:u,branch-load-misses:u -c 10 stdbuf -oL ./test cache_miss --cache_miss_mode single_mom --prefix /home/Downloads/5year/hus/ --dataset hus11_01_4wise --abs_error 0.${i} --hashing 0 --query_num 10 --initial_size ${j} --predicate 42 0.33333333333333331483 > ${TimeHus11}31hus11SelectCol${col}AbsError0${i}_initSize${j}Mom.txt
#             fi
#         done
#     done
# done