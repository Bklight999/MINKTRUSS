# MINKTRUSS
This repository contains the codebase for the following manuscript: 
***Efficient Minimum k-Truss Search: A Decomposition-based Approach***

# Compile the code

```cpp
g++ main.cpp -o minktruss
```
This will generate an executable file **minktruss**

# Run the code

```cpp
./minktruss $alogirthm $input_file_path $output_file_path $k_lower $k_upper $query_id / the value of c (optional)
```
Note that, the algorithms mentioned in our paper are transfered to a series of numbers:

 - 0: MINKTRUSS-DSA
 - 1: MINKTRUSS-DSA-topc
 - 2: MINKTRUSS-DSA-query-based
 - 3: MINKTRUSS-ENUM
 - 4: MINKTRUSS-ENUM without MINKTRUSS-HEURISTIC

**For example:**
 if you want to run *MINKTRUSS--DSA* on *RGG* from *k=5 to k=15* and *c=10*, you can use the command:

```cpp
./minktruss 0 path/to/your_input_file path/to/your_output_file 5 15 10
```

 if you want to run *MINKTRUSS--DSA* on *Skitter* with *k=5* and *query_id= 537*, you can use the command:

```cpp
./minktruss 0 path/to/your_input_file path/to/your_output_file 5 5 537
```
