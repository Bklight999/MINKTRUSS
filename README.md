# MINKTRUSS
This repository contains the codebase for the following manuscript: 
***Efficient Minimum k-Truss Search: A Decomposition-based Approach***

# Dataset Preparation
The dataset dataset.zip can be unzipped directly.

| **Dataset** | **$n$**   | **$m$**    | **$d_{avg}$** | **$k_{max}$** |
|-------------|-----------|------------|---------------|---------------|
| Slashdot    | 70,068    | 358,647    | 10.24         | 35            |
| Epinions    | 75,879    | 405,740    | 10.6          | 33            |
| UCLA        | 20,453    | 747,604    | 73            | 66            |
| Harvard     | 15,126    | 824,617    | 109           | 42            |
| Youtube     | 1,134,890 | 2,987,624  | 2.6           | 19            |
| Google      | 875,713   | 4,322,051  | 9.87          | 44            |
| Wiki        | 2,394,385 | 4,659,565  | 1.9           | 53            |
| Flixster    | 2,523,386 | 7,918,801  | 6.28          | 30            |
| Skitter     | 1,696,415 | 11,095,298 | 13.08         | 68            |
| RGG         | 16,777,216| 132,557,200| 15            | 21            |



Other datasets can be downloaded from the following links:

 - **SNAP**: http://snap.stanford.edu
 - **Networkrepository**: https://networkrepository.com/

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
./minktruss 0 path/to/RGG.txt path/to/your_output_file 5 15 10
```

 if you want to run *MINKTRUSS--DSA* on *Skitter* with *k=5* and *query_id= 537*, you can use the command:

```cpp
./minktruss 0 path/to/Skitter.txt path/to/your_output_file 5 5 537
```
# Long-term Preservation
This repository is maintained by Qifan Zhang (bklight999@gmail.com).
