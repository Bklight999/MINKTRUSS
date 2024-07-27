// Branch and bound for min_k_truss.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include"Utility.h"
#include"Graph.h"
#include"graph_enum.h"
#include"BnB1.h"
#include"ktb_enum.h"
#pragma warning(disable: 4996)

vector<std::string> readtxt_name(string file)
{
    ifstream infile;
    infile.open(file);
    cout << "infile.is_open() is : " << infile.is_open() << endl;

    if (infile.is_open() == 0) {
        cout << "file not existed" << endl;
    }
    vector<std::string> txt_name;
    string str;

    while (getline(infile, str))
    {
        str.pop_back();//getline reads line breaks, whereas file names in linux will report an error even if there is just one extra line break
        if (str == "EOF")break;
        txt_name.push_back(str);

    }
    infile.close();
    return txt_name;
}

int K;

vector<string> fileset;


int main(int argc, char* argv[]) 
{
    vector<int> kset;
    int k_l = stoi(argv[4]);
    int k_u = stoi(argv[5]);

    for (int i = k_l; i <= k_u; i += 1)
        kset.push_back(i);

    string filename = argv[2];// e.g. path/to/Skitter.txt"
    fileset.clear();
    fileset.push_back(filename);

    int def_alg = stoi(argv[1]); // 0:DSA 1:DSA-topc 2:DSA-query_search 3:ENUM+Heu_Add 4:ENUM+Heu_Del 5:ENUM  

    ofstream outfile;
    string out_file_name = argv[3]; // e.g. path/to/output.txt"

    if (def_alg == 0)
    {
       outfile.open(out_file_name);
        //outfile.open(argv[2]);

        for (int i = 0; i < fileset.size(); i++)
        {
            for (int turn = 0; turn < kset.size(); turn++)
            {
                K = kset[turn];
                string filename = fileset[i];

                clock_t start, finish;

                Graph* graph = new Graph(filename.c_str(), K);  

                bool timeflag = true;
                start = clock();
                int cnt = graph->find_min_k_truss(start, timeflag, "BKEPLEX");//IE_ENUM BNB BKEPLEX
                if (!timeflag) printf("Time exceeds!\n");
                else
                {
                    finish = clock();
                    printf("The running time is: %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
                }

                if (timeflag)
                    outfile << filename.substr(23) << " " << "BNB" << " " << "K = " << K << " " << (double)(finish - start) / CLOCKS_PER_SEC << "s" << " " << "size is" << " " << cnt << endl;
                else
                    outfile << filename.substr(23) << " " << "K = " << K << " " << "size is" << " " << "Time exceeds\n" << endl;
            }
            outfile << endl << endl;
        }
        outfile.close();
    }


    if (def_alg == 1)
    {
        outfile.open(out_file_name);
        int c = stoi(argv[6]); 

        for (int i = 0; i < fileset.size(); i++)
        {
            for (int turn = 0; turn < kset.size(); turn++)
            {
                K = kset[turn];
                string filename = fileset[i];

                clock_t start, finish;
                start = clock();
                Graph* graph = new Graph(filename.c_str(), K);  

                bool timeflag = true;
                int cnt = graph->find_c_min_k_truss("IE_ENUM", start, timeflag, c);//IE_ENUM

                if (!timeflag) printf("Time exceeds\n");
                else
                {
                    finish = clock();
                    printf("The running time is: %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
                }

                if (timeflag)
                    outfile << filename.substr(23) << " " << "BNB" << " " << "K = " << K << " " << (double)(finish - start) / CLOCKS_PER_SEC << "s" << " " << "size is" << " " << cnt << endl;
                else
                    outfile << filename.substr(23) << " " << "K = " << K << " " << "size is" << " " << "Time exceeds\n" << endl;
            }
            outfile << endl << endl;
        }
        outfile.close();
    }

    if (def_alg == 2)
    {
        outfile.open(out_file_name);
        int q_vertex = stoi(argv[6]);

        for (int i = 0; i < fileset.size(); i++)
        {
            for (int turn = 0; turn < kset.size(); turn++)
            {
                K = kset[turn];
                string filename = fileset[i];

                clock_t start, finish;
                start = clock();
                Graph* graph = new Graph(filename.c_str(), K); 

                bool timeflag = true;

                //srand((unsigned)time(NULL));
                int cnt = graph->find_min_k_truss_query(start, timeflag, "BNB", q_vertex);
                if (!timeflag) printf("Time exceeds\n");
                else
                {
                    finish = clock();
                    printf("The running time is: %lf s\n\n", (double)(finish - start) / CLOCKS_PER_SEC);
                }

                outfile << "Query vertex = " << q_vertex << endl;


                if (timeflag)
                    outfile << filename.substr(23) << " " << "BNB" << " " << "K = " << K << " " << (double)(finish - start) / CLOCKS_PER_SEC << "s" << " " << "size is" << " " << cnt << endl;
                else
                    outfile << filename.substr(23) << " " << "K = " << K << " " << "size is" << " " << "Time exceeds\n" << endl;
            }

            outfile << endl << endl;
        }
        outfile.close();
    }

    if (def_alg == 3)
    {
        outfile.open(out_file_name);

        for (int i = 0; i < fileset.size(); i++)
        {
            for (int turn = 0; turn < kset.size(); turn++)
            {
                K = kset[turn];
                string filename = fileset[i];

                clock_t start, finish;
                start = clock();
                GraphE* graphe = new GraphE(filename.c_str(), K);

                bool timeflag = true;

                int cnt = graphe->find_min_k_truss_add(start, timeflag);

                if (!timeflag) printf("Time exceeds\n");
                else
                {
                    finish = clock();
                    printf("The running time is: %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
                }

                if (timeflag)
                    outfile << filename.substr(23) << " " << "BNB" << " " << "K = " << K << " " << (double)(finish - start) / CLOCKS_PER_SEC << "s" << " " << "size is" << " " << cnt << endl;
                else
                    outfile << filename.substr(23) << " " << "K = " << K << " " << "size is" << " " << "Time exceeds\n" << endl;
            }
            outfile << endl << endl;
        }
        outfile.close();
    }

    
    if (def_alg == 4)
    {
        outfile.open(out_file_name);

        for (int i = 0; i < fileset.size(); i++)
        {
            for (int turn = 0; turn < kset.size(); turn++)
            {
                K = kset[turn];
                string filename = fileset[i];

                clock_t start, finish;
                start = clock();
                GraphE* graphe = new GraphE(filename.c_str(), K);

                bool timeflag = true;

                int cnt = graphe->find_min_k_truss_enum(start, timeflag);

                if (!timeflag) printf("Time exceeds\n");
                else
                {
                    finish = clock();
                    printf("The running time is: %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
                }

                if (timeflag)
                    outfile << filename.substr(23) << " " << "BNB" << " " << "K = " << K << " " << (double)(finish - start) / CLOCKS_PER_SEC << "s" << " " << "size is" << " " << cnt << endl;
                else
                    outfile << filename.substr(23) << " " << "K = " << K << " " << "size is" << " " << "Time exceeds" << endl;
            }
            outfile << endl << endl;
        }
        outfile.close();
    }

    return 0;
}
