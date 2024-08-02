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
        str.pop_back();
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
    for (int i = stoi(argv[4]); i <= stoi(argv[5]); i += 1)
        kset.push_back(i);


    string dataset = argv[2];
    string filename = "/home/zqf/dataset/" + dataset + ".txt";
    fileset.push_back(filename); 


    int def_alg = stoi(argv[1]); // 0:DSA 1:DSA-topc 2:DSA-query_search 3:MTEnum 4:MTEnum \ MTHeur

    ofstream outfile;
    outfile.open(argv[3]);

    if (def_alg == 0)
    {
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
                int cnt = graph->find_min_k_truss(start, timeflag, "IE_ENUM");//IE_ENUM
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
                int cnt = graph->find_c_min_k_truss("IE_ENUM", start, timeflag, c);

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
                int cnt = graph->find_min_k_truss_query(start, timeflag, "IE_ENUM", q_vertex);
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
