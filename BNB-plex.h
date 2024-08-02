﻿#pragma once
#ifndef _BNBPLEX_H_
#define _BNBPLEX_H_

#include"Utility.h"
#include"Timer.h"

class BnBplex
{
private:
	int n;
	int m;
	int S;
	int must_include;
	int* pstart;
	int* pend;
	int* edges;
	int* tri_cnt;
	int* edge_list;
	int* edgelist_pointer;
	ui* oldID;
	int* newID;
	int* rid_vertex;
	int* rid_edges;
	int* sid_vertex;
	int* sid_edges;

	int* deleted_vertex;
	int* deleted_edges;
	int* v_del_edges;
	int* inS_v;
	int* inS_e;
	vector<int> del_v;
	vector<int> del_e;
	vector<int> svex;
	vector<int> sedge;
	vector<int> nodes;

	int* notadj_v;
	int* notadj_e;
	int* adj_v;
	int LB;
	int s_size;


	int* degree;
	int* deleted_flag;

public:
	BnBplex(int _S, int lb)
	{
		pstart = nullptr;
		pend = nullptr;
		v_del_edges = nullptr;
		edges = nullptr;
		tri_cnt = nullptr;
		edge_list = nullptr;
		edgelist_pointer = nullptr;
		oldID = nullptr;
		newID = nullptr;
		rid_vertex = nullptr;
		rid_edges = nullptr;
		sid_vertex = nullptr;
		sid_edges = nullptr;
		deleted_vertex = nullptr;
		deleted_edges = nullptr;
		inS_v = nullptr;
		inS_e = nullptr;
		degree = nullptr;
		S = _S;
		LB = lb;
		s_size = 0;
	}
	~BnBplex()
	{
		if (deleted_flag != nullptr) {
			delete[] deleted_flag;
			deleted_flag = nullptr;
		}
		if (pstart != nullptr) {
			delete[] pstart;
			pstart = nullptr;
		}
		if (pend != nullptr) {
			delete[] pend;
			pend = nullptr;
		}
		if (edgelist_pointer != NULL) {
			delete[] edgelist_pointer;
			edgelist_pointer = NULL;
		}
		if (edge_list != NULL) {
			delete[] edge_list;
			edge_list = NULL;
		}
		if (tri_cnt != NULL) {
			delete[] tri_cnt;
			tri_cnt = NULL;
		}
		if (oldID != nullptr) {
			delete[] oldID;
			oldID = nullptr;
		}
		if (newID != nullptr) {
			delete[] newID;
			newID = nullptr;
		}
		if (edges != nullptr) {
			delete[] edges;
			edges = nullptr;
		}
		if (rid_vertex != nullptr) {
			delete[] rid_vertex;
			rid_vertex = nullptr;
		}
		if (rid_edges != nullptr) {
			delete[] rid_edges;
			rid_edges = nullptr;
		}
		if (sid_vertex != nullptr)
		{
			delete[] sid_vertex;
			sid_vertex = nullptr;
		}
		if (sid_edges != nullptr)
		{
			delete[] sid_edges;
			sid_edges = nullptr;
		}
		if (deleted_vertex != nullptr)
		{
			delete[] deleted_vertex;
			deleted_vertex = nullptr;
		}
		if (deleted_edges != nullptr)
		{
			delete[] deleted_edges;
			deleted_edges = nullptr;
		}
		if (inS_v != nullptr)
		{
			delete[] inS_v;
			inS_v = nullptr;
		}
		if (inS_e != nullptr)
		{
			delete[] inS_e;
			inS_e = nullptr;
		}
		if (notadj_v != nullptr)
		{
			delete[] notadj_v;
			notadj_v = nullptr;
		}
		if (notadj_e != nullptr)
		{
			delete[] notadj_e;
			notadj_e = nullptr;
		}
		if (adj_v != nullptr)
		{
			delete[] adj_v;
			adj_v = nullptr;
		}
		if (degree != nullptr)
		{
			delete[] degree;
			degree = nullptr;
		}

		del_v.clear();
		del_e.clear();
		svex.clear();
		sedge.clear();
	}


	void load_graph(int must_include_, int rid_n, ui* rid, vector< pair<int, int> >& vp, ui* old_id)
	{
		n = rid_n;
		m = vp.size() * 2;
		pstart = new int[n + 1];
		pend = new int[n];
		edges = new int[m];
		tri_cnt = new int[m >> 1];
		edge_list = new int[m];
		edgelist_pointer = new int[m];
		oldID = new ui[n];
		deleted_vertex = new int[n];
		deleted_edges = new int[m >> 1];
		v_del_edges = new int[m >> 1];
		degree = new int[n];
		must_include = must_include_;

		for (int i = 0; i < n; i++)
		{
			oldID[i] = old_id[i];
		}

		notadj_e = new int[m >> 1];
		notadj_v = new int[n];
		adj_v = new int[n];

		inS_v = new int[n];
		inS_e = new int[m >> 1];
		deleted_flag = new int[m >> 1];

		memset(deleted_flag, 0, sizeof(int) * (m >> 1));
		memset(deleted_vertex, 0, sizeof(int) * n);
		memset(deleted_edges, 0, sizeof(int) * (m >> 1));
		memset(v_del_edges, -1, sizeof(int) * (m >> 1));
		memset(tri_cnt, 0, sizeof(int) * (m >> 1));
		memset(degree, 0, sizeof(int) * n);
		memset(notadj_e, 0, sizeof(int) * (m >> 1));
		memset(notadj_v, 0, sizeof(int) * n);
		memset(adj_v, 0, sizeof(int) * n);
		memset(inS_v, 0, sizeof(int) * n);
		memset(inS_e, 0, sizeof(int) * (m >> 1));



		int size = vp.size();
		for (int i = 0; i < size; i++)
		{
			vp.push_back(mp(vp[i].second, vp[i].first));
		}

		sort(vp.begin(), vp.end());

		int eid = 0;
		pstart[0] = 0;
		for (int i = 0; i < rid_n; i++)
		{
			pstart[i + 1] = pstart[i];
			while (eid < vp.size() && vp[eid].first == i) edges[pstart[i + 1]++] = vp[eid++].second;
		}
		printf("Before the ENUM stage, there are %u vertices and %lu edges\n", n, m >> 1);
	}
	void rearrange_graph()
	{
		oriented_triangle_counting(n, m, pstart, pend, edges, edgelist_pointer);
		int* pend2 = new int[n];
		reorganize_oriented_graph(n, tri_cnt, edge_list, pstart, pend, pend2, edges, edgelist_pointer);
		for (int i = 0; i < n; i++) pend[i] = pstart[i + 1];
		for (int i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];

		addtoS(must_include);
	}
	void oriented_triangle_counting(int n, int m, int* pstart, int* pend, int* edges, int* tri_cnt)
	{
		for (int i = 0; i < n; i++)
		{
			int& end = pend[i] = pstart[i];
			for (int j = pstart[i]; j < pstart[i + 1]; j++)
			{
				if (edges[j] > i) edges[end++] = edges[j];
			}
		}


#ifndef  NDEBUG
		long long sum = 0;
		for (int i = 0; i < n; i++) sum += pend[i] - pstart[i];
		assert(sum * 2 == m);
#endif // ! NDEBUG

		int* adj = new int[n];
		int cnt = 0;
		memset(tri_cnt, 0, m * sizeof(int));
		memset(adj, 0, n * sizeof(int));

		for (int i = 0; i < n; i++)
		{
			for (int u = pstart[i]; u < pend[i]; u++) adj[edges[u]] = u + 1;

			for (int u = pstart[i]; u < pend[i]; u++)
			{
				int v = edges[u];
				for (int j = pstart[v]; j < pend[v]; j++)
				{
					if (adj[edges[j]])
					{
						tri_cnt[j]++;
						tri_cnt[u]++;
						tri_cnt[adj[edges[j]] - 1]++;
						cnt++;
					}
				}
			}

			for (int u = pstart[i]; u < pend[i]; u++) adj[edges[u]] = 0;
		}

#ifndef  NDEBUG

		//for (int i = 0; i < n; i++)
		//{
		//	for (int u = pstart[i]; u < pend[i]; u++)
		//		cout << tri_cnt[u] << endl;


		//}

		printf("Total number of triangles:%d\n", cnt);

#endif // ! NDEBUG
	}
	void reorganize_oriented_graph(int n, int* tri_cnt, int* edge_list, int* pstart, int* pend, int* pend2, int* edges, int* edgelist_pointer)
	{
		memset(tri_cnt, 0, (m / 2) * sizeof(int));
		for (int i = 0; i < n; i++) pend2[i] = pend[i];


		int pos = 0;

		for (int i = 0; i < n; i++)
		{
			for (int j = pstart[i]; j < pend[i]; j++)
			{
				tri_cnt[pos >> 1] = edgelist_pointer[j];
				edge_list[pos++] = i, edge_list[pos++] = edges[j];

				int& k = pend2[edges[j]];
				edgelist_pointer[j] = edgelist_pointer[k] = (pos >> 1) - 1;
				edges[k++] = i;
			}
		}



#ifndef NDEBUG
		for (int i = 0; i < n; i++) assert(pend2[i] == pstart[i + 1]);
#endif // !NDEBUG


		for (int i = 0; i < n; i++)
		{
			pend2[i] = pend[i];
			pend[i] = pstart[i];
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = pend2[i]; j < pstart[i + 1]; j++)
			{
				int& k = pend[edges[j]];
				edgelist_pointer[k] = edgelist_pointer[j];
				edges[k++] = i;
			}
		}

		int* ids = new int[m];
		int* buf = new int[m];

		for (int i = 0; i < n; i++)
		{
			if (pend[i] == pstart[i] || pend[i] == pstart[i + 1]) continue;
			int k = pend[i], j = pstart[i], pos = 0;
			while (j < pend[i] && k < pstart[i + 1])
			{
				if (edges[j] < edges[k])
				{
					ids[pos] = edges[j];
					buf[pos++] = edgelist_pointer[j++];
				}
				else
				{
					ids[pos] = edges[k];
					buf[pos++] = edgelist_pointer[k++];
				}
			}

			while (j < pend[i])
			{
				ids[pos] = edges[j];
				buf[pos++] = edgelist_pointer[j++];
			}

			while (k < pstart[i + 1])
			{
				ids[pos] = edges[k];
				buf[pos++] = edgelist_pointer[k++];
			}

			for (int j = 0; j < pos; j++)
			{
				edges[pstart[i] + j] = ids[j];
				edgelist_pointer[pstart[i] + j] = buf[j];
			}
		}

	}

	void cal_notadj(int e)
	{
		notadj_e[e] = 0;
		for (auto v : svex)
		{
			if (!check_common_neighbour(v, e))
				notadj_e[e]++;
		}
	}


	void delfrD_vertex(int line, int u)
	{
		if (deleted_vertex[u])
		{
			printf("%d has been deleted\n", u);
			printf("������õ������: %d\n", line);
		}
		//if (n == 31) printf("del %d\n", u);
		assert(!deleted_vertex[u]);
		deleted_vertex[u] = 1;
		del_v.push_back(u);

		int* exist = new int[n];
		memset(exist, 0, sizeof(int) * n);
		for (int j = pstart[u]; j < pend[u]; j++)
		{
			int v = edges[j];
			//if (!deleted_vertex[v] && !deleted_edges[e])����Ҫ
			exist[v] = 1;
		}

		for (int j = pstart[u]; j < pend[u]; j++)
		{
			int v = edges[j];
			if (!deleted_vertex[v] && !deleted_edges[edgelist_pointer[j]])
			{
				for (int i = pstart[v]; i < pend[v]; i++)
				{
					int w = edges[i], e = edgelist_pointer[i];
					if (!deleted_vertex[w] && !deleted_edges[e] && w > v && exist[w] == 1)
						tri_cnt[e]--;
				}
				deleted_edges[edgelist_pointer[j]] = 2;
				del_e.push_back(edgelist_pointer[j]);
				assert(degree[v] > 0);
				degree[v]--;
			}
		}
		delete[] exist;
	}

	void addtoD_vertex(int line, int u)
	{
		//if (!deleted_vertex[u]) return;
		if (!deleted_vertex[u])
		{
			printf("������õ������: %d\n", line);
		}
		assert(deleted_vertex[u]);
		deleted_vertex[u] = 0;
		del_v.pop_back();

		int* exist = new int[n];
		memset(exist, 0, sizeof(int) * n);
		for (int j = pstart[u]; j < pend[u]; j++)
		{
			int v = edges[j];
			int e = edgelist_pointer[j];
			//if (!deleted_vertex[v] && !deleted_edges[e])//����ȷ��e�ǲ�����Ϊɾ��v��ɾ����
			exist[v] = 1;
		}

		for (int j = pend[u] - 1; j >= (int)pstart[u]; j--)
		{
			int v = edges[j];
			if (deleted_edges[edgelist_pointer[j]] == 2 && !deleted_vertex[v])//
			{
				for (int i = pend[v] - 1; i >= (int)pstart[v]; i--)
				{
					int w = edges[i], e = edgelist_pointer[i];
					if (!deleted_vertex[w] && !deleted_edges[e] && w > v && exist[w] == 1)
						tri_cnt[e]++;
				}
				deleted_edges[edgelist_pointer[j]] = 0;
				del_e.pop_back();
				degree[v]++;

				//cal_notadj(edgelist_pointer[j]);
			}
		}
		delete[] exist;
	}

	void delfrD_edge(int line, int e)
	{
		//printf("try to del edge %d-%d\n", edge_list[e << 1], edge_list[(e << 1) + 1]);

		if (deleted_edges[e])
		{
			printf("%d-%d has been deleted\n", edge_list[e << 1], edge_list[(e << 1) + 1]);
			printf("������õ������: %d\n", line);
		}
		assert(!deleted_edges[e]);

		deleted_edges[e] = 1;
		del_e.push_back(e);
		int u = edge_list[e << 1], v = edge_list[(e << 1) + 1];
		if (degree[u] <= 0 || degree[v] <= 0) printf("degree u = %d and degree v = %d\n", degree[u], degree[v]);
		assert(degree[u] > 0);
		assert(degree[v] > 0);
		degree[u]--, degree[v]--;

		int ii = pstart[u], jj = pstart[v];

		while (true)
		{
			while (ii < pend[u] && deleted_edges[edgelist_pointer[ii]]) ii++;
			while (jj < pend[v] && deleted_edges[edgelist_pointer[jj]]) jj++;

			if (ii >= pend[u] || jj >= pend[v]) break;
			if (edges[ii] == edges[jj])
			{
				tri_cnt[edgelist_pointer[ii]]--;
				tri_cnt[edgelist_pointer[jj]]--;

				ii++;
				jj++;
			}
			else if (edges[ii] < edges[jj]) ii++;
			else if (edges[ii] > edges[jj]) jj++;
		}
	}


	void addtoD_edge(int e)
	{
		assert(deleted_edges[e]);
		deleted_edges[e] = 0;
		del_e.pop_back();
		int u = edge_list[e << 1], v = edge_list[(e << 1) + 1];
		degree[u]++, degree[v]++;

		int ii = pstart[u], jj = pstart[v];

		while (true)
		{
			while (ii < pend[u] && deleted_edges[edgelist_pointer[ii]]) ii++;
			while (jj < pend[v] && deleted_edges[edgelist_pointer[jj]]) jj++;

			if (ii >= pend[u] || jj >= pend[v]) break;
			if (edges[ii] == edges[jj])
			{
				tri_cnt[edgelist_pointer[ii]]++;
				tri_cnt[edgelist_pointer[jj]]++;

				ii++;
				jj++;
			}
			else if (edges[ii] < edges[jj]) ii++;
			else if (edges[ii] > edges[jj]) jj++;
		}

		cal_notadj(e);
	}

	bool check_common_neighbour(int v, int e)
	{
		/*assert(!deleted_edges[e]);*/
		//printf("%d\n", e);
		int flag = 0;
		int u = edge_list[e << 1], w = edge_list[(e << 1) + 1];
		if (u == v || w == v) return true;//我靠，这句话太重要了
		//printf("%d\n", e);
		for (int i = pstart[u]; i < pend[u]; i++)
		{
			if (edges[i] == v && !deleted_edges[edgelist_pointer[i]])
			{
				flag = 1;
				break;
			}
		}

		//printf("%d %d\n", u, w);
		for (int i = pstart[w]; i < pend[w]; i++)
		{

			if (edges[i] == v && flag == 1 && !deleted_edges[edgelist_pointer[i]])
			{
				flag = 2;
				break;
			}
		}

		if (flag == 2) return true;
		return false;
	}

	bool check_notanyone_inS(int v)
	{
		if (svex.empty()) return false;
		for (int i = pstart[v]; i < pend[v]; i++)
		{
			int v = edges[i];
			int e = edgelist_pointer[i];
			if (inS_v[v] && !deleted_edges[e]) return false;
		}
		return true;
	}



	bool canadd(int v)
	{
		if (svex.empty()) return true;
		

		if (cal_notadj_v(v) >= S + 1) return false;

		if (check_notanyone_inS(v))
		{
			if (svex.size() >= S + 1) return false;
			else return true;
		}



		bool ret = true;
		int* hash = new int[m];
		memset(hash, 0, sizeof(int) * m);

		vector<int> sedge_temp = sedge; 


		svex.push_back(v);
		inS_v[v] = true;

		for (int k = pstart[v]; k < pend[v]; k++)
		{
			int u = edges[k];
			int e = edgelist_pointer[k];
			if (!deleted_edges[e] && inS_v[u])
				sedge_temp.push_back(e);
		}



		vector<int> todel_e;
		queue<int> Q_del;




		for (auto it = sedge_temp.begin(); it != sedge_temp.end();)
		{
			int e = *it;
			cal_notadj(e);
			if (notadj_e[e] > S)
			{
				Q_del.push(e);

				it = sedge_temp.erase(it);
			}
			else ++it;
		}


		while (!Q_del.empty())
		{
			int e = Q_del.front();
			Q_del.pop();
			todel_e.push_back(e);
			assert(hash[e] == 0);
			hash[e] = 1;
			//printf("del %d\n", e);
			delfrD_edge(__LINE__, e);
			for (auto it = sedge_temp.begin(); it != sedge_temp.end();)
			{
				e = *it;
				cal_notadj(e);
				if (notadj_e[e] > S)
				{
					Q_del.push(e);

					it = sedge_temp.erase(it);
				}
				else ++it;
			}
		}



		for (auto u : svex)
		{
			if (check_notanyone_inS(u))
				ret = false;
		}

		for (int e = todel_e.size() - 1; e >= 0; e--)
			addtoD_edge(todel_e[e]);


		svex.pop_back();
		inS_v[v] = false;

		free(hash);

		return ret;
	}


	void addtoS(int u)
	{
		inS_v[u] = true;
		svex.push_back(u);
		s_size++;

		for (int k = pstart[u]; k < pend[u]; k++)
		{
			int v = edges[k];
			int e = edgelist_pointer[k];
			if (!deleted_edges[e] && inS_v[v])
			{
				sedge.push_back(e);
				inS_e[e] = true;
			}
		}

		queue<int> Q_del;

		for (auto it = sedge.begin(); it != sedge.end();)
		{
			int e = *it;
			cal_notadj(e);
			if (notadj_e[e] > S)
			{
				Q_del.push(e);
				it = sedge.erase(it);
			}
			else ++it;
		}

		while (!Q_del.empty())
		{
			int e = Q_del.front();
			Q_del.pop();
			delfrD_edge(__LINE__, e);

			v_del_edges[e] = u;
			for (auto it = sedge.begin(); it != sedge.end();)
			{
				e = *it;
				cal_notadj(e);
				if (notadj_e[e] > S)
				{
					Q_del.push(e);
					it = sedge.erase(it);
				}
				else ++it;
			}
		}

	}

	void delfrS(int u)
	{
		inS_v[u] = false;
		svex.pop_back();
		s_size--;


		for (int e = 0; e < (m >> 1); e++)
		{

			if (deleted_edges[e] && v_del_edges[e] == u)
			{
				addtoD_edge(e);
				int x = edge_list[e << 1], y = edge_list[(e << 1) + 1];

				if (inS_v[x] && inS_v[y])
				{
					sedge.push_back(e);
					inS_e[e] = true;
				}
			}
		}

		for (int k = pstart[u]; k < pend[u]; k++)
		{
			int v = edges[k];
			int e = edgelist_pointer[k];
			if (!deleted_edges[e] && inS_v[v])
			{
				inS_e[e] = false;
				sedge.erase(remove(sedge.begin(), sedge.end(), e), sedge.end());
			}
		}
	}


	bool check_adj_vertex(int u, int v)
	{
		for (int i = pstart[u]; i < pend[u]; i++)
		{
			int e = edgelist_pointer[i];
			if (!deleted_vertex[edges[i]] && edges[i] == v) return true;//!deleted_edges[e]&& 
		}
		return false;
	}

	int cal_notadj_v(int v)
	{
		int cnt = 0;
		for (auto u : svex)
			if (!check_adj_vertex(u, v) && u != v) cnt++;
		return cnt;
	}

	int combination(int n, int k) {
		vector<vector<int>> dp(n + 1, vector<int>(k + 1, 0));

		for (int i = 0; i <= n; ++i) {
			for (int j = 0; j <= min(i, k); ++j) {
				if (j == 0 || j == i)
					dp[i][j] = 1;
				else
					dp[i][j] = dp[i - 1][j - 1] + dp[i - 1][j];
			}
		}

		return dp[n][k];
	}



	void findSubsets(vector<vector<int>>& result, vector<int>& subset, vector<int> remain_undel_v, int start,int n, int s, int target) 
	{
		if (result.size() >= target) return;

		if (s == 0) {
			result.push_back(subset);
			printf("Remain subset:%d\n", target - result.size());
			return;
		}

		for (int i = start; i < remain_undel_v.size(); ++i) {
			int v = remain_undel_v[i];
			subset.push_back(v);
			findSubsets(result, subset, remain_undel_v, i + 1, n, s - 1, target);
			subset.pop_back();
		}
	}


	bool bb_s_plex(int curS, clock_t start, bool& timeflag)
	{
		for (auto u : svex)
			inS_v[u] = 100;
		for (int i = 0; i < n; i++)
		{
			if (inS_v[i] != 100)
				inS_v[i] = 0;
		}
		for (auto u : svex)
			inS_v[u] = 1;


		//pruning rules

		if (svex.size() >= LB)
		{
			cout << "规模大于LB，返回正确，时间是：" << endl;
			return true;
		}

		clock_t finish;
		finish = clock(); 
		double time = (double)(finish - start) / CLOCKS_PER_SEC;
		if (time > 86400)
		{
			timeflag = false;
			return false;
		}

		if (curS < LB) return false;

		int minID = -1;
		for (int i = 0; i < n; ++i) if (!deleted_vertex[i]) {
			if (minID == -1 || degree[i] < degree[minID])
				minID = i;
		}

		int min_tri = n, min_eid = -1;
		for (int i = 0; i < (m >> 1); i++)
		{
			if (!deleted_edges[i] && tri_cnt[i] < min_tri)
			{
				min_tri = tri_cnt[i];
				min_eid = i;
			}
		}


		if (tri_cnt[min_eid] >= curS - 2 - S && degree[minID] >= curS - S - 1)
		{
			//printf("min tri_cnt = %d min_e = %d-%d\n", tri_cnt[min_eid], oldID[edge_list[min_eid << 1]], oldID[edge_list[(min_eid << 1) + 1]]);
			//output_min_k_truss();
			printf("The size of current vertex set is %d\n\n", svex.size());
			printf("There are %d vertices remaining\n\n", curS);
			return true;
		}
		

		if (degree[minID] <= LB - 2 - S) //说明MinID不可能加入最终解，第9行
		{
			if (inS_v[minID])
				return false;
			delfrD_vertex(__LINE__, minID);
			bool ret = bb_s_plex(curS - 1, start, timeflag);
			addtoD_vertex(__LINE__, minID);
			return ret;
		}

		int maxID = -1;//maxID指在s中非邻居数量最多的节点
		for (int i = 0; i < n; ++i) if (!deleted_vertex[i]) {
			if (maxID == -1 || cal_notadj_v(i) > cal_notadj_v(maxID))
				maxID = i;
		}

		//如果maxID不在s中，且他在s中非邻居大于S+1，那他必被删，第5行，lemma2

		if (!inS_v[maxID] && cal_notadj_v(maxID) >= S + 1) {
			delfrD_vertex(__LINE__, maxID);
			bool ret = bb_s_plex(curS - 1, start, timeflag);
			addtoD_vertex(__LINE__, maxID);
			return ret;
		}

		//如果maxID在s中，且非邻居大于S(不包括他自己，所以加上她自己就是S+1），该节点是临界点，可以借助临界点来删除c中和临界点不相邻的节点，第5行，lemma2

		if (inS_v[maxID] && cal_notadj_v(maxID) >= S) {
			if (cal_notadj_v(maxID) >= S + 1) return false;
			std::vector <int> todel;
			for (int i = 0; i < n; ++i) if (!deleted_vertex[i] && !inS_v[i] && !check_adj_vertex(maxID, i)) {
				todel.push_back(i);
			}
			if (todel.size()) {
				for (auto x : todel) delfrD_vertex(__LINE__,x);
				bool ret = bb_s_plex(curS - todel.size(), start, timeflag);
				for (auto x : todel) addtoD_vertex(__LINE__, x);
				return ret;
			}
		}


		//开始分支（以minID为基准，在这里指的是k-unsatisfied节点，因为若整个图不满足780行的情况，那么minID肯定是k-unsatisfied)
		std::vector <int> branch;
		for (int x = 0; x < n; ++x) if (!deleted_vertex[x] && x != minID && !inS_v[x]) {
			if (!check_adj_vertex(minID, x))
				branch.push_back(x);
		}

		random_shuffle(branch.begin(), branch.end());

		//若minID在s中
		if (inS_v[minID]) {
			int canselect = S - cal_notadj_v(minID), pos = -1;
			bool ret = false;
			for (int i = 0; !ret && i < canselect; ++i) {
				delfrD_vertex(__LINE__, branch[i]);
				if (i && !canadd(branch[i - 1])) {
					addtoD_vertex(__LINE__, branch[i]);
					break;
				}
				if (i) {
					addtoS(branch[i - 1]);
					pos = i - 1;
				}
				ret |= bb_s_plex(curS - 1, start, timeflag);
				addtoD_vertex(__LINE__, branch[i]);
			}
			if (ret) {
				for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
				return true;
			}
			//后续与minID不相连的节点都得删除
			for (int i = canselect; i < (int)branch.size(); ++i) {
				delfrD_vertex(__LINE__, branch[i]);
			}
			//最后一个分支
			if (canselect == 0 || canadd(branch[canselect - 1])) {
				if (canselect) addtoS(branch[canselect - 1]);
				ret |= bb_s_plex(curS - branch.size() + canselect, start, timeflag);
				if (canselect) delfrS(branch[canselect - 1]);
			}
			for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
			for (int i = canselect; i < (int)branch.size(); ++i) {
				addtoD_vertex(__LINE__, branch[i]);
			}
			return ret;
		}

		//若minID不在s中
		else {
			delfrD_vertex(__LINE__, minID);

			bool ret = bb_s_plex(curS - 1, start, timeflag);
			addtoD_vertex(__LINE__, minID);

			if (ret) return true;
			int canselect = S - cal_notadj_v(minID);
			if (!canadd(minID)) return false;
			addtoS(minID);
			int pos = -1;
			for (int i = 0; !ret && i < canselect; ++i) {
				delfrD_vertex(__LINE__, branch[i]);
				if (i && !canadd(branch[i - 1])) {
					addtoD_vertex(__LINE__, branch[i]);
					break;
				}
				if (i) {
					addtoS(branch[i - 1]);
					pos = i - 1;
				}
				ret |= bb_s_plex(curS - 1, start,  timeflag);
				addtoD_vertex(__LINE__, branch[i]);
			}
			if (ret) {
				for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
				delfrS(minID);
				return true;
			}
			for (int i = canselect; i < (int)branch.size(); ++i) {
				delfrD_vertex(__LINE__, branch[i]);
			}
			if (canselect == 0 || canadd(branch[canselect - 1])) {
				if (canselect) addtoS(branch[canselect - 1]);
				ret |= bb_s_plex(curS - branch.size() + canselect, start, timeflag);
				if (canselect) delfrS(branch[canselect - 1]);
			}
			for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
			for (int i = canselect; i < (int)branch.size(); ++i) {
				addtoD_vertex(__LINE__, branch[i]);
			}
			delfrS(minID);
			return ret;
		}

		return false;
	}

	void output_min_k_truss()
	{
		vector<int> ans_v;
		for (int i = 0; i < n; i++)
			if (!deleted_vertex[i]) ans_v.push_back(i);
		ofstream outfile;
		string outname = "/home/zhangqifan/min_k_truss/dataset/min_k_truss_temp" + to_string(oldID[must_include]) + ".txt";
		outfile.open(outname);
		set<pair<int, int>> ans_e;
		int cnt = 0;
		for (auto v : ans_v)
		{
			cnt++;
			printf("ans %d = %d ", cnt, oldID[v]);
			int deg = 0;
			for (int j = pstart[v]; j < pend[v]; j++)
			{
				int u = edges[j];
				if (!deleted_vertex[u])deg++;
			}
			printf("degree = %d\n", deg);



			for (int j = pstart[v]; j < pend[v]; j++)
			{
				int u = edges[j];
				if (!deleted_vertex[u])
				{
					if (v < u)
					{
						ans_e.insert(mp(oldID[v], oldID[u]));
					}
					else
					{
						ans_e.insert(mp(oldID[u], oldID[v]));
					}
				}
			}
		}
		for (auto e : ans_e)
			outfile << e.first << " " << e.second << " " << endl;
		outfile.close();
	}



};








#endif // !_BNB_H_
