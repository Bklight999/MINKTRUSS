#ifndef _KTB_ENUM_TOPC_H
#define _KTB_ENUM_TOPC_H

#include"Utility.h"
#include"Timer.h"

class ENUMC
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
	vector<vector<int>> top_ans;

	int* notadj_v;
	int* notadj_e;
	int* adj_v;
	int LB;
	int s_size;
	int c;


	int* degree;
	int* deleted_flag;

public:
	ENUMC(int _S, int lb)
	{
		pstart = nullptr;
		pend = nullptr;
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
	~ENUMC()
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
	void load_graph_IE(int must_include_, int rid_n, ui* rid, vector< pair<int, int> >& vp, ui* old_id, int c, vector<vector<int>> topc_ans)
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
		this->c = c;

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

	void load_graph_GLOBAL(int rid_n, ui* rid, vector< pair<int, int> >& vp)
	{
		must_include = -1;
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
		degree = new int[n];

		notadj_e = new int[m >> 1];
		notadj_v = new int[n];
		adj_v = new int[n];

		inS_v = new int[n];
		inS_e = new int[m >> 1];
		deleted_flag = new int[m >> 1];

		memset(deleted_flag, 0, sizeof(int) * (m >> 1));
		memset(deleted_vertex, 0, sizeof(int) * n);
		memset(deleted_edges, 0, sizeof(int) * (m >> 1));
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

	bool check_adj_vertex(int u, int v)
	{
		if (u == v) return false;
		for (int i = pstart[u]; i < pend[u]; i++)
			if (!deleted_vertex[edges[i]] && !deleted_edges[edgelist_pointer[i]] && edges[i] == v) return true;
		return false;
	}

	int cal_notadj_v(int v)
	{
		int cnt = 0;
		for (auto u : svex)
			if (!check_adj_vertex(u, v)) cnt++;
		return cnt;
	}

	bool check_common_neighbour(int v, int e)
	{
		//printf("%d\n", e);
		int flag = 0;
		int u = edge_list[e << 1], w = edge_list[(e << 1) + 1];
		if (u == v || w == v) return true;
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

	bool check_neighbour_inS(int u, int v)
	{
		for (int ii = pstart[u]; ii < pend[u]; ii++)
		{
			int x = edges[ii], e = edgelist_pointer[ii];
			if (!deleted_edges[e] && x == v) return true;
		}

		return false;
	}

	void delfrD_vertex(int line, int u)
	{
		if (deleted_vertex[u])
		{
			printf("%d has been deleted\n", u);
			printf("错误调用的入口在: %d\n", line);
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
			//if (!deleted_vertex[v] && !deleted_edges[e])不需要
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
			printf("错误调用的入口在: %d\n", line);
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
			//if (!deleted_vertex[v] && !deleted_edges[e])//不能确定e是不是因为删除v才删除的
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

	void cal_notadj(int e)
	{
		notadj_e[e] = 0;
		for (auto v : svex)
		{
			if (!check_common_neighbour(v, e))
				notadj_e[e]++;
		}
	}

	void delfrD_edge(int line, int e)
	{
		//printf("try to del edge %d-%d\n", edge_list[e << 1], edge_list[(e << 1) + 1]);

		if (deleted_edges[e])
		{
			printf("%d-%d has been deleted\n", edge_list[e << 1], edge_list[(e << 1) + 1]);
			printf("错误调用的入口在: %d\n", line);
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
			//printf("del %d-%d\n", edge_list[e << 1], edge_list[(e << 1) + 1]);
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
		//printf("stage check can add %d!!!!!!!!!\n\n", v);
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

	void output_svex()
	{
		string outname = "/home/zhangqifan/min_k_truss/dataset/min_k_truss_temp" + to_string(oldID[must_include]) + ".txt";
		ofstream outfile;
		outfile.open(outname);
		set<pair<int, int>> ans_e;
		int cnt = 0;
		vector<int> ans;
		for (auto u : svex)
		{
			if (cnt == LB) break;
			if (!check_notanyone_inS(u))
			{
				ans.push_back(u);
				cnt++;
			}
		}

		for (int i = 0; i < ans.size(); i++)
		{
			for (int j = i + 1; j < ans.size(); j++)
			{
				int u = ans[i], v = ans[j];
				if (check_adj_vertex(u, v))
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

	void backtrack(vector<int>& nums, vector<int>& cur, int r, int start, int target_num)
	{
		if (top_ans.size() >= target_num) return;

		if (cur.size() == r)
		{
			for (int i = 0; i < top_ans.size(); i++)
			{
				if (cur == top_ans[i])
					return;
			}
			top_ans.push_back(cur);
			/*for (auto x : cur)
				cout << x << " ";
			cout << endl;*/
			return;
		}

		for (int i = start; i < nums.size(); i++)
		{
			cur.push_back(nums[i]);
			backtrack(nums, cur, r, i + 1, target_num);
			cur.pop_back();
		}
	}


	bool find_bounded_stb(clock_t start, bool& timeflag, int &C)
	{
		oriented_triangle_counting(n, m, pstart, pend, edges, edgelist_pointer);
		int* pend2 = new int[n];
		reorganize_oriented_graph(n, tri_cnt, edge_list, pstart, pend, pend2, edges, edgelist_pointer);

		if (must_include != -1) addtoS(must_include); 


		for (int i = 0; i < n; i++) pend[i] = pstart[i + 1];
		for (int i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];


		bool ret;
		if (must_include != -1)
		{
			printf("must_include = %d\n", must_include);
			nodes.resize(n - 1);

			for (int i = 1; i < n; i++)
				nodes[i - 1] = i;
			
			bool flag = false;
			ret = enumeration(flag, start, timeflag, C);
		}
		else
		{
			nodes.resize(n);
			for (int i = 0; i < n; i++)
				nodes[i] = i;
			std::mt19937 gen(12345);  
			std::shuffle(nodes.begin(), nodes.end(), gen); 
			bool flag = false;

			ret = enumeration(flag, start, timeflag, C);
		}

		return ret;
	}

	bool canadd_old(int v)
	{
		if (svex.empty()) return true;
		int* flag = new int[m >> 1];
		memset(flag, 0, sizeof(int) * (m >> 1));

		for (int k = pstart[v]; k < pend[v]; k++)
		{
			ui u = edges[k], e1 = edgelist_pointer[k];
			cal_notadj(e1);

			if (!deleted_edges[e1] && inS_v[u] && notadj_e[e1] > S)
			{
				delfrD_edge(__LINE__, e1);
				flag[e1] = 1;
			}
		}
		//int tot = 0;
		for (auto e : sedge)
		{
			if (!check_common_neighbour(v, e))
			{
				cal_notadj(e);
				if (notadj_e[e] > S) printf("S = %d, e = %d - %d, notadj_e[e] = %d\n", S, edge_list[e << 1], edge_list[(e << 1) + 1], notadj_e[e]);
				assert(notadj_e[e] <= S);
				//if ((++tot) > S) return false;
				if (notadj_e[e] == S)
				{
					for (int k = pend[v] - 1; k >= pstart[v]; k--)
					{
						ui u = edges[k], e1 = edgelist_pointer[k];
						if (flag[e1] == 1)
						{
							addtoD_edge(e1);
							flag[e1] = 0;
						}
					}
					return false;
				}
			}
		}

		if (check_notanyone_inS(v))
		{
			for (int k = pend[v] - 1; k >= pstart[v]; k--)
			{
				ui u = edges[k], e1 = edgelist_pointer[k];
				if (flag[e1] == 1)
				{
					addtoD_edge(e1);
					flag[e1] = 0;
				}
			}
			return false;
		}

		for (int k = pend[v] - 1; k >= pstart[v]; k--)
		{
			ui u = edges[k], e1 = edgelist_pointer[k];
			if (flag[e1] == 1)
			{
				addtoD_edge(e1);
				flag[e1] = 0;
			}
		}
		delete[] flag;
		return true;
	}

	void addtoS_old(int u)
	{

		for (int k = pstart[u]; k < pend[u]; k++)
		{
			int v = edges[k];
			int e = edgelist_pointer[k];
			cal_notadj(e);

			if (!deleted_edges[e] && inS_v[v] && notadj_e[e] > S)
			{
				delfrD_edge(__LINE__, e);
				deleted_edges[e] = 10;
			}

		}

		for (int k = pstart[u]; k < pend[u]; k++)
		{
			int v = edges[k];
			int e = edgelist_pointer[k];
			cal_notadj(e);
			if (!deleted_edges[e] && inS_v[v] && notadj_e[e] <= S)
			{
				inS_e[e] = true;
				sedge.push_back(e);
			}
		}

	

		inS_v[u] = true;
		svex.push_back(u);
		s_size++;



	}

	void delfrS_old(int u)
	{
		inS_v[u] = false;
		svex.pop_back();
		s_size--;

		for (int k = pend[u] - 1; k >= pstart[u]; k--)
		{
			int e = edgelist_pointer[k];
			int v = edges[k];
			if (inS_e[e])
			{
				inS_e[e] = false;
				sedge.pop_back();
			}
			else if (deleted_edges[e] == 10 && inS_v[v]) 
				addtoD_edge(e);
		}


	}

private:


	bool enumeration(bool& flag, clock_t start, bool& timeflag, int& C) 
	{
		if (flag) return true;

		clock_t finish;
		finish = clock();
		double time = (double)(finish - start) / CLOCKS_PER_SEC;
		if (time > 86400)
		{
			timeflag = false;
			return false;
		}


		bool ret = false, empty = false;


		if (svex.size() >= LB)
		{
			int count = 0;
			for (auto u : svex)
			{
				if (check_notanyone_inS(u))
				{
					count++;
				}
			}

			if (svex.size() - count < LB)
			{
				empty = true;
			}

			if (!empty)
			{
				printf("\n\n\noptimal solution = %d\n\n\n", svex.size());
				printf("time = %lf\n", time);
				//exit(1000);
				printf("\n");
				flag = true;
				return true;
			}
		}

		

		int curS = nodes.size() + svex.size();

		int curSS = 0;
		for (int i = 0; i < n; i++)
		{
			if (!deleted_vertex[i])
				curSS++;
		}
		assert(curS == curSS);
		

		if (curS < LB)
		{
			//printf("curS = %d < LB = %d\n", curS, LB);
			return false;
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

		if (tri_cnt[min_eid] >= curS - 2 - S)
		{
			printf("The size of current vertex set is %d\n\n", svex.size());
			printf("There are %d vertices remaining\n\n", curS);
			vector<int> nums, cur;
			for (int i = 0; i < n; i++)
			{
				if (!deleted_vertex[i] && i != must_include)
					nums.push_back(oldID[i]);
			}
			/*	for (auto x : nums)
					cout << x << " ";
				cout << endl;*/

			cur.push_back(oldID[must_include]);
			//printf("CurS = %d, nums.size() = %d\n", curS, nums.size());
			int x = curS - 1, y = LB - 1;
			printf("Extract %d from %d, current size = %d\n", y, x, top_ans.size());
			backtrack(nums, cur, LB, 0, c);
			C = top_ans.size();
			finish = clock();
			double temp_time = (double)(finish - start) / CLOCKS_PER_SEC;
			printf("Already find %d ans! Remain %d\n", C, (c - C) > 0 ? c - C : 0);
			printf("Time cost: %lf\n\n", temp_time);
			if (C >= c)return true;
			return false;
		}

		vector<int> del;
		for (auto i : nodes)
		{
			if (!deleted_vertex[i] && cal_notadj_v(i) >= S + 1) //!inS_v[i] && 
				del.push_back(i);
		}

	

		if (del.size()) {

			for (int x = 0; x < del.size(); x++)
			{
				delfrD_vertex(__LINE__, del[x]);
				nodes.erase(remove(nodes.begin(), nodes.end(), del[x]), nodes.end());
			}
			bool ret1 = enumeration(flag, start, timeflag, C);
			//printf("%d\n", curS - todel.size());
			for (int x = del.size() - 1; x >= 0; x--)
			{
				addtoD_vertex(__LINE__, del[x]);
				nodes.push_back(del[x]);
			}
			return ret1;
		}



		
		int max_v = -1, max_ins_degree = m;
		for (auto x : nodes)
		{
			//if (deleted_vertex[x]) continue;
			int cnt = 0;
	
			for (auto u : svex)
			{
				if (check_neighbour_inS(u, x))
					cnt++;
			}
			if (cnt < max_ins_degree)
			{
				max_ins_degree = cnt;
				max_v = x;
			}
			//printf("v = %d and degree = %d\n", x, cnt);
		}
		if (max_v == -1) return false;

		nodes.erase(remove(nodes.begin(), nodes.end(), max_v), nodes.end());


		int current_node = max_v;//nodes[v+1]
		if (canadd(current_node))
		{
			addtoS(current_node);
			ret |= enumeration(flag, start, timeflag, C);
			delfrS(current_node);
		}

		delfrD_vertex(__LINE__, current_node);
	
		ret |= enumeration(flag, start, timeflag, C);
		addtoD_vertex(__LINE__, current_node);

		nodes.push_back(max_v);

		return ret;
	}


};



#endif
