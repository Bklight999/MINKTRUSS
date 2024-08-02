#ifndef _BNB_ITERALL_H_
#define _BNB_ITERALL_H_

#include"Utility.h"
#include"Timer.h"

class BnBiterall
{
private:
	int n;
	int m;
	int S;
	int start_from_zero;
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
	int* inS_v;
	int* inS_e;
	vector<int> del_v;
	vector<int> del_e;
	vector<int> svex;
	vector<int> sedge;

	int* notadj_v;
	int* notadj_e;
	int* adj_v;
	int LB;
	int s_size;
	int must_include;

	int* degree;
	int* deleted_flag;

public:
	BnBiterall(int _S, int lb, int startfromzero)
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
		start_from_zero = startfromzero;
		S = _S;
		LB = lb;
		s_size = 0;
	}
	~BnBiterall()
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

	void load_graph(int must_include_, int rid_n, ui* rid, vector< pair<int, int> >& vp, ui* ids)
	{
		must_include = must_include_;
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

		//printf("exchange\n");


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
			//vp[i].first = newID[vp[i].first];
			//vp[i].second = newID[vp[i].second];
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


		printf("Before the BnB stage, there are %u vertices and %lu edges\n", n, m >> 1);
	}
	void rearrange_graph()
	{
		oriented_triangle_counting(n, m, pstart, pend, edges, edgelist_pointer);
		int* pend2 = new int[n];
		reorganize_oriented_graph(n, tri_cnt, edge_list, pstart, pend, pend2, edges, edgelist_pointer);
		for (int i = 0; i < n; i++) pend[i] = pstart[i + 1];
		for (int i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];

		addtoS(must_include);

		//for (int i = 0; i < n; i++)
		//{
		//	for (int j = pstart[i]; j < pend[i]; j++)
		//	{
		//		printf("%d-%d\n", i, edges[j]);
		//	}
		//}

		//for (int e = 0; e < (m >> 1); e++)
		//	printf("%d-%d\n", edge_list[e << 1], edge_list[(e << 1) + 1]);
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
	}

	//限定只能删除S外的边
	void delfrD_edge(int e)
	{
		if (deleted_edges[e]) printf("%d-%d has been deleted\n", edge_list[e << 1], edge_list[(e << 1) + 1]);
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

				//改变相关联边在S中的notadj_e

				/*int score1 = inS_v[u] + inS_v[v];
				if (score1 == 1)
				{
					int e = inS_v[u] ? edgelist_pointer[jj] : edgelist_pointer[ii];
					notadj_e[e]++;
					deleted_flag[e] = 1;
				}*/

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

				/*		int score1 = inS_v[u] + inS_v[v];

						int e = edgelist_pointer[ii];
						if (deleted_flag[e])
						{
							deleted_flag[e] = 0;
							notadj_e[e]--;
						}
						e = edgelist_pointer[jj];
						if (deleted_flag[e])
						{
							deleted_flag[e] = 0;
							notadj_e[e]--;
						}*/
						//if (score1 == 1)
						//{
						//	//改变相关联边在S中的notadj_e
						//	int e = inS_v[u] ? edgelist_pointer[jj] : edgelist_pointer[ii];
						//	notadj_e[e]--;
						//}

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



		//printf("stage check can add %d\n\n", v);
		int* flag = new int[m >> 1];
		memset(flag, 0, sizeof(int) * (m >> 1));

		for (int k = pstart[v]; k < pend[v]; k++)
		{
			ui u = edges[k], e1 = edgelist_pointer[k];
			cal_notadj(e1);

			if (!deleted_edges[e1] && inS_v[u] && notadj_e[e1] > S)
			{
				delfrD_edge(e1);
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
		return true;
	}

	void addtoS(int u)
	{
		//if (inS_v[u]) return;
		//printf("Stage adding %d: \n\n", u);

		for (int k = pstart[u]; k < pend[u]; k++)
		{
			int v = edges[k];
			int e = edgelist_pointer[k];
			cal_notadj(e);

			if (!deleted_edges[e] && inS_v[v] && notadj_e[e] > S)
			{
				//printf("!!!!%d %d-%d\n", notadj_e[369], edge_list[369 << 1], edge_list[(369 << 1) + 1]);
				delfrD_edge(e);
				deleted_edges[e] = 10;
			}
			//else if (!deleted_edges[e] && inS_v[v] && notadj_e[e] <= S)
			//{
			//	inS_e[e] = true;
			//	sedge.push_back(e);
			//}不能放在一起，因为有可能删掉一些边以后才>S，要在后面另起一行
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

		//for (int e = 0; e < (m >> 1); e++)
		//{
		//	if (!deleted_edges[e] && !check_common_neighbour(u, e))
		//		++notadj_e[e];
		//}

		inS_v[u] = true;
		svex.push_back(u);
		s_size++;


		//printf("svex: ");
		//for (ui i = 0; i < svex.size(); i++)
		//	printf("%d ", svex[i]);
		//printf("\n");

		//printf("sedge: ");
		//for (ui i = 0; i < sedge.size(); i++)
		//	printf("%d - %d and notadj_e = %d ", edge_list[sedge[i] << 1], edge_list[(sedge[i] << 1) + 1], notadj_e[sedge[i]]);
		//printf("\n\n");


	}

	void delfrS(int u)
	{
		//if (!inS_v[u]) return;
		//printf("Stage deling %d:\n\n", u);
		inS_v[u] = false;
		svex.pop_back();
		s_size--;

		//for (int e = 0; e < (m >> 1); e++)
		//{
		//	if (!deleted_edges[e] && !check_common_neighbour(u, e))
		//		--notadj_e[e];

		//	//if (deleted_flag[e] == 1)
		//	//{
		//	//	--notadj_e[e];
		//	//	deleted_flag[e] = 0;
		//	//}
		//}

		for (int k = pend[u] - 1; k >= pstart[u]; k--)//一定要倒着写，因为涉及到push_back和pop_back的顺序
		{
			int e = edgelist_pointer[k];
			int v = edges[k];
			if (inS_e[e])
			{
				inS_e[e] = false;
				sedge.pop_back();
			}
			else if (deleted_edges[e] == 10 && inS_v[v]) //u有可能有一些边，这些边的另外一个节点在S里面，但不是因为u的加入才删除的。所以把因为u的加入才删除的边标为10
				addtoD_edge(e);
		}



		//printf("svex: ");
		//for (ui i = 0; i < svex.size(); i++)
		//	printf("%d ", svex[i]);
		//printf("\n\n");

		//printf("sedge: ");
		//for (ui i = 0; i < sedge.size(); i++)
		//	printf("%d - %d and notadj_e = %d ", edge_list[sedge[i] << 1], edge_list[(sedge[i] << 1) + 1], notadj_e[sedge[i]]);
		//printf("\n\n");

	}

	bool check_adj_vertex(int u, int v)
	{
		for (int i = pstart[u]; i < pend[u]; i++)
			if (!deleted_vertex[edges[i]] && edges[i] == v) return true;
		return false;
	}

	int cal_notadj_v(int v)
	{
		int cnt = 0;
		for (auto u : svex)
			if (!check_adj_vertex(u, v)) cnt++;
		return cnt;
	}

	bool bb_s_truss_bar(int curS, int& size)
	{
		//pruning rules

		if (svex.size() >= LB)
		{
			size = svex.size();
			cout << "规模大于LB，返回正确，时间是：" << endl;
			return true;
		}

		if (curS < LB) return false;

		int min_tri = n, min_eid = -1, min_tri_s = n, min_eid_s = -1;
		for (int i = 0; i < (m >> 1); i++)
		{
			if (!deleted_edges[i] && tri_cnt[i] < min_tri)
			{
				min_tri = tri_cnt[i];
				min_eid = i;
			}
			if (inS_e[i] && tri_cnt[i] < min_tri_s)
			{
				min_tri_s = tri_cnt[i];
				min_eid_s = i;
			}
		}
		if (min_eid == -1) return false;

		if (min_eid_s != -1 && min_tri_s + 2 + S < LB) return false;

		//if (min_tri + 2 + S < LB) return false;

		//printf("tri[min_eid] = %d and curS = %d\n", tri_cnt[min_eid], curS);
		if (tri_cnt[min_eid] >= curS - 2 - S)
		{
			//output_min_k_truss();
			printf("The size of current vertex set is %d\n\n", svex.size());
			printf("There are %d vertices remaining\n\n", curS);
			size = curS;
			return true;
		}

		vector<int> del;
		for (int i = 0; i < n; i++)
		{
			if (!inS_v[i] && !deleted_vertex[i] && cal_notadj_v(i) >= S + 1)
				del.push_back(i);
		}

		//基于非邻居小于S+1的剪枝

		if (del.size()) {

			for (int x = 0; x < del.size(); x++)
			{
				delfrD_vertex(__LINE__, del[x]);
			}
			bool ret = bb_s_truss_bar(curS - del.size(), size);
			//printf("%d\n", curS - todel.size());
			for (int x = del.size() - 1; x >= 0; x--)
			{
				addtoD_vertex(__LINE__, del[x]);
			}
			return ret;
		}



		int max_tri = -1, max_eid = -1;
		for (int i = 0; i < (m >> 1); i++)
		{
			if (!deleted_edges[i] && inS_e[i])
			{
				cal_notadj(i);
				if (notadj_e[i] > max_tri)
				{
					max_tri = notadj_e[i];
					max_eid = i;
				}
			}
		}

		//基于vertex的剪枝
		if (max_eid != -1) cal_notadj(max_eid);
		if (max_eid != -1 && inS_e[max_eid] && notadj_e[max_eid] >= S) {
			if (notadj_e[max_eid] >= S + 1)
			{
				printf("%d %d\n", max_eid, notadj_e[max_eid]);
				return false;
			}
			std::vector <int> todel;
			for (int v = 0; v < n; v++) if (!deleted_vertex[v] && !inS_v[v] && !check_common_neighbour(v, max_eid)) {
				todel.push_back(v);
			}
			if (todel.size()) {

				for (int x = 0; x < todel.size(); x++)
				{
					delfrD_vertex(__LINE__, todel[x]);
				}
				bool ret = bb_s_truss_bar(curS - todel.size(), size);
				//printf("%d\n", curS - todel.size());
				for (int x = todel.size() - 1; x >= 0; x--)
				{
					addtoD_vertex(__LINE__, todel[x]);
				}
				return ret;
			}
		}
		//重新计算min_eid和max_eid，因为之前有可能被prune

		min_tri = n, min_eid = -1;
		for (int i = 0; i < (m >> 1); i++)
		{
			int u = edge_list[i << 1], v = edge_list[(i << 1) + 1];
			if (!deleted_edges[i] && !deleted_vertex[u] && !deleted_vertex[v] && tri_cnt[i] < min_tri)
			{
				min_tri = tri_cnt[i];
				min_eid = i;
			}
		}
		if (min_eid == -1) return false;

		//max_tri = -1, max_eid = -1;
		//for (int i = 0; i < (m >> 1); i++)
		//{
		//	if (!deleted_edges[i])
		//	{
		//		cal_notadj(i);
		//		if (notadj_e[i] > max_tri)
		//		{
		//			max_tri = notadj_e[i];
		//			max_eid = i;
		//		}
		//	}
		//}

		////基于edge的剪枝

		//if (tri_cnt[min_eid] < LB - 2 - S && !deleted_edges[min_eid]) //注意，不是LB+1 - 2 -S，因为如果有解到达了LB，我们就成功了
		//{
		//	if (inS_e[min_eid]) return false;//基于UB的剪枝

		//	delfrD_edge(min_eid);
		//	bool ret = bb_s_truss_bar(curS);
		//	addtoD_edge(min_eid);
		//	return ret;
		//}

  //      //基于edge的剪枝

  //       if (max_eid == -1) return false;

  //       cal_notadj(max_eid);
		// if (max_eid != -1 && !inS_e[max_eid] && notadj_e[max_eid] > S)
		// {

		//	 delfrD_edge(max_eid);
		//	 bool ret = bb_s_truss_bar(curS);
		//	 addtoD_edge(max_eid);
		//	 return ret;
		// }

		//开始分支

		vector<int> branch;
		for (int i = 0; i < n; i++)
		{
			if (!inS_v[i] && !deleted_vertex[i] && !check_common_neighbour(i, min_eid))
				branch.push_back(i);
		}

		if (branch.size() == 0) return false;

		random_shuffle(branch.begin(), branch.end());
		cal_notadj(min_eid);
		int canselect = S - notadj_e[min_eid];

		if (inS_e[min_eid])
		{
			cal_notadj(min_eid);
			//int canselect = S - notadj_e[min_eid];
			if (canselect < 0) return false;//检查一下
			if (canselect < 0)
			{
				printf("min_eid = %d\n", min_eid);
				printf("s_size = %d\n", s_size);
				printf("S = %d, notadj[min_eid] = %d\n", S, notadj_e[min_eid]);
				printf("svex: ");
				for (ui i = 0; i < svex.size(); i++)
					printf("%d ", svex[i]);
				printf("\n");

				printf("sedge: ");
				for (ui i = 0; i < sedge.size(); i++)
					printf("%d - %d and notadj_e = %d ", edge_list[sedge[i] << 1], edge_list[(sedge[i] << 1) + 1], notadj_e[sedge[i]]);
				printf("\n min_edge:");
				printf("%d - %d\n", edge_list[min_eid << 1], edge_list[(min_eid << 1) + 1]);
				printf("deleted? %d\n", deleted_edges[min_eid]);
			}
			assert(canselect >= 0);
			int pos = -1;
			bool ret = false;
			for (int i = 0; i < canselect && i < branch.size() && !ret; i++)
			{
				//printf("%d %d\n", branch.size(), i);
				delfrD_vertex(__LINE__, branch[i]);
				if (i && !canadd(branch[i - 1])) {
					addtoD_vertex(__LINE__, branch[i]);
					break;
				}
				if (i) {
					addtoS(branch[i - 1]);
					pos = i - 1;
				}
				ret |= bb_s_truss_bar(curS - 1, size);
				addtoD_vertex(__LINE__, branch[i]);
			}

			if (ret) {
				for (int i = pos; i >= 0; --i)
				{
					// printf("%d %d\n", branch.size(), i);
					delfrS(branch[i]);
				}
				return true;
			}

			//test
			//printf("svex size = %d\n", svex.size());
			//for (ui i = 0; i < svex.size(); i++) printf("%d ", svex[i]);
			//printf("\n");
			//for (ui i = 0; i < sedge.size(); i++) printf("%d-%d ", edge_list[sedge[i] << 1], edge_list[(sedge[i] << 1) + 1]);
			//printf("\n");
			//for (ui i = 0; i < sedge.size(); i++)
			//      printf("notadj = %d\n", notadj_e[sedge[i]]);

			//后续与minID不相连的节点都得删除
			for (int i = canselect; i < (int)branch.size(); ++i) {
				// printf("%d %d %d\n", branch.size(), i, canselect);//00000000000000000000000000000
				delfrD_vertex(__LINE__, branch[i]);
			}

			//最后一个分支

			//printf("%d %d\n", branch.size(), canselect);
			if (branch.size() >= canselect)
			{
				if (canselect == 0 || canadd(branch[canselect - 1])) {
					if (canselect) addtoS(branch[canselect - 1]);
					ret |= bb_s_truss_bar(curS - branch.size() + canselect, size);
					if (canselect) delfrS(branch[canselect - 1]);
				}
			}
			else
			{
				if (canselect == 0 || canadd(branch[branch.size() - 1])) {
					if (canselect) addtoS(branch[branch.size() - 1]);
					ret |= bb_s_truss_bar(curS, size);
					if (canselect) delfrS(branch[branch.size() - 1]);
				}
			}
			for (int i = pos; i >= 0; --i)
			{
				// printf("%d %d\n", branch.size(), i);
				delfrS(branch[i]);
			}
			for (int i = (int)branch.size() - 1; i >= (int)canselect; --i)
			{
				// printf("%d %d\n", branch.size(), i);
				addtoD_vertex(__LINE__, branch[i]);
			}
			return ret;
		}

		else
		{
			int u = edge_list[min_eid << 1], v = edge_list[(min_eid << 1) + 1];

			if ((inS_v[u] && !inS_v[v]) || (inS_v[v] && !inS_v[u]))
			{
				v = (inS_v[u] == 1) ? v : u;
				// printf("%d %d\n", min_eid, deleted_edges[min_eid]);
				delfrD_vertex(__LINE__, v);
				bool ret = bb_s_truss_bar(curS - 1, size);
				addtoD_vertex(__LINE__, v);

				//printf("s_size = %d\n", s_size);
				if (ret) return true;
				cal_notadj(min_eid);
				//int canselect = S - notadj_e[min_eid]; //canselect<=0，还要分支吗
				if (canselect < 0) return false;//我的评价是，不要了。相当于!canadd(v)
				if (!canadd(v)) return false;
				addtoS(v);

				/*		printf("%d-%d\n", u, v);
						for (ui i = 0; i < svex.size(); i++)
							printf("%d ", svex[i]);
						printf("\n");*/

				int pos = -1;
				for (int i = 0; !ret && i < canselect && i < branch.size(); ++i) {
					// printf("%d %d\n", branch.size(), i);
					delfrD_vertex(__LINE__, branch[i]);
					if (i && !canadd(branch[i - 1])) {
						addtoD_vertex(__LINE__, branch[i]);
						break;
					}
					if (i) {
						addtoS(branch[i - 1]);
						pos = i - 1;
					}
					ret |= bb_s_truss_bar(curS - 1, size);
					addtoD_vertex(__LINE__, branch[i]);
				}

				if (ret) {
					for (int i = pos; i >= 0; --i)
					{
						// printf("%d %d\n", branch.size(), i);
						delfrS(branch[i]);
					}
					delfrS(v);
					return true;
				}
				//printf("S = %d notadj = %d canselect = %d s_size = %d\n", S, notadj_e[min_eid], canselect, s_size);
				for (int i = canselect; i < (int)branch.size(); ++i) {
					//printf("%d %d\n", branch.size(), i);
					delfrD_vertex(__LINE__, branch[i]);
				}

				if (branch.size() >= canselect)
				{
					if (canselect == 0 || canadd(branch[canselect - 1])) {
						if (canselect) addtoS(branch[canselect - 1]);
						ret |= bb_s_truss_bar(curS - branch.size() + canselect, size);
						if (canselect) delfrS(branch[canselect - 1]);
					}
				}
				else
				{
					if (canselect == 0 || canadd(branch[branch.size() - 1])) {
						if (canselect) addtoS(branch[branch.size() - 1]);
						ret |= bb_s_truss_bar(curS, size);
						if (canselect) delfrS(branch[branch.size() - 1]);
					}
				}

				for (int i = 0; i <= pos; ++i)
				{
					//if (n == 42) printf("%d %d\n", branch.size(), i);
					delfrS(branch[i]);
				}
				for (int i = (int)branch.size() - 1; i >= (int)canselect; --i)
				{
					//if (n == 47) printf("%d %d\n", branch.size(), i);
					addtoD_vertex(__LINE__, branch[i]);
				}
				delfrS(v);
				return ret;
			}
			else if (!inS_v[u] && !inS_v[v])
			{
				//Branch 1，delete u and v
				delfrD_vertex(__LINE__, v);
				delfrD_vertex(__LINE__, u);
				bool ret = bb_s_truss_bar(curS - 2, size);
				addtoD_vertex(__LINE__, u);
				addtoD_vertex(__LINE__, v);
				if (ret) return true;

				//Branch 2, delete u, add v to S
				if (!canadd(v)) ret = false;
				else
				{
					delfrD_vertex(__LINE__, u);
					addtoS(v);
					ret = bb_s_truss_bar(curS - 1, size);
					delfrS(v);
					addtoD_vertex(__LINE__, u);
					if (ret) return true;
				}
				//Branch 3, delete v, add u to S

				if (!canadd(u)) ret = false;
				else
				{
					delfrD_vertex(__LINE__, v);
					addtoS(u);
					ret = bb_s_truss_bar(curS - 1, size);
					delfrS(u);
					addtoD_vertex(__LINE__, v);
					if (ret) return true;
				}
				//Branch 4, add u and v to S

				if (!canadd(u) || !canadd(v)) ret = false;
				else
				{
					cal_notadj(min_eid);
					//int canselect = S - notadj_e[min_eid];
					if (canselect < 0) return false;
					if (canadd(u)) addtoS(u);
					else return false;
					if (canadd(v)) addtoS(v);
					else
					{
						delfrS(u);
						return false;
					}

					int pos = -1;
					for (int i = 0; !ret && i < canselect && i < branch.size(); ++i) {
						//if (n == 42) printf("S = %d %d %d %d\n", S, branch.size(), i, canselect);
						delfrD_vertex(__LINE__, branch[i]);
						if (i && !canadd(branch[i - 1])) {
							addtoD_vertex(__LINE__, branch[i]);
							break;
						}
						if (i) {
							addtoS(branch[i - 1]);
							pos = i - 1;
						}
						ret |= bb_s_truss_bar(curS - 1, size);
						addtoD_vertex(__LINE__, branch[i]);
					}

					if (ret) {
						for (int i = pos; i >= 0; --i)
						{
							delfrS(branch[i]);
						}
						delfrS(v);
						delfrS(u);
						return true;
					}
					for (int i = canselect; i < (int)branch.size(); ++i) {
						delfrD_vertex(__LINE__, branch[i]);
					}

					//branch的大小有可能小于canselect=S-notadj[min_eid]
					//printf("%d %d S = %d notadj = %d\n", branch.size(), canselect, S, notadj_e[min_eid]);
					if (branch.size() >= canselect)
					{
						if (canselect == 0 || canadd(branch[canselect - 1])) {
							if (canselect) addtoS(branch[canselect - 1]);
							ret |= bb_s_truss_bar(curS - branch.size() + canselect, size);
							if (canselect) delfrS(branch[canselect - 1]);
						}
					}
					else
					{
						if (canselect == 0 || canadd(branch[branch.size() - 1])) {
							if (canselect) addtoS(branch[branch.size() - 1]);
							ret |= bb_s_truss_bar(curS, size);
							if (canselect) delfrS(branch[branch.size() - 1]);
						}
					}
					for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
					for (int i = (int)branch.size() - 1; i >= (int)canselect; --i)
					{
						//printf("%d %d\n", branch.size(), i);
						addtoD_vertex(__LINE__, branch[i]);
					}

					delfrS(v);
					delfrS(u);
					return ret;
				}
			}
		}
		return false;
	}

private:
	void output_min_k_truss()
	{
		vector<int> ans_v;
		for (int i = 0; i < n; i++)
			if (!deleted_vertex[i]) ans_v.push_back(i);
		ofstream outfile;
		outfile.open("C:\\Users\\zqf\\Desktop\\code_reachieve\\graph_visualize\\min_k_truss.txt");
		set<pair<int, int>> ans_e;
		for (auto v : ans_v)
		{
			for (int j = pstart[v]; j < pend[v]; j++)
			{
				int u = edges[j];
				if (!deleted_vertex[u])
				{
					if (v < u)
					{
						if (start_from_zero) ans_e.insert(mp(oldID[v], oldID[u]));
						else ans_e.insert(mp(oldID[v] + 1, oldID[u] + 1));

					}
					else
					{
						if (start_from_zero) ans_e.insert(mp(oldID[u], oldID[v]));
						else ans_e.insert(mp(oldID[u] + 1, oldID[v] + 1));
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

