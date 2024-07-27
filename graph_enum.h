#ifndef _GRAPHENUM_H_
#define _GRAPHENUM_H_

#include"Utility.h"
#include "Timer.h"

ui* vex_tri_cnt;

bool compare(ui u, ui v)
{
	return vex_tri_cnt[u] < vex_tri_cnt[v];
}

int timelimit = 86400;
class GraphE {
private:
	ui n;
	ept m;
	ui K;
	ui start_from_zero;//节点编号是否从0开始
	ui LB;//下界大小
	string dir;

	
	vector<int> nodes;
	queue<ui> Qe;//用于存放支持度小于K-2的边

	ui* inS_v;//数组，=1说明节点在已选集内，=0表示不在
	ui* inS_e;//=1表示边在已选集内

	//存储边的数据结构
	ept* pstart;
	ept* pend;
	ept* pend_buf;
	ui* edges;
	ui* edgelist_pointer;
	ui* edge_list;
	ui* tri_cnt;
	ui* degree;
	ui* newid;
	ui* old_id;
	ui* exist_oldid;

	//=1表示节点/边已经删除，=0表示未删除
	ui* deleted_vertex;
	ui* deleted_edge;

	vector<ui> svex;
	vector<ui> sedge;
	vector<ui> common_neighbors;
	ui* inS_vex;
	ui* inS_edge;
	ui* inS_degree;


	ui* tri_cnt_inS;
	ui* visit;

	//用于保存图
	ept* temp_pstart;
	ept* temp_pend;
	ept temp_m;
	ui temp_n;
	ui* temp_edgelist_pointer;
	ui* temp_edge_list;
	ui* temp_tri_cnt;
	ui* temp_degree;
	ui* temp_deleted_vertex;
	ui* temp_deleted_edge;
	ui* temp_edges;

	ui* s_tri_cnt;//边在已选集内的支持度
	ui* s_degree;//边在已选集内的度数



public:
	GraphE(const char* filename, const int _K)
	{
		K = _K;
		dir = string(filename);

		pstart = nullptr;
		pend = pend_buf = nullptr;
		edges = nullptr;

		edgelist_pointer = NULL;
		edge_list = NULL;
		tri_cnt = NULL;

	}
	~GraphE()
	{
		if (pstart != nullptr) {
			delete[] pstart;
			pstart = nullptr;
		}
		if (pend != nullptr) {
			delete[] pend;
			pend = nullptr;
		}
		if (pend_buf != NULL) {
			delete[] pend_buf;
			pend_buf = NULL;
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
		if (edges != nullptr) {
			delete[] edges;
			edges = nullptr;
		}

	}

	void heuristic_adding(clock_t start, bool &timeflag)
	{
		//删除孤立节点
		for (ui i = 0; i < n; i++)
			if (degree[i] == 0) deleted_vertex[i] = 1;

		clean(0);

		if (checkempty(deleted_vertex))
		{
			printf("图删完了");
			return;
		}

		svex.clear();
		sedge.clear();
		common_neighbors.clear();

		inS_vex = new ui[n];
		inS_edge = new ui[m >> 1];
		vex_tri_cnt = new ui[n];
		tri_cnt_inS = new ui[m >> 1];
		inS_degree = new ui[n];
		visit = new ui[n];
		memset(inS_vex, 0, sizeof(ui) * n);
		memset(inS_edge, 0, sizeof(ui) * (m >> 1));
		memset(vex_tri_cnt, 0, sizeof(ui) * n);
		memset(visit, 0, sizeof(ui) * n);
		memset(tri_cnt_inS, 0, sizeof(ui) * (m >> 1));
		memset(inS_degree, 0, sizeof(ui) * (n));


		int max_tri = -1, max_eid = -1;

		for (ui i = 0; i < (m >> 1); i++)
		{
			if ((int)(tri_cnt[i]) > max_tri)
			{
				max_tri = tri_cnt[i];
				max_eid = i;
			}
		}

		if (max_tri != -1)
		{
			ui u = edge_list[max_eid << 1], v = edge_list[(max_eid << 1) + 1];
			svex.push_back(u);
			svex.push_back(v);
			sedge.push_back(max_eid);
			inS_vex[u] = inS_vex[v] = 1;
			inS_edge[max_eid] = 1;
			tri_cnt_inS[max_eid] = 0;
			inS_degree[u] = inS_degree[v] = 1;
			push_common_neighbour(max_eid);
		}


		while (1)
		{
			clock_t finish;
			finish = clock();
			double time = (double)(finish - start) / CLOCKS_PER_SEC;
			if ((double)(finish - start) / CLOCKS_PER_SEC > 1800)
			{
				timeflag = false;
				return;
			}
			for (auto v : common_neighbors)
			{
				ui cnt = 0;
				for (auto e : sedge)
					if (check_common_neighbour(v, e)) cnt++;
				vex_tri_cnt[v] = cnt;
			}

			sort(common_neighbors.begin(), common_neighbors.end(), compare);

			ui u = common_neighbors.back();
			common_neighbors.pop_back();
			if (inS_vex[u]) continue;
			svex.push_back(u);
			//inS_degree[u] = 0;
			inS_vex[u] = 1;


			//加入u在S中的相关联的边，并计算支持度

			for (ui j = pstart[u]; j < pend[u]; j++)
			{
				ui e = edgelist_pointer[j];
				if (!deleted_edge[e] && inS_vex[edges[j]])
				{
					inS_degree[u]++;
					inS_degree[edges[j]]++;
					sedge.push_back(e);
					inS_edge[e] = 1;
					for (auto v : svex)
					{
						if (check_common_neighbour(v, e))
							tri_cnt_inS[e]++;
					}
					push_common_neighbour(e);
				}
			}

			//加入u之后，增加S中与之关联的边的支持度
			for (auto e : sedge)
			{
				if (check_common_neighbour(u, e))
					tri_cnt_inS[e]++;
			}
			//printf("111\n");
			//保存S，看看能否通过删边来得到一个ktruss

			vector<ui> temp_svex, temp_sedge;
			ui* temp_inS_degree = new ui[n];
			ui* temp_inS_vex = new ui[n];
			ui* temp_inS_edge = new ui[m >> 1];
			ui* temp_tri_cnt_inS = new ui[m >> 1];

			for (ui i = 0; i < n; i++)
			{
				temp_inS_degree[i] = inS_degree[i];
				temp_inS_vex[i] = inS_vex[i];
			}

			for (ui i = 0; i < (m >> 1); i++)
			{
				temp_inS_edge[i] = inS_edge[i];
				temp_tri_cnt_inS[i] = tri_cnt_inS[i];
			}


			temp_svex = svex, temp_sedge = sedge;
			queue<ui> del_queue;

			//for (auto e :temp_sedge)
			//{
			//	printf("%d\n", temp_tri_cnt_inS[e]);
			//}

			for (auto e : temp_sedge)
			{
				if (temp_tri_cnt_inS[e] < K - 2)
				{
					del_queue.push(e);
				}
			}

			while (!del_queue.empty())
			{
				ui e = del_queue.front();
				del_queue.pop();
				delfrD_edge_in_heu(e, temp_inS_edge, temp_inS_vex, temp_inS_degree, temp_tri_cnt_inS, del_queue);
			}

			/*		for (auto e : temp_sedge)
					{
						printf("%d\n", temp_tri_cnt_inS[e]);
					}*/

			for (auto v : temp_svex)
				if (temp_inS_vex[v] && temp_inS_degree[v] == 0) temp_inS_vex[v] = 0;

			//for (auto e : sedge)
			//{
			//	if (temp_inS_edge[e])
			//	{
			//		int u = edge_list[e << 1], v = edge_list[(e << 1) + 1];
			//		printf("%d %d\n", temp_inS_degree[u], temp_inS_degree[v]);
			//	}
			//}

			ui cnt = 0;
			for (auto v : temp_svex)
				if (temp_inS_vex[v]) cnt++;
			if (cnt != 0)
			{
				LB = cnt;
				printf("LB = %d\n", LB);
				return;
			}



			//if (exam_truss(sedge))
			//{
			//	for (auto v : svex) printf("%d ", v);
			//	printf("\n");
			//	LB = svex.size();
			//	printf("LB = %d\n", LB);
			//	return;
			//}
		}
	}

	void delfrD_edge_in_heu(int e, ui* temp_inS_edge, ui* temp_inS_vex, ui* temp_inS_degree, ui* temp_inS_tri_cnt, queue<ui>& del_queue)
	{
		if (!temp_inS_edge[e]) printf("%d-%d has been deleted\n", edge_list[e << 1], edge_list[(e << 1) + 1]);
		assert(temp_inS_edge[e]);

		temp_inS_edge[e] = 0;
		ui u = edge_list[e << 1], v = edge_list[(e << 1) + 1];
		if (temp_inS_degree[u] <= 0 || temp_inS_degree[v] <= 0) printf("degree u = %d and degree v = %d\n", temp_inS_degree[u], temp_inS_degree[v]);
		assert(temp_inS_degree[u] > 0);
		assert(temp_inS_degree[v] > 0);
		temp_inS_degree[u]--, temp_inS_degree[v]--;
		temp_inS_tri_cnt[e] = 0;

		ept ii = pstart[u], jj = pstart[v];
		ept u_n = pstart[u], v_n = pstart[v];

		while (true)
		{
			while (ii < pend[u] && !temp_inS_edge[edgelist_pointer[ii]]) ii++;
			while (jj < pend[v] && !temp_inS_edge[edgelist_pointer[jj]]) jj++;

			if (ii >= pend[u] || jj >= pend[v]) break;

			if (edges[ii] == edges[jj])
			{
				if ((temp_inS_tri_cnt[edgelist_pointer[ii]]--) == K - 2) del_queue.push(edgelist_pointer[ii]);
				if ((temp_inS_tri_cnt[edgelist_pointer[jj]]--) == K - 2) del_queue.push(edgelist_pointer[jj]);
				ii++;
				jj++;
			}
			else if (edges[ii] < edges[jj])
			{
				ii++;
			}
			else if (edges[ii] > edges[jj])
			{
				jj++;
			}

		}

		if (temp_inS_degree[u] == 0) temp_inS_vex[u] = 0;
		if (temp_inS_degree[v] == 0) temp_inS_vex[v] = 0;
	}

	bool check_common_neighbour(ui v, ui e)
	{
		//printf("%d\n", e);
		ui flag = 0;
		ui u = edge_list[e << 1], w = edge_list[(e << 1) + 1];
		if (u == v || w == v) return false;//我靠，这句话太重要了
		//printf("%d\n", e);
		for (ept i = pstart[u]; i < pend[u]; i++)
		{
			if (edges[i] == v && !deleted_edge[edgelist_pointer[i]])
			{
				flag = 1;
				break;
			}
		}

		//printf("%d %d\n", u, w);
		for (ept i = pstart[w]; i < pend[w]; i++)
		{

			if (edges[i] == v && flag == 1 && !deleted_edge[edgelist_pointer[i]])
			{
				flag = 2;
				break;
			}
		}

		if (flag == 2) return true;
		return false;
	}

	void push_common_neighbour(ui e)
	{
		ui u = edge_list[e << 1], v = edge_list[(e << 1) + 1];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted_edge[edgelist_pointer[j]] && !inS_vex[edges[j]]) visit[edges[j]] = 1;

		for (ept j = pstart[v]; j < pend[v]; j++)
		{
			if (!deleted_edge[edgelist_pointer[j]] && visit[edges[j]] && !inS_vex[edges[j]])
				common_neighbors.push_back(edges[j]);
		}

		for (ept j = pstart[u]; j < pend[u]; j++) visit[edges[j]] = 0;
	}


	int find_min_k_truss_add(clock_t start, bool timeflag)
	{
		printf("start algorithm\n");
		//读图
		read_graph();

		//对输入的图进行处理（计算三角形等）
		process_graph();

		//存储图
		store_graph();

		LB = n;
		//用启发式算法计算一个下界LB
		//heuristic();

		heuristic_adding(start, timeflag);
		//return LB;

		//return LB;
		if (LB == n)
		{
			printf("当K = %d时，图已经被剪完了。", K);
			return -1;
		}

		if (LB == K)
		{
			return K;
		}
		//else return timelimit;


		//恢复图
		recover_graph();

		//预处理，删去支持度小于K-2的边
		clean(0);

		//重新安排节点和边的编号
		rebulid_graph(tri_cnt, edge_list, edgelist_pointer);
		pstart[n] = pend[n - 1];
		//LB = n;



		//开始枚举
		inS_v = new ui[n];
		inS_e = new ui[m >> 1];
		s_tri_cnt = new ui[m >> 1];
		s_degree = new ui[n];
		//printf("111\n");
		memset(inS_v, 0, sizeof(ui) * n);
		memset(inS_e, 0, sizeof(ui) * (m >> 1));
		memset(s_tri_cnt, 0, sizeof(ui) * (m >> 1));
		memset(s_degree, 0, sizeof(ui) * n);


		vector<ui> cur_vset;
		cur_vset.push_back(0);
		inS_v[0] = 1;
		enumeration(start, n, m, 0, LB, cur_vset);
		cur_vset.pop_back();
		inS_v[0] = 0;
		enumeration(start, n, m, 0, LB, cur_vset);

		clock_t finish;
		finish = clock();
		double time = (double)(finish - start) / CLOCKS_PER_SEC;
		if (time > timelimit)
		{
			timeflag = false;
			return -1;
		}
		cout << "min size is: " << LB << endl;
		return LB;
	}

	int find_min_k_truss_del(clock_t start, bool timeflag)
	{
		//读图
		read_graph();

		//对输入的图进行处理（计算三角形等）
		process_graph();

		//存储图
		store_graph();

		LB = n;
		//用启发式算法计算一个下界LB
		heuristic();

		//heuristic_adding(start, timeflag);

		return LB;

		if (LB == n)
		{
			printf("当K = %d时，图已经被剪完了。", K);
			return -1;
		}

		if (LB == K)
		{
			return K;
		}

		//恢复图
		recover_graph();

		//预处理，删去支持度小于K-2的边
		clean(0);

		//重新安排节点和边的编号
		rebulid_graph(tri_cnt, edge_list, edgelist_pointer);
		pstart[n] = pend[n - 1];
		//LB = n;

		//开始枚举
		inS_v = new ui[n];
		inS_e = new ui[m >> 1];
		s_tri_cnt = new ui[m >> 1];
		s_degree = new ui[n];
		//printf("111\n");
		memset(inS_v, 0, sizeof(ui) * n);
		memset(inS_e, 0, sizeof(ui) * (m >> 1));
		memset(s_tri_cnt, 0, sizeof(ui) * (m >> 1));
		memset(s_degree, 0, sizeof(ui) * n);


		vector<ui> cur_vset;
		cur_vset.push_back(0);
		inS_v[0] = 1;
		enumeration(start, n, m, 0, LB, cur_vset);
		cur_vset.pop_back();
		inS_v[0] = 0;
		enumeration(start, n, m, 0, LB, cur_vset);

		clock_t finish;
		finish = clock();
		double time = (double)(finish - start) / CLOCKS_PER_SEC;
		if (time > timelimit)
		{
			timeflag = false;
			return -1;
		}
		cout << "min size is: " << LB << endl;
		return LB;
	}

	int find_min_k_truss_enum(clock_t start, bool timeflag)
	{
		//读图
		read_graph();

		//对输入的图进行处理（计算三角形等）
		process_graph();
		store_graph();

		LB = n;
		//预处理，删去支持度小于K-2的边

		recover_graph();
		clean(0);

		//重新安排节点和边的编号
		rebulid_graph(tri_cnt, edge_list, edgelist_pointer);
		pstart[n] = pend[n - 1];
		//LB = n;

		//cout << n << endl << endl;



		//开始枚举
		inS_v = new ui[n];
		inS_e = new ui[m >> 1];
		s_tri_cnt = new ui[m >> 1];
		s_degree = new ui[n];
		//printf("111\n");
	    memset(inS_v, 0, sizeof(ui) * n);
		memset(inS_e, 0, sizeof(ui) * (m >> 1));
		memset(s_tri_cnt, 0, sizeof(ui) * (m >> 1));
		memset(s_degree, 0, sizeof(ui) * n);

		//按随机顺序对顶点排序
		nodes.resize(n);
		for (int i = 0; i < n; i++)
			nodes[i] = i;
		std::random_device rd;  // 用于获取随机种子
		std::mt19937 gen(rd()); // 使用Mersenne Twister引擎作为随机数生成器
		std::shuffle(nodes.begin(), nodes.end(), gen);  // 对节点编号打乱
		/*for (int i = 0; i < n; i++)
			cout << nodes[i] << endl;*/

		

		vector<ui> cur_vset;
		cur_vset.push_back(nodes[0]);
		inS_v[nodes[0]] = 1;
		
		enumeration(start, n, m, 0, LB, cur_vset);
		cur_vset.pop_back();
		inS_v[nodes[0]] = 0;
		enumeration(start, n, m, 0, LB, cur_vset);

		clock_t finish;
		finish = clock();
		double time = (double)(finish - start) / CLOCKS_PER_SEC;
		if (time > timelimit)
		{
			timeflag = false;
			return -1;
		}
		cout << "min size is: " << LB << endl;
		return LB;
	}

	int find_min_k_truss_query(clock_t start, bool timeflag, int q)
	{
		//读图
		read_graph();

		//对输入的图进行处理（计算三角形等）
		process_graph();
		store_graph();

		LB = n;
		//cout << n << endl << endl;
		//heuristic_adding(start, timeflag);
		//预处理，删去支持度小于K-2的边
		//heuristic_adding(start, timeflag);

		//heuristic_query(start, timeflag, q);

		recover_graph();
		clean(0);

		//重新安排节点和边的编号
		rebulid_graph(tri_cnt, edge_list, edgelist_pointer);
		pstart[n] = pend[n - 1];
		//LB = n;

		int count = 0;
		for (int i = 0; i < n; i++)
		{
			if (count == 5) break;
			if (exist_oldid[i])
			{
				cout << i << " ";
				count++;
			}
		}
		cout << endl;

		//cout << n << endl << endl;

		//判断q是否还在图中
		if (!exist_oldid[q])
		{
			printf("q has been deleted\n");
			return -1;
		}

		

		//赋予q新的编号
		q = newid[q];


		//开始枚举
		inS_v = new ui[n];
		inS_e = new ui[m >> 1];
		s_tri_cnt = new ui[m >> 1];
		s_degree = new ui[n];
		//printf("111\n");
		memset(inS_v, 0, sizeof(ui) * n);
		memset(inS_e, 0, sizeof(ui) * (m >> 1));
		memset(s_tri_cnt, 0, sizeof(ui) * (m >> 1));
		memset(s_degree, 0, sizeof(ui) * n);




		vector<ui> cur_vset;
		cur_vset.push_back(q);
		inS_v[q] = 1;
		
		/*cur_vset.push_back(0);
		inS_v[0] = 1;*/

		enumeration_query(start, n, m, 0, LB, cur_vset, q);
		//cur_vset.pop_back();
		//inS_v[0] = 0;
		//enumeration(start, n, m, 0, LB, cur_vset);

		clock_t finish;
		finish = clock();
		double time = (double)(finish - start) / CLOCKS_PER_SEC;
		if (time > timelimit)
		{
			timeflag = false;
			return -1;
		}
		cout << "min size is: " << LB << endl;
		return LB;
	}



private:
	void oriented_triangle_counting(ui n, ui m, ept* pstart, ept* pend, ui* edges, ui* tri_cnt)
	{
		for (ui i = 0; i < n; i++)
		{
			ept& end = pend[i] = pstart[i];
			for (ui j = pstart[i]; j < pstart[i + 1]; j++)
			{
				if (edges[j] > i) edges[end++] = edges[j];
			}
		}

		ui* adj = new ui[n];
		ui cnt = 0;
		memset(tri_cnt, 0, m * sizeof(ui));
		memset(adj, 0, n * sizeof(ui));

		for (ui i = 0; i < n; i++)
		{
			for (ui u = pstart[i]; u < pend[i]; u++) adj[edges[u]] = u + 1;

			for (ui u = pstart[i]; u < pend[i]; u++)
			{
				ui v = edges[u];
				for (ui j = pstart[v]; j < pend[v]; j++)
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

			for (ui u = pstart[i]; u < pend[i]; u++) adj[edges[u]] = 0;
		}

#ifndef  NDEBUG

		//for (ui i = 0; i < n; i++)
		//{
		//	for (ui u = pstart[i]; u < pend[i]; u++)
		//		cout << tri_cnt[u] << endl;


		//}

		//printf("Total number of triangles:%d\n", cnt);
#endif // ! NDEBUG
	}
	void reorganize_oriented_graph(ui n, ui* tri_cnt, ui* edge_list, ept* pstart, ept* pend, ept* pend2, ui* edges, ui* edgelist_pointer)
	{
		memset(tri_cnt, 0, (m / 2) * sizeof(ui));
		for (ui i = 0; i < n; i++) pend2[i] = pend[i];


		ept pos = 0;

		for (ui i = 0; i < n; i++)
		{
			for (ept j = pstart[i]; j < pend[i]; j++)
			{
				tri_cnt[pos >> 1] = edgelist_pointer[j];
				edge_list[pos++] = i, edge_list[pos++] = edges[j];

				ept& k = pend2[edges[j]];
				edgelist_pointer[j] = edgelist_pointer[k] = (pos >> 1) - 1;
				edges[k++] = i;
			}
		}



#ifndef NDEBUG
		for (ui i = 0; i < n; i++) assert(pend2[i] == pstart[i + 1]);
#endif // !NDEBUG


		for (ui i = 0; i < n; i++)
		{
			pend2[i] = pend[i];
			pend[i] = pstart[i];
		}

		for (ui i = 0; i < n; i++)
		{
			for (ept j = pend2[i]; j < pstart[i + 1]; j++)
			{
				ept& k = pend[edges[j]];
				edgelist_pointer[k] = edgelist_pointer[j];
				edges[k++] = i;
			}
		}

		ept* ids = new ept[m];
		ui* buf = new ui[m];

		for (ui i = 0; i < n; i++)
		{
			if (pend[i] == pstart[i] || pend[i] == pstart[i + 1]) continue;
			ept k = pend[i], j = pstart[i], pos = 0;
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

			for (ept j = 0; j < pos; j++)
			{
				edges[pstart[i] + j] = ids[j];
				edgelist_pointer[pstart[i] + j] = buf[j];
			}
		}

	}
	void triangle_counting_inS(ui n, ui m, ept* pstart, ept* pend, ui* edges, ui* s_tri_cnt, vector<ui> svex)
	{
		ui* tmp_edges = new ui[m];
		ui* tmp_edgelist_pointer = new ui[m];
		for (ui i = 0; i < m; i++)
		{
			tmp_edgelist_pointer[i] = edgelist_pointer[i];
			tmp_edges[i] = edges[i];
		}

		for (ui i = 0; i < svex.size(); i++)
		{
			ui w = svex[i];
			ept& end = pend[w] = pstart[w];
			for (ui j = pstart[w]; j < pstart[w + 1]; j++)
			{
				if (edges[j] > w)
				{
					edgelist_pointer[end] = edgelist_pointer[j];
					edges[end++] = edges[j];
				}
			}
		}

	

		ui* adj = new ui[n];
		ui cnt = 0;
		memset(s_tri_cnt, 0, (m >> 1) * sizeof(ui));
		memset(adj, 0, n * sizeof(ui));

		for (ui i = 0; i < svex.size(); i++)
		{
			ui w = svex[i];
			for (ui u = pstart[w]; u < pend[w]; u++)
				if (inS_v[edges[u]]) adj[edges[u]] = u + 1;

			for (ui u = pstart[w]; u < pend[w]; u++)
			{
				ui v = edges[u];
				
				if (!inS_v[v]) continue;
				for (ui j = pstart[v]; j < pend[v]; j++)
				{
					//printf("%d %d %d %d\n", n, v, m, j);
					if (!inS_v[edges[j]]) continue;
					//printf("%d %d\n", n, edges[j]);
					if (adj[edges[j]])
					{
						s_tri_cnt[edgelist_pointer[j]]++;
						s_tri_cnt[edgelist_pointer[u]]++;
						s_tri_cnt[edgelist_pointer[adj[edges[j]] - 1]]++;
						cnt++;

						/*ui e1 = edgelist_pointer[j], e2 = edgelist_pointer[u], e3 = edgelist_pointer[adj[edges[j]] - 1];
						printf("%d-%d %d-%d %d-%d\n", edge_list[e1 << 1], edge_list[(e1 << 1) + 1], edge_list[(e2 << 1)], edge_list[(e2 << 1) + 1], edge_list[(e3 << 1)], edge_list[(e3 << 1) + 1]);*/
					}
				}
			}

			for (ui u = pstart[w]; u < pend[w]; u++) adj[edges[u]] = 0;
		}

		for (ui i = 0; i < svex.size(); i++)
		{
			ui w = svex[i];
			pend[w] = pstart[w + 1];
			//for (ui j = pstart[w]; j < pend[w]; j++)
			//{
			//	edges[j]= tmp_edges[j];
			//	edgelist_pointer[j] = tmp_edgelist_pointer[j];
			//}
		}

		for (ui i = 0; i < m; i++)
		{
			edgelist_pointer[i] = tmp_edgelist_pointer[i];
			edges[i] = tmp_edges[i];
		}

		delete[] tmp_edges;
		delete[] tmp_edgelist_pointer;
		delete[] adj;

#ifndef  NDEBUG
		//printf("Total number of triangles:%d\n", cnt);
#endif // ! NDEBUG

	}
	void enumeration(clock_t start, ui n, ui m, ui v, ui& lb, vector<ui> cur_vset) //ui* inS_v, ui* inS_e
	{
		//cout << 1111111111 << endl;
		clock_t finish;
		finish = clock();
		double time = (double)(finish - start) / CLOCKS_PER_SEC;
		if (time > timelimit)
		{
			return;
		}
	/*	cout << endl << "cur_vset.size " << cur_vset.size() << endl;
		cout << "They are: ";
		for (ui i = 0; i < cur_vset.size(); i++)
			cout << cur_vset[i] << " ";
		cout << endl << endl;*/

		//当前已选集大小超过LB，剪枝
		if (cur_vset.size() >= lb)
		{
			//printf("Large Size\n");
			return;
		}

		if (lb == K)
		{
			printf("Already find the optimal solution with size %d\n", K);
			return;
		}

		memset(s_tri_cnt, 0, sizeof(ui) * (m >> 1));

		//计算S内的边的支持度
		triangle_counting_inS(n, m, pstart, pend, edges, s_tri_cnt, cur_vset);




		//保存已选集的图
		ui* temp_inS_e = new ui[m >> 1];
		ui* temp_s_tri_cnt = new ui[m >> 1];
		ui* temp_inS_v = new ui[n];
		ui* temp_s_degree = new ui[n];
		for (ui i = 0; i < n; i++)
		{
			temp_inS_v[i] = inS_v[i];
			temp_s_degree[i] = s_degree[i];
		}
		for (ui i = 0; i < (m >> 1); i++)
		{
			temp_inS_e[i] = inS_e[i];
			temp_s_tri_cnt[i] = s_tri_cnt[i];
		}

		//for (ui i = 0; i < (m >> 1); i++)
		//{
		//	if (inS_e[i])
		//	{
		//		printf("K - 2 = %d, edge %d-%d. tri_cnt = %d\n",K-2, edge_list[i << 1], edge_list[(i << 1) + 1], s_tri_cnt[i]);
		//	}
		//}

		//删去支持度小于K-2的边
		while (!Qe.empty()) Qe.pop();

		for (ui i = 0; i < (m >> 1); i++)
		{
			if (temp_inS_e[i] && temp_s_tri_cnt[i] < K - 2)
			{
				//printf("%d-%d inque\n", edge_list[i << 1], edge_list[(i << 1) + 1]);
				Qe.push(i);
			}
		}

		while (!Qe.empty())
		{
			ui e = Qe.front();
			Qe.pop();
			delfrD_edge_S(e, temp_inS_e, temp_inS_v, temp_s_degree, temp_s_tri_cnt);
		}

		for (auto v : cur_vset)
		{
			if (temp_inS_v[v] && temp_s_degree[v] == 0)
				temp_inS_v[v] = 0;
			//printf("vertex %d. inS= %d temp_s_degree = %d\n", v, temp_inS_v[v], temp_s_degree[v]);
		}

		//若图没被删完，说明此时已经得到一个ktruss，保存大小，返回
		if (check_truss(cur_vset, temp_inS_v))
		{
			bool clique = true;
			ui min_size = 0;

			for (auto v : cur_vset)
			{
				if (temp_inS_v[v])
				{
					cout << temp_s_degree[v] << " ";
					min_size++;
				}
			}
			cout << endl;

			for (auto v : cur_vset)
			{
				if (temp_inS_v[v] && temp_s_degree[v] < (min_size - 1))
					clique = false;
			}

			
			lb = min_size;
			printf("Get smaller size: %d\n", lb);

			if (lb == K || clique == true)
			{
				ofstream outfile;
				outfile.open("/home/zhangqifan/min_k_truss/dataset/result_enum_Email.txt");
				finish = clock();
				printf("The running time is: %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
				outfile << " " << "ENUM" << " " << "K = " << K << " " << (double)(finish - start) / CLOCKS_PER_SEC << "s" << " " << "size is" << " " << K << endl;
				outfile.close();
				exit(0);
			}
			return;
		}

		delete[] temp_inS_e;
		delete[] temp_s_tri_cnt;
		delete[] temp_inS_v;
		delete[] temp_s_degree;


		//枚举
		for (ui i = v + 1; i < n; i++)
		{
			cur_vset.push_back(nodes[i]);
			inS_v[nodes[i]] = true;

			//加入S内和v相连的边，并更新s_degree
			for (ept j = pstart[nodes[i]]; j < pend[nodes[i]]; j++)
			{
				if (inS_v[edges[j]])
				{
					s_degree[nodes[i]]++;
					s_degree[edges[j]]++;
					inS_e[edgelist_pointer[j]] = 1;
				}
			}

			if (cur_vset.size() < lb)
				enumeration(start, n, m, i, lb, cur_vset);


			cur_vset.pop_back();
			inS_v[nodes[i]] = false;
			s_degree[nodes[i]] = 0;
			//删除S内和v相连的边，并更新s_degree
			for (ept j = pstart[nodes[i]]; j < pend[nodes[i]]; j++)
			{
				ui e = edgelist_pointer[j];
				inS_e[edgelist_pointer[j]] = 0;
				if (inS_v[edges[j]]) s_degree[edges[j]]--;
			}

			if (cur_vset.size() >= lb) return;
		}

		return;
	}

	void enumeration_query(clock_t start, ui n, ui m, ui v, ui& lb, vector<ui> cur_vset, int q) //ui* inS_v, ui* inS_e
	{
		//cout << 1111111111 << endl;
		clock_t finish;
		finish = clock();
		double time = (double)(finish - start) / CLOCKS_PER_SEC;
		if (time > timelimit)
		{
			return;
		}
		/*	cout << endl << "cur_vset.size " << cur_vset.size() << endl;
			cout << "They are: ";
			for (ui i = 0; i < cur_vset.size(); i++)
				cout << cur_vset[i] << " ";
			cout << endl << endl;*/

			//当前已选集大小超过LB，剪枝
		if (cur_vset.size() >= lb)
		{
			//printf("Large Size\n");
			return;
		}

		if (lb == K)
		{
			printf("Already find the optimal solution with size %d\n", K);
			return;
		}

		memset(s_tri_cnt, 0, sizeof(ui) * (m >> 1));

		//计算S内的边的支持度
		triangle_counting_inS(n, m, pstart, pend, edges, s_tri_cnt, cur_vset);




		//保存已选集的图
		ui* temp_inS_e = new ui[m >> 1];
		ui* temp_s_tri_cnt = new ui[m >> 1];
		ui* temp_inS_v = new ui[n];
		ui* temp_s_degree = new ui[n];
		for (ui i = 0; i < n; i++)
		{
			temp_inS_v[i] = inS_v[i];
			temp_s_degree[i] = s_degree[i];
		}
		for (ui i = 0; i < (m >> 1); i++)
		{
			temp_inS_e[i] = inS_e[i];
			temp_s_tri_cnt[i] = s_tri_cnt[i];
		}

		//for (ui i = 0; i < (m >> 1); i++)
		//{
		//	if (inS_e[i])
		//	{
		//		printf("K - 2 = %d, edge %d-%d. tri_cnt = %d\n",K-2, edge_list[i << 1], edge_list[(i << 1) + 1], s_tri_cnt[i]);
		//	}
		//}

		//删去支持度小于K-2的边
		while (!Qe.empty()) Qe.pop();

		for (ui i = 0; i < (m >> 1); i++)
		{
			if (temp_inS_e[i] && temp_s_tri_cnt[i] < K - 2)
			{
				//printf("%d-%d inque\n", edge_list[i << 1], edge_list[(i << 1) + 1]);
				Qe.push(i);
			}
		}

		while (!Qe.empty())
		{
			ui e = Qe.front();
			Qe.pop();
			delfrD_edge_S(e, temp_inS_e, temp_inS_v, temp_s_degree, temp_s_tri_cnt);
		}

		for (auto v : cur_vset)
		{
			if (temp_inS_v[v] && temp_s_degree[v] == 0)
				temp_inS_v[v] = 0;
			//printf("vertex %d. inS= %d temp_s_degree = %d\n", v, temp_inS_v[v], temp_s_degree[v]);
		}

		//若图没被删完，说明此时已经得到一个ktruss，保存大小，返回
		if (check_truss(cur_vset, temp_inS_v) && temp_inS_v[q])
		{
			bool clique = true;
			ui min_size = 0;

			for (auto v : cur_vset)
			{
				if (temp_inS_v[v])
				{
					cout << temp_s_degree[v] << " ";
					min_size++;
				}
			}
			cout << endl;

			for (auto v : cur_vset)
			{
				if (temp_inS_v[v] && temp_s_degree[v] < (min_size - 1))
					clique = false;
			}


			lb = min_size;
			printf("Get smaller size: %d\n", lb);





			if (lb == K || clique == true)
			{
				ofstream outfile;
				outfile.open("/home/zhangqifan/min_k_truss/dataset/result_enum_Email.txt");
				finish = clock();
				printf("The running time is: %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
				outfile << " " << "ENUM" << " " << "K = " << K << " " << (double)(finish - start) / CLOCKS_PER_SEC << "s" << " " << "size is" << " " << K << endl;
				outfile.close();
				exit(0);
			}
			return;
		}

		delete[] temp_inS_e;
		delete[] temp_s_tri_cnt;
		delete[] temp_inS_v;
		delete[] temp_s_degree;


		//枚举
		for (ui i = v; i < n; i++)
		{
			if (i == q) continue;

			cur_vset.push_back(i);
			inS_v[i] = true;

			//加入S内和v相连的边，并更新s_degree
			for (ept j = pstart[i]; j < pend[i]; j++)
			{

				if (inS_v[edges[j]])
				{
					s_degree[i]++;
					s_degree[edges[j]]++;
					inS_e[edgelist_pointer[j]] = 1;
				}
			}

			if (cur_vset.size() < lb && (i + 1 < n))
				enumeration_query(start, n, m, i + 1, lb, cur_vset, q);


			cur_vset.pop_back();
			inS_v[i] = false;
			s_degree[i] = 0;
			//删除S内和v相连的边，并更新s_degree
			for (ept j = pstart[i]; j < pend[i]; j++)
			{
				ui e = edgelist_pointer[j];
				inS_e[edgelist_pointer[j]] = 0;
				if (inS_v[edges[j]]) s_degree[edges[j]]--;
			}

			if (cur_vset.size() >= lb) return;
		}

		return;
	}

	void delfrD_edge_S(ui e, ui* inS_e, ui* inS_v, ui* degree, ui* tri_cnt)
	{
		//printf("%d-%d has been deleted\n", edge_list[e << 1], edge_list[(e << 1) + 1]);
		if (!inS_e[e]) printf("%d-%d has been deleted\n", edge_list[e << 1], edge_list[(e << 1) + 1]);
		assert(inS_e[e]);

		inS_e[e] = 0;
		ui u = edge_list[e << 1], v = edge_list[(e << 1) + 1];
		if (degree[u] <= 0 || degree[v] <= 0) printf("degree u = %d and degree v = %d\n", degree[u], degree[v]);
		assert(degree[u] > 0);
		assert(degree[v] > 0);
		degree[u]--, degree[v]--;

		int ii = pstart[u], jj = pstart[v];

		while (true)
		{
			while (ii < pend[u] && !inS_e[edgelist_pointer[ii]]) ii++;
			while (jj < pend[v] && !inS_e[edgelist_pointer[jj]]) jj++;

			if (ii >= pend[u] || jj >= pend[v]) break;
			if (edges[ii] == edges[jj])
			{
				if ((tri_cnt[edgelist_pointer[ii]]--) == K - 2) Qe.push(edgelist_pointer[ii]);
				if ((tri_cnt[edgelist_pointer[jj]]--) == K - 2) Qe.push(edgelist_pointer[jj]);
				ii++;
				jj++;
			}
			else if (edges[ii] < edges[jj]) ii++;
			else if (edges[ii] > edges[jj]) jj++;
		}

		if (degree[u] == 0) inS_v[u] = 0;
		if (degree[v] == 0) inS_v[v] = 0;
	}
	ui count_vertices(ui n, ept m, ui* edge_list, bool* visit)//填写单向m
	{
		vector<int> cnt;
		ui* exist = new ui[n];
		memset(exist, false, sizeof(ui) * n);
		for (ept idx = 0; idx < m; idx++)
		{
			if (!visit[idx])
			{
				ui u = edge_list[idx << 1], v = edge_list[(idx << 1) + 1];
				if (!exist[u])
				{
					exist[u] = true;
					cnt.push_back(u);
				}
				if (!exist[v])
				{
					exist[v] = true;
					cnt.push_back(v);
				}
			}
		}
		return cnt.size();
	}

	bool check_neighbour(ui u, ui v)
	{
		for (ept i = pstart[u]; i < pend[u]; i++)
		{
			if (edges[i] == v) return true;
		}
		return false;
	}

	bool check_truss(vector<ui> svex, ui* temp_inS_v)
	{
		for (auto v : svex)
		{
			if (temp_inS_v[v]) return true;
		}
		return false;
	}

	void read_graph()
	{
		FILE* f = Utility::open_file((dir).c_str(), "r");
		fscanf(f, "%u%lu%u", &n, &m, &start_from_zero);

		printf("n = %d, m = %d, K = %d\n", n, m, K);
		m = m * 2;
		vector<pair<ui, ui>> vp;

		for (ui i = 0; i < m / 2; i++)
		{
			ui a, b;
			fscanf(f, "%u%u", &a, &b);
			if (!start_from_zero)
			{
				a = a - 1;
				b = b - 1;
			}
			vp.pb(mp(a, b));
			vp.pb(mp(b, a));
		}

		sort(vp.begin(), vp.end());

		if (pstart != nullptr) delete[] pstart;
		pstart = new ept[n + 1];
		if (edges != nullptr) delete[] edges;
		edges = new ui[m];
		pstart[0] = 0;
		ept idx = 0;

		for (ui i = 0; i < n; i++)
		{
			pstart[i + 1] = pstart[i];
			while (idx < vp.size() && vp[idx].first == i) edges[pstart[i + 1]++] = vp[idx++].second;
		}

		fclose(f);

		printf("Finish reading graph!\n");
	}

	void process_graph()
	{
		edgelist_pointer = new ui[m];
		pend = new ept[n];
		memset(edgelist_pointer, 0, m * sizeof(ui));
		oriented_triangle_counting(n, m, pstart, pend, edges, edgelist_pointer);

		edge_list = new ui[m];
		tri_cnt = new ui[m >> 1];
		ept* pend_buf = new ept[n];

		reorganize_oriented_graph(n, tri_cnt, edge_list, pstart, pend, pend_buf, edges, edgelist_pointer);

		degree = new ui[n];
		deleted_vertex = new ui[n];
		deleted_edge = new ui[m >> 1];

		memset(degree, 0, sizeof(ui) * n);
		memset(deleted_vertex, 0, sizeof(ui) * n);
		memset(deleted_edge, 0, sizeof(ui) * (m >> 1));

		for (ui i = 0; i < n; i++)
		{
			pend[i] = pstart[i + 1];
			degree[i] = pend[i] - pstart[i];
			//printf("i = %d degree[i] = %d\n", i, degree[i]);
		}
	}

	void delfrD_edge(int e)
	{
		//printf("%d-%d has been deleted\n", edge_list[e << 1], edge_list[(e << 1) + 1]);
		if (deleted_edge[e]) printf("%d-%d has been deleted\n", edge_list[e << 1], edge_list[(e << 1) + 1]);
		assert(!deleted_edge[e]);

		deleted_edge[e] = 1;
		ui u = edge_list[e << 1], v = edge_list[(e << 1) + 1];
		if (degree[u] <= 0 || degree[v] <= 0) printf("degree u = %d and degree v = %d\n", degree[u], degree[v]);
		assert(degree[u] > 0);
		assert(degree[v] > 0);
		degree[u]--, degree[v]--;
		//printf("After--, degree u = %d and degree v = %d\n", degree[u], degree[v]);

		ept ii = pstart[u], jj = pstart[v];
		ept u_n = pstart[u], v_n = pstart[v];

		//for (ept j = pstart[3]; j < pend[3]; j++)
		//{
		//	ui e = edgelist_pointer[j];
		//	printf("%d - %d\n", edge_list[e << 1], edge_list[(e << 1) + 1]);
		//}

		while (true)
		{
			while (ii < pend[u] && deleted_edge[edgelist_pointer[ii]]) ii++;
			while (jj < pend[v] && deleted_edge[edgelist_pointer[jj]]) jj++;

			if (ii >= pend[u] || jj >= pend[v]) break;

			if (edges[ii] == edges[jj])
			{
				edges[u_n] = edges[ii], edgelist_pointer[u_n++] = edgelist_pointer[ii];
				edges[v_n] = edges[jj], edgelist_pointer[v_n++] = edgelist_pointer[jj];
				if ((tri_cnt[edgelist_pointer[ii]]--) == K - 2) Qe.push(edgelist_pointer[ii]);
				if ((tri_cnt[edgelist_pointer[jj]]--) == K - 2) Qe.push(edgelist_pointer[jj]);
				ii++;
				jj++;
			}
			else if (edges[ii] < edges[jj])
			{
				edges[u_n] = edges[ii], edgelist_pointer[u_n++] = edgelist_pointer[ii];
				ii++;
			}
			else if (edges[ii] > edges[jj])
			{
				edges[v_n] = edges[jj], edgelist_pointer[v_n++] = edgelist_pointer[jj];
				jj++;
			}

		}

		while (ii < pend[u])
		{
			if (!deleted_edge[edgelist_pointer[ii]])
				edges[u_n] = edges[ii], edgelist_pointer[u_n++] = edgelist_pointer[ii];
			ii++;
		}

		while (jj < pend[v])
		{
			if (!deleted_edge[edgelist_pointer[jj]])
				edges[v_n] = edges[jj], edgelist_pointer[v_n++] = edgelist_pointer[jj];
			jj++;
		}
		pend[u] = u_n;
		pend[v] = v_n;


		//printf("degree u %d\n", pend[u] - pstart[u]);
		//printf("degree v %d\n", pend[v] - pstart[v]);
		assert(degree[u] == pend[u] - pstart[u]);
		assert(degree[v] == pend[v] - pstart[v]);

		if (degree[u] == 0) deleted_vertex[u] = 1;
		if (degree[v] == 0) deleted_vertex[v] = 1;
	}

	bool checkempty(ui* deleted_vertex)
	{
		for (ui i = 0; i < n; i++)
			if (!deleted_vertex[i])
				return false;
		return true;
	}

	void clean(ui flag)
	{
		//设定flag是为了防止重复入队
		if (!flag)
		{
			for (ui i = 0; i < (m >> 1); i++)
			{
				if (!deleted_edge[i] && tri_cnt[i] < K - 2)
					Qe.push(i);
			}
		}

		while (!Qe.empty())
		{
		ui e = Qe.front();
		Qe.pop();
		delfrD_edge(e);
		}

		for (ui i = 0; i < n; i++)
			if (degree[i] == 0) deleted_vertex[i] = 1;

		//ui remain_n = 0, remain_m = 0;
		//for (ui i = 0; i < n; i++)
		//	if (!deleted_vertex[i]) remain_n++;
		//for (ui i = 0; i < (m >> 1); i++)
		//	if (!deleted_edge[i]) remain_m++;
		//printf("预处理后，还剩%d个节点和%d条边\n", remain_n, remain_m);
	}

	void heuristic()
	{
		ui* temp_deleted_vertex = new ui[n];
		memset(temp_deleted_vertex, 0, sizeof(ui) * n);
		//删除孤立节点
		for (ui i = 0; i < n; i++)
			if (degree[i] == 0) deleted_vertex[i] = 1;

		clean(0);
		while (!checkempty(deleted_vertex))
		{
			//store graph
			for (ui i = 0; i < n; i++)
				temp_deleted_vertex[i] = deleted_vertex[i];

			int min_eid = -1, min_tri = n;

			for (ui i = 0; i < (m >> 1); i++)
			{
				if (!deleted_edge[i] && tri_cnt[i] < min_tri)
				{
					min_eid = i;
					min_tri = tri_cnt[i];
				}
			}

			delfrD_edge(min_eid);
			clean(1);

			//删除孤立节点
			for (ui i = 0; i < n; i++)
				if (degree[i] == 0) deleted_vertex[i] = 1;
		}

		ui cnt = 0;

		for (ui i = 0; i < n; i++)
			if (!temp_deleted_vertex[i]) cnt++;

		LB = cnt;

		printf("min_k_truss的下界 LB = %u\n", LB);

	}

	void heuristic_query(clock_t start, bool& timeflag, int q)
	{
		//删除孤立节点
		for (ui i = 0; i < n; i++)
			if (degree[i] == 0) deleted_vertex[i] = 1;

		clean(0);

		if (checkempty(deleted_vertex))
		{
			printf("图删完了");
			return;
		}

		if (deleted_vertex[q])
		{
			printf("q has been deleted\n");
			return;
		}


		svex.clear();
		sedge.clear();
		common_neighbors.clear();

		inS_vex = new ui[n];
		inS_edge = new ui[m >> 1];
		vex_tri_cnt = new ui[n];
		tri_cnt_inS = new ui[m >> 1];
		inS_degree = new ui[n];
		visit = new ui[n];
		memset(inS_vex, 0, sizeof(ui) * n);
		memset(inS_edge, 0, sizeof(ui) * (m >> 1));
		memset(vex_tri_cnt, 0, sizeof(ui) * n);
		memset(visit, 0, sizeof(ui) * n);
		memset(tri_cnt_inS, 0, sizeof(ui) * (m >> 1));
		memset(inS_degree, 0, sizeof(ui) * (n));


		int max_tri = -1, max_eid = -1;

		for (ept i = pstart[q]; i < pend[q]; i++)
		{
			ui eid = edgelist_pointer[i];
			if (!deleted_edge[eid] && (int)(tri_cnt[eid]) > max_tri)
			{
			   max_tri = tri_cnt[eid];
			   max_eid = eid;
            }
		}

		/*for (ui i = 0; i < (m >> 1); i++)
		{
			if ((int)(tri_cnt[i]) > max_tri)
			{
				max_tri = tri_cnt[i];
				max_eid = i;
			}
		}*/

		if (max_tri != -1)
		{
			ui u = edge_list[max_eid << 1], v = edge_list[(max_eid << 1) + 1];
			svex.push_back(u);
			svex.push_back(v);
			sedge.push_back(max_eid);
			inS_vex[u] = inS_vex[v] = 1;
			inS_edge[max_eid] = 1;
			tri_cnt_inS[max_eid] = 0;
			inS_degree[u] = inS_degree[v] = 1;
			push_common_neighbour(max_eid);
		}


		while (1)
		{
			clock_t finish;
			finish = clock();
			double time = (double)(finish - start) / CLOCKS_PER_SEC;
			if ((double)(finish - start) / CLOCKS_PER_SEC > 1800)
			{
				timeflag = false;
				return;
			}
			for (auto v : common_neighbors)
			{
				ui cnt = 0;
				for (auto e : sedge)
					if (check_common_neighbour(v, e)) cnt++;
				vex_tri_cnt[v] = cnt;
			}

			sort(common_neighbors.begin(), common_neighbors.end(), compare);

			ui u = common_neighbors.back();
			common_neighbors.pop_back();
			if (inS_vex[u]) continue;
			svex.push_back(u);
			//inS_degree[u] = 0;
			inS_vex[u] = 1;


			//加入u在S中的相关联的边，并计算支持度

			for (ui j = pstart[u]; j < pend[u]; j++)
			{
				ui e = edgelist_pointer[j];
				if (!deleted_edge[e] && inS_vex[edges[j]])
				{
					inS_degree[u]++;
					inS_degree[edges[j]]++;
					sedge.push_back(e);
					inS_edge[e] = 1;
					for (auto v : svex)
					{
						if (check_common_neighbour(v, e))
							tri_cnt_inS[e]++;
					}
					push_common_neighbour(e);
				}
			}

			//加入u之后，增加S中与之关联的边的支持度
			for (auto e : sedge)
			{
				if (check_common_neighbour(u, e))
					tri_cnt_inS[e]++;
			}
			//printf("111\n");
			//保存S，看看能否通过删边来得到一个ktruss

			vector<ui> temp_svex, temp_sedge;
			ui* temp_inS_degree = new ui[n];
			ui* temp_inS_vex = new ui[n];
			ui* temp_inS_edge = new ui[m >> 1];
			ui* temp_tri_cnt_inS = new ui[m >> 1];

			for (ui i = 0; i < n; i++)
			{
				temp_inS_degree[i] = inS_degree[i];
				temp_inS_vex[i] = inS_vex[i];
			}

			for (ui i = 0; i < (m >> 1); i++)
			{
				temp_inS_edge[i] = inS_edge[i];
				temp_tri_cnt_inS[i] = tri_cnt_inS[i];
			}


			temp_svex = svex, temp_sedge = sedge;
			queue<ui> del_queue;

			//for (auto e :temp_sedge)
			//{
			//	printf("%d\n", temp_tri_cnt_inS[e]);
			//}

			for (auto e : temp_sedge)
			{
				if (temp_tri_cnt_inS[e] < K - 2)
				{
					del_queue.push(e);
				}
			}

			while (!del_queue.empty())
			{
				ui e = del_queue.front();
				del_queue.pop();
				delfrD_edge_in_heu(e, temp_inS_edge, temp_inS_vex, temp_inS_degree, temp_tri_cnt_inS, del_queue);
			}

			/*		for (auto e : temp_sedge)
					{
						printf("%d\n", temp_tri_cnt_inS[e]);
					}*/

			for (auto v : temp_svex)
				if (temp_inS_vex[v] && temp_inS_degree[v] == 0) temp_inS_vex[v] = 0;

			//for (auto e : sedge)
			//{
			//	if (temp_inS_edge[e])
			//	{
			//		int u = edge_list[e << 1], v = edge_list[(e << 1) + 1];
			//		printf("%d %d\n", temp_inS_degree[u], temp_inS_degree[v]);
			//	}
			//}

			ui cnt = 0;
			for (auto v : temp_svex)
				if (temp_inS_vex[v]) cnt++;
			if (cnt != 0)
			{
				LB = cnt;
				printf("LB = %d\n", LB);
				return;
			}



			//if (exam_truss(sedge))
			//{
			//	for (auto v : svex) printf("%d ", v);
			//	printf("\n");
			//	LB = svex.size();
			//	printf("LB = %d\n", LB);
			//	return;
			//}
		}
	}

	void store_graph()
	{
		temp_m = m;
		temp_n = n;
		temp_pstart = new ept[n];
		temp_pend = new ept[n];
		temp_edgelist_pointer = new ui[m];
		temp_edge_list = new ui[m];
		temp_tri_cnt = new ui[m >> 1];
		temp_degree = new ui[n];
		temp_edges = new ui[m];

		temp_deleted_edge = new ui[m >> 1];
		temp_deleted_vertex = new ui[n];

		for (ui i = 0; i < temp_n; i++)
		{
			temp_pstart[i] = pstart[i];
			temp_pend[i] = pend[i];
			temp_degree[i] = degree[i];
			temp_deleted_vertex[i] = deleted_vertex[i];
		}

		for (ept i = 0; i < temp_m; i++)
		{
			temp_edges[i] = edges[i];
			temp_edge_list[i] = edge_list[i];
			temp_edgelist_pointer[i] = edgelist_pointer[i];

			if (i < (temp_m >> 1))
			{
				temp_tri_cnt[i] = tri_cnt[i];
				temp_deleted_edge[i] = deleted_edge[i];
			}
		}
	}

	void recover_graph()
	{
		m = temp_m;
		n = temp_n;

		for (ui i = 0; i < temp_n; i++)
		{
			pstart[i] = temp_pstart[i];
			pend[i] = temp_pend[i];
			degree[i] = temp_degree[i];
			deleted_vertex[i] = temp_deleted_vertex[i];
		}

		for (ept i = 0; i < temp_m; i++)
		{
			edges[i] = temp_edges[i];
			edge_list[i] = temp_edge_list[i];
			edgelist_pointer[i] = temp_edgelist_pointer[i];

			if (i < (temp_m >> 1))
			{
				tri_cnt[i] = temp_tri_cnt[i];
				deleted_edge[i] = temp_deleted_edge[i];
			}
		}
	}

	void rebulid_graph(ui* tri_cnt, ui* edge_list, ui* edgelist_pointer)
	{
		ept* pstart1 = new ept[n];
		ept* pend1 = new ept[n];

		ui n_cnt = 0;
		ept m_cnt = 0;
		newid = new ui[n];
		old_id = new ui[n];
		exist_oldid = new ui[n];
		memset(newid, 0, sizeof(ui) * n);
		memset(old_id, 0, sizeof(ui) * n);
		memset(exist_oldid, 1, sizeof(ui) * n);

		for (ui i = 0; i < n; i++)
		{
			if (pend[i] - pstart[i] > 0)
			{
				newid[i] = n_cnt;
				old_id[n_cnt] = i;
				pstart1[n_cnt] = pstart[i];
				pend1[n_cnt] = pend[i];
				n_cnt++;
				m_cnt += (pend[i] - pstart[i]);
			}
			else exist_oldid[i] = 0;
		}

		for (ui i = 0; i < (m >> 1); i++)
			if (deleted_edge[i]) tri_cnt[i] = 0;


		if (K > 2)
		{
			ui cnt = 0;
			for (ept idx = 0; idx < (m >> 1); idx++)
			{
				if (tri_cnt[idx] > 0)
				{
					tri_cnt[cnt] = tri_cnt[idx];
					ui u = edge_list[idx << 1], v = edge_list[(idx << 1) + 1];
					edge_list[cnt << 1] = newid[u];
					edge_list[(cnt << 1) + 1] = newid[v];

					for (ept i = pstart[u]; i < pend[u]; i++)
					{
						if (edgelist_pointer[i] == idx)
						{
							edgelist_pointer[i] = cnt;
							break;
						}
					}

					for (ept i = pstart[v]; i < pend[v]; i++)
					{
						if (edgelist_pointer[i] == idx)
						{
							edgelist_pointer[i] = cnt;
							break;
						}
					}
					cnt++;
				}
			}
		}

		for (ui i = 0; i < n; i++)
		{
			for (ept j = pstart[i]; j < pend[i]; j++)
				edges[j] = newid[edges[j]];
		}

		n = n_cnt;
		m = m_cnt;
		pstart = new ept[n + 1];
		pend = new ept[n];
		for (ui i = 0; i < n; i++)
		{
			pstart[i] = pstart1[i];
			pend[i] = pend1[i];
		}

		ept count = 0;
		ui t_degree = 0;
		for (ui i = 0; i < n; i++)
		{
			//t_degree = pend[i] - pstart[i];
			//printf("degree %d = %d ", i, t_degree);
			ui t_n = pstart[i];
			pstart[i] = count;
			for (ept j = t_n; j < pend[i]; j++)
			{
				edges[count] = edges[j];
				edgelist_pointer[count++] = edgelist_pointer[j];
			}
			pend[i] = count;
		}
		cout << endl;
		printf("In the remaining graph: n = %d, m = %d\n", n_cnt, m / 2);
	

	}

};



#endif