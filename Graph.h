#ifndef _GRAPH_H_
#define _GRAPH_H_

#include"Utility.h"
#include "Timer.h"
#include"BnB1.h"
#include"ktb_enum.h"
#include"BNB_GLOBAL.h"
#include"BnB_iterall.h"
#include"BNB-Topc.h"
#include"ktb_enum_topc.h"
#include"BNB-plex.h"
#include"BK-eplex.h"


class Graph {
private:
	ui n;
	ept m;
	ui K;
	ui start_from_zero;
	string dir;

	std::vector<ui> ktruss_bar;

	ept* pstart;
	ept* pend; 
	ept* pend_buf;
	ui* edges;

	ui* edgelist_pointer; 
	ui* edge_list; 
	ui* tri_cnt; 
	ui* old_id;
	ui* old_old_id;
	ui* deleted_oldid;
	ui* vertex_support;


	ui* s_degree;
	ept* s_pstart;
	ept* s_pend;
	ui* s_edges;
	ui* s_peel_sequence;
	ui* s_core;
	char* s_vis;

	ui* s_edgelist_pointer;
	ui* s_tri_cnt;
	ui* s_edge_list;
	ui* s_active_edgelist;
	char* s_deleted;

public:
	Graph(const char* filename, const int _K)
	{

		K = _K;
		dir = string(filename);

		pstart = nullptr;
		pend = pend_buf = nullptr;
		edges = nullptr;

		edgelist_pointer = NULL;
		edge_list = NULL;
		tri_cnt = NULL;

		s_edgelist_pointer = NULL;
		s_tri_cnt = s_edge_list = NULL;
		s_active_edgelist = NULL;
		s_deleted = NULL;

		s_degree = s_edges = NULL;
		s_pstart = s_pend = NULL;
		s_peel_sequence = s_core = NULL;
		s_vis = NULL;
	}

	~Graph()
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
		if (s_degree != NULL) {
			delete[] s_degree;
			s_degree = NULL;
		}
		if (s_pstart != NULL) {
			delete[] s_pstart;
			s_pstart = NULL;
		}
		if (s_pend != NULL) {
			delete[] s_pend;
			s_pend = NULL;
		}
		if (s_edges != NULL) {
			delete[] s_edges;
			s_edges = NULL;
		}
		if (s_peel_sequence != NULL) {
			delete[] s_peel_sequence;
			s_peel_sequence = NULL;
		}
		if (s_core != NULL) {
			delete[] s_core;
			s_core = NULL;
		}
		if (s_vis != NULL) {
			delete[] s_vis;
			s_vis = NULL;
		}
		if (s_edgelist_pointer != NULL) {
			delete[] s_edgelist_pointer;
			s_edgelist_pointer = NULL;
		}
		if (s_active_edgelist != NULL) {
			delete[] s_active_edgelist;
			s_active_edgelist = NULL;
		}
		if (s_deleted != NULL) {
			delete[] s_deleted;
			s_deleted = NULL;
		}
	}


	void read_graph_ratio(double ratio) {

		FILE* f = Utility::open_file((dir).c_str(), "r");
		fscanf(f, "%u%lu%u", &n, &m, &start_from_zero);

		printf("n = %d, m = %d, K = %d\n", n, m, K);
		ui ub_n = (ui)n * ratio;
		ept ub_e = 0;
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
			if (a <= ub_n && b <= ub_n)
			{
				vp.pb(mp(a, b));
				vp.pb(mp(b, a));
				ub_e++;
			}
		}
		n = ub_n, m = ub_e * 2;

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

		printf("After extracting %lf vertices, n = %d, m = %d, K = %d\n", ratio, n, m/2, K);

		fclose(f);

		printf("Finish reading graph!\n");
	}

	//read graph
	void read_graph()
	{
		FILE* f = Utility::open_file((dir).c_str(), "r");
		fscanf(f, "%u%lu%u", &n, &m, &start_from_zero);

		printf("n = %d, m = %d, K = %d\n", n, m, K);
		deleted_oldid = new ui[n];
		memset(deleted_oldid, 0, sizeof(int) * n);
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

	void test_truss()
	{
		read_graph();
		edgelist_pointer = new ui[m];
		pend = new ept[n];
		memset(edgelist_pointer, 0, m * sizeof(ui));
		oriented_triangle_counting(n, m, pstart, pend, edges, edgelist_pointer);



		edge_list = new ui[m];
		tri_cnt = new ui[m / 2];
		ept* pend_buf = new ept[n];

		reorganize_oriented_graph(n, tri_cnt, edge_list, pstart, pend, pend_buf, edges, edgelist_pointer);

		cout << "m = " << m / 2 << endl;
		for (int i = 0; i < m / 2; i++)
		{
			if (tri_cnt[i] < K - 2)
			{
				cout << "!!!!!!!!!!wrong!!!!!!!!!!!" << endl;
				return;
			}
			//cout << tri_cnt[i] << endl;
		}
		cout << "truss!!!" << endl;
		return;
	}

	//max_k_truss preprocess
	void max_truss()
	{
		edgelist_pointer = new ui[m];
		pend = new ept[n];
		memset(edgelist_pointer, 0, m * sizeof(ui));
		oriented_triangle_counting(n, m, pstart, pend, edges, edgelist_pointer);

		

		edge_list = new ui[m];
		tri_cnt = new ui[m / 2];
		ept* pend_buf = new ept[n];

		reorganize_oriented_graph(n, tri_cnt, edge_list, pstart, pend, pend_buf, edges, edgelist_pointer);


		for (ui i = 0; i < n; i++) pend[i] = pstart[i + 1];


		queue<ept> Qe;
		bool* inq = new bool[m];
		bool* deleted = new bool[m];
		memset(inq, false, m);
		memset(deleted, false, m);



		for (ui i = 0; i < n; i++)
		{
			for (ept j = pstart[i]; j < pend[i]; j++)
			{
				ept idx = edgelist_pointer[j];
				ui cnt = tri_cnt[idx];
				if (cnt < K - 2 && !inq[idx])
				{
					//cout <<cnt << endl;
					inq[idx] = true;
					Qe.push(idx);
				}
			}
		}


		while (!Qe.empty())
		{
			ept idx = Qe.front();
			Qe.pop();
			deleted[idx] = true;
			tri_cnt[idx] = 0;
			ui u = edge_list[idx << 1], v = edge_list[(idx << 1) + 1];
			//printf("remove %d-%d\n", u, v);
			ept ii = pstart[u], jj = pstart[v];
			ept u_n = pstart[u], v_n = pstart[v];

			while (true)
			{
				while (ii < pend[u] && deleted[edgelist_pointer[ii]]) ++ii;
				while (jj < pend[v] && deleted[edgelist_pointer[jj]]) ++jj;
				if (ii >= pend[u] || jj >= pend[v]) break;

				if (edges[ii] == edges[jj])
				{
					edges[u_n] = edges[ii], edgelist_pointer[u_n++] = edgelist_pointer[ii];
					edges[v_n] = edges[jj], edgelist_pointer[v_n++] = edgelist_pointer[jj];

					if ((--tri_cnt[edgelist_pointer[ii]]) < (K - 2) && !inq[edgelist_pointer[ii]])
					{
						inq[edgelist_pointer[ii]] = true;
						Qe.push(edgelist_pointer[ii]);
					}

					if ((--tri_cnt[edgelist_pointer[jj]]) < (K - 2) && !inq[edgelist_pointer[jj]])
					{
						inq[edgelist_pointer[jj]] = true;
						Qe.push(edgelist_pointer[jj]);
					}

					++ii;
					++jj;
				}
				else if (edges[ii] < edges[jj])
				{
					edges[u_n] = edges[ii], edgelist_pointer[u_n++] = edgelist_pointer[ii];
					++ii;
				}
				else if (edges[jj] < edges[ii])
				{
					edges[v_n] = edges[jj], edgelist_pointer[v_n++] = edgelist_pointer[jj];
					++jj;
				}
			}

			while (ii < pend[u])
			{
				if (!deleted[edgelist_pointer[ii]])
				{
					edges[u_n] = edges[ii];
					edgelist_pointer[u_n++] = edgelist_pointer[ii];

				}
				++ii;
			}

			while (jj < pend[v])
			{
				if (!deleted[edgelist_pointer[jj]])
				{
					edges[v_n] = edges[jj];
					edgelist_pointer[v_n++] = edgelist_pointer[jj];

				}
				++jj;
			}
			pend[u] = u_n; pend[v] = v_n;
		}


		ui cnt = 0;
		for (ui i = 0; i < n; i++)
		{
			if (pend[i] - pstart[i] > 0)
				cnt++;

		}
		cout << "remaining vertices: " << cnt << endl;

		rebulid_graph(tri_cnt, edge_list, edgelist_pointer);

	}



	ept core_truss_co_pruning(ui& deleted_nodes, queue<ui>& Qv, ui d_threshold, queue<ept>& Qe, ui t_threshold, ui* tri_cnt, ui* active_edgelist, ui& active_edgelist_n, ui* edge_list, ui* edgelist_pointer, bool* deleted, ui* degree, ept* pstart, ept* pend, ui* edges, bool* exists)
	{
		ept new_n = 0;
		for (ept i = 0; i < active_edgelist_n; i++)
		{
			if (!deleted[active_edgelist[i]])
			{
				if (tri_cnt[active_edgelist[i]] < t_threshold) Qe.push(active_edgelist[i]);
				else active_edgelist[new_n++] = active_edgelist[i];
			}
		}
		active_edgelist_n = new_n;

		ui* outgraph = new ui[n];
		memset(outgraph, 0, sizeof(ui) * n);

		ui deleted_n = 0;

		while (!Qe.empty() || !Qv.empty())
		{
			if (Qe.empty())
			{
				ui v = Qv.front();
				Qv.pop();
				if (!outgraph[v]) deleted_nodes++;
				outgraph[v] = 1;
				ept v_n = pstart[v];

				for (ept i = pstart[v]; i < pend[v]; i++)
				{
					if (!deleted[edgelist_pointer[i]])
					{
						edgelist_pointer[v_n] = edgelist_pointer[i];
						edges[v_n++] = edges[i];
						exists[edges[i]] = 1;
					}
				}

				pend[v] = v_n;

				for (ept i = pstart[v]; i < pend[v]; i++) deleted[edgelist_pointer[i]] = true;
				deleted_n += pend[v] - pstart[v];
				degree[v] = 0;

				for (ept i = pstart[v]; i < pend[v]; i++)
				{
					ui u = edges[i];
					ept u_n = pstart[u];
					for (ept j = pstart[u]; j < pend[u]; j++)
					{
						if (!deleted[edgelist_pointer[j]])
						{
							edges[u_n] = edges[j];
							edgelist_pointer[u_n++] = edgelist_pointer[j];
							if (edges[j] > u && exists[edges[j]])
								if ((tri_cnt[edgelist_pointer[j]]--) <= t_threshold) Qe.push(edgelist_pointer[j]);
						}
					}
					pend[u] = u_n;

					if ((degree[u]--) <= d_threshold) Qv.push(u);
				}
				for (ept k = pstart[v]; k < pend[v]; k++) exists[edges[k]] = 0;
			}

			while (!Qe.empty())
			{
				ept e = Qe.front();
				Qe.pop();
				ui u = edge_list[e << 1], v = edge_list[(e << 1) + 1];
				//printf("u=%d v=%d %u %u\n", u, v, degree[u] - 1, degree[v] - 1);
				if ((degree[u]--) <= d_threshold) Qv.push(u);
				if ((degree[v]--) <= d_threshold) Qv.push(v);



				deleted[e] = true;
				deleted_n++;

				ept u_n = pstart[u], v_n = pstart[v];
				ept ii = pstart[u], jj = pstart[v];

				while (true)
				{
					while (ii < pend[u] && !deleted[edgelist_pointer[ii]]) ii++;
					while (jj < pend[v] && !deleted[edgelist_pointer[jj]]) jj++;

					if (ii >= pend[u] || jj >= pend[u]) break;

					if (ii == jj)
					{
						edges[u_n] = edges[ii], edgelist_pointer[u_n++] = edgelist_pointer[ii];
						edges[v_n] = edges[jj], edgelist_pointer[v_n++] = edgelist_pointer[jj];
						if ((tri_cnt[edgelist_pointer[ii]]--) <= t_threshold) Qe.push(edgelist_pointer[ii]);
						if ((tri_cnt[edgelist_pointer[jj]]--) <= t_threshold) Qe.push(edgelist_pointer[jj]);

						ii++;
						jj++;
					}
					else if (ii < jj)
					{
						edges[u_n] = edges[ii], edgelist_pointer[u_n++] = edgelist_pointer[ii];
						ii++;
					}
					else if (jj < ii)
					{
						edges[v_n] = edges[jj], edgelist_pointer[v_n++] = edgelist_pointer[jj];
						jj++;
					}
				}

				while (ii < pend[u])
				{
					if (!deleted[edgelist_pointer[ii]])
						edges[u_n] = edges[ii], edgelist_pointer[u_n++] = edgelist_pointer[ii];
					ii++;
				}
				pend[u] = u_n;

				while (jj < pend[v])
				{
					if (!deleted[edgelist_pointer[jj]])
						edges[v_n] = edges[jj], edgelist_pointer[v_n++] = edgelist_pointer[jj];
					jj++;
				}
				pend[v] = v_n;
			}
		}

		return deleted_n;
	}


	ept CTCP(bool* inq, ui& deleted_nodes, queue<ui>& Qv, ui d_threshold, queue<ept>& Qe, ui t_threshold, ui* tri_cnt, ui* active_edgelist, ui& active_edgelist_n, ui* edge_list, ui* edgelist_pointer, bool* deleted, ui* degree, ept* pstart, ept* pend, ui* edges, bool* exists)
	{
		for (ept i = 0; i < m / 2; i++)
		{
			if (!deleted[i] && tri_cnt[i] < t_threshold) Qe.push(i);
		}

		if (Qv.empty()) deleted_nodes = 0;
		else deleted_nodes = 1;

		ui deleted_edges = 0;
		peeling(deleted_edges, Qv, Qe, t_threshold, tri_cnt, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);

		bool flag = false;

		while (true)
		{
			flag = false;
			for (ui i = 0; i < n; i++)
			{
				if (degree[i] < d_threshold && !inq[i])
				{
					Qv.push(i);
					inq[i] = true;
					flag = true;
					peeling(deleted_edges, Qv, Qe, t_threshold, tri_cnt, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
					/*for (ui i = 0; i < n; i++)
						printf("degree %d = %d\n", i + 1, degree[i]);
					printf("\n");*/
					deleted_nodes++;
				}
			}
			if (!flag) break;
		}

		return deleted_edges;
	}
	void peeling(ui& deleted_edges, queue<ui>& Qv, queue<ept>& Qe, ui t_threshold, ui* tri_cnt, ui* edge_list, ui* edgelist_pointer, bool* deleted, ui* degree, ept* pstart, ept* pend, ui* edges, bool* exists)
	{
		while (!Qv.empty() || !Qe.empty())
		{
			while (!Qe.empty())
			{
				ui e = Qe.front();
				Qe.pop();
				deleted_edges++;
				deleted[e] = true;
				ui u = edge_list[e << 1], v = edge_list[(e << 1) + 1];
				//printf("remove %d - %d\n", u+1, v+1);
				if (degree[u] <= 0) printf("%d\n", degree[u]);
				assert(degree[u] > 0);
				assert(degree[v] > 0);

				degree[u]--;
				degree[v]--;

				ept ii = pstart[u], jj = pstart[v];
				ept u_n = pstart[u], v_n = pstart[v];

				while (true)
				{
					while (ii < pend[u] && deleted[edgelist_pointer[ii]]) ii++;
					while (jj < pend[v] && deleted[edgelist_pointer[jj]]) jj++;
					if (ii >= pend[u] || jj >= pend[v]) break;

					if (edges[ii] == edges[jj])
					{
						edges[u_n] = edges[ii], edgelist_pointer[u_n++] = edgelist_pointer[ii];
						edges[v_n] = edges[jj], edgelist_pointer[v_n++] = edgelist_pointer[jj];
						if ((tri_cnt[edgelist_pointer[ii]]--) == t_threshold) Qe.push(edgelist_pointer[ii]);
						if ((tri_cnt[edgelist_pointer[jj]]--) == t_threshold) Qe.push(edgelist_pointer[jj]);
						ii++;
						jj++;
					}
					else if (edges[ii] < edges[jj])
					{
						edges[u_n] = edges[ii], edgelist_pointer[u_n++] = edgelist_pointer[ii];
						ii++;
					}
					else
					{
						edges[v_n] = edges[jj], edgelist_pointer[v_n++] = edgelist_pointer[jj];
						jj++;
					}
				}

				while (ii < pend[u])
				{
					if (!deleted[edgelist_pointer[ii]])
						edges[u_n] = edges[ii], edgelist_pointer[u_n++] = edgelist_pointer[ii];
					ii++;
				}

				while (jj < pend[v])
				{
					if (!deleted[edgelist_pointer[jj]])
						edges[v_n] = edges[jj], edgelist_pointer[v_n++] = edgelist_pointer[jj];
					jj++;
				}


				pend[u] = u_n;
				pend[v] = v_n;

				//printf("%d: %d %d\n", u + 1, pstart[u], pend[u]);
			}

			if (!Qv.empty())
			{
				ui u = Qv.front();
				Qv.pop();
				//printf("\nremove %d\n", u + 1);
				ept u_n = pstart[u];

				//test
			/*	for (ept k = pstart[u]; k < pend[u]; k++) printf("%d ", edges[k]);
				printf("\n");*/

				for (ept k = pstart[u]; k < pend[u]; k++) if (!deleted[edgelist_pointer[k]]) {
					edges[u_n] = edges[k]; edgelist_pointer[u_n++] = edgelist_pointer[k];
					exists[edges[k]] = 1;
				}
				pend[u] = u_n;

				for (ept k = pstart[u]; k < pend[u]; k++) deleted[edgelist_pointer[k]] = 1;
				deleted_edges += pend[u] - pstart[u];
				degree[u] = 0;

				for (ept k = pstart[u]; k < pend[u]; k++) {
					ui v = edges[k];
					ept v_n = pstart[v];
					for (ept x = pstart[v]; x < pend[v]; x++) if (!deleted[edgelist_pointer[x]]) {
						edges[v_n] = edges[x]; edgelist_pointer[v_n++] = edgelist_pointer[x];
						if (edges[x] > v && exists[edges[x]]) {
							if ((tri_cnt[edgelist_pointer[x]]--) == t_threshold) Qe.push(edgelist_pointer[x]);
						}
					}
					pend[v] = v_n;
					degree[v]--;
				}
				for (ept k = pstart[u]; k < pend[u]; k++) exists[edges[k]] = 0;
			}
		}
	}

	void rebulid_graph(ui* tri_cnt, ui* edge_list, ui* edgelist_pointer)
	{
		ept* pstart1 = new ept[n];
		ept* pend1 = new ept[n];

		ui n_cnt = 0;
		ept m_cnt = 0;
		ui* newid = new ui[n];
		old_old_id = new ui[n];

		for (ui i = 0; i < n; i++)
		{
			if (pend[i] - pstart[i] > 0)
			{
				newid[i] = n_cnt;
				pstart1[n_cnt] = pstart[i];
				old_old_id[n_cnt] = i;
				pend1[n_cnt] = pend[i];
				n_cnt++;
				m_cnt += (pend[i] - pstart[i]);
			}
			else deleted_oldid[i] = 1;
		}


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
		pstart = new ept[n];
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

	int find_min_k_truss_query(clock_t start, bool& timeflag, string method, int choosed_vertex)
	{
		read_graph();
		max_truss();

		if (deleted_oldid[choosed_vertex])
		{
			printf("Query vertex has been deleted by the maximum ktruss\n");
			return -1;
		}
		const char* filename;
		filename = "/home/zhangqifan/min_k_truss/dataset/ktruss1.txt";
		ofstream outfile;
		outfile.open(filename);
		for (int i = 0; i < (m >> 1); i++)
		{
			int u = edge_list[i << 1], v = edge_list[(i << 1) + 1];
			if (start_from_zero)
				outfile << old_old_id[u] << " " << old_old_id[v] << endl;
			else
				outfile << old_old_id[u] + 1 << " " << old_old_id[v] + 1 << endl;
		}
		outfile.close();



		ui temp_n = n;
		ept temp_m = m;

		ept* temp_pstart = new ept[n];
		ept* temp_pend = new ept[n];
		ui* temp_edges = new ui[m];
		ui* temp_edgelist_pointer = new ui[m];
		ui* temp_tri_cnt = new ui[m >> 1];
		ui* temp_edge_list = new ui[m];

		memset(temp_pstart, 0, sizeof(ept) * n);
		memset(temp_pstart, 0, sizeof(ept) * n);
		memset(temp_edges, 0, sizeof(ui) * m);
		memset(temp_edgelist_pointer, 0, sizeof(ui) * m);
		memset(temp_tri_cnt, 0, sizeof(ui) * (m >> 1));
		memset(temp_edge_list, 0, sizeof(ui) * m);


		//store graph

		for (ui i = 0; i < n; i++)
		{
			temp_pstart[i] = pstart[i];
			temp_pend[i] = pend[i];
			for (ept j = pstart[i]; j < pend[i]; j++)
			{
				ui e = edgelist_pointer[j];
				temp_edges[j] = edges[j];
				temp_edgelist_pointer[j] = edgelist_pointer[j];
				temp_edge_list[e << 1] = edge_list[e << 1];
				temp_edge_list[(e << 1) + 1] = edge_list[(e << 1) + 1];
				//printf("j = %d and %d-%d %d\n",j, i, edges[j], temp_edgelist_pointer[j]);
			}
		}

		for (ui i = 0; i < (m >> 1); i++) temp_tri_cnt[i] = tri_cnt[i];


		bool ans = false;

		for (ui i = 0; i <= n - K && !ans; i++)
		{

			if ((int)(n - K) < 0) break;
			printf("n = %d, m = %d\n", n, m / 2);
			ui S = i;


			ui* active_edgelist = new ui[m >> 1];
			ui active_edgelist_n = m >> 1;
			for (ui i = 0; i < (m >> 1); i++) active_edgelist[i] = i;

			bool* deleted = new bool[m >> 1];
			memset(deleted, false, sizeof(bool) * (m >> 1));
			bool* exists = new bool[n];
			memset(exists, false, sizeof(bool) * n);

			ui* degree = new ui[n];
			for (ui i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];

			queue<ui> Qv;
			queue<ept> Qe;



			ui lb = S + K;
			printf("S = %d and LB = %d\n", S, lb);
			ui deleted_nodes = 0;

			bool* inq = new bool[n];
			memset(inq, false, sizeof(bool) * n);

			ept x = CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);//ֻҪ����lb�ͺ��ˣ�����Ҫ����
			//ept x = CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - 2 * S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);

			printf("deleted %lu edges\n", x);
			printf("deleted %d nodes\n", deleted_nodes);
			m -= 2 * x;
			//printf("After CTCP, n = %u, m = %lu\n\n", n - deleted_nodes, m / 2);

			for (ui i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];


			
			//Choose vertex. Find 2-hop neighbours. Branch and Bound.
			ui* choosed = new ui[n];
			memset(choosed, 0, sizeof(ui) * n);

			int query_v;
			for (int ii = 0; ii < n; ii++)
			{
				if (old_old_id[ii] == choosed_vertex)
				{
					query_v = ii;
					break;
				}
			}

			for (int j = 0; j < n && m; j++) //&&ktruss_bar.size()<=UB
			{

				int min_key = n, min_id = -1;
				min_id = query_v;
				min_key = degree[min_id];

				if (m == 0) break;

				assert(degree[min_id] == min_key);
				ui* ids = new ui[n];
				ui ids_n = 0;
				ui* rid = new ui[n];
				vector<pair<int, int> > vp; vp.reserve(m / 2);
				ui* Q = new ui[n];
				ui* t_degree = new ui[n];
				memset(t_degree, 0, sizeof(ui) * n);
				ui* exist = new ui[n];
				memset(exist, 0, sizeof(ui) * n);


				//�Ҷ����ھ�
				extract_subgraph_and_prune(lb, S + 1, min_id, ids, ids_n, rid, vp, Q, t_degree, exist, pend, deleted, edgelist_pointer);

				if (ids_n != 0)
					printf("\nids_n = %d\n", ids_n);

				//for (ui i = 0; i < ids_n; i++) printf(" %u", ids[i]);

				//for (ui i = 0; i < ids_n; i++) printf("\n%d", rid[ids[i]]);
				//for (ui i = 0; i < vp.size(); i++) printf("\n%d-%d", vp[i].first, vp[i].second);


				if (ids_n >= lb)
				{
					old_id = new ui[ids_n];
					for (int i = 0; i < ids_n; i++)
					{
						old_id[rid[ids[i]]] = old_old_id[ids[i]];
					}

					BnB1* bnb = new BnB1(S, lb);
					bnb->load_graph(rid[min_id], ids_n, rid, vp, old_id);

					bnb->rearrange_graph();   ///��ԭ

					//ɾ
					/*bnb->test();
					exit(1);*/


					bool ret = bnb->bb_s_truss_bar(ids_n, start, timeflag);
					bnb->~BnB1();


					//test!!!!!!!!!!!!!!�ǵ�ɾ
					//ret = true;


					if (ret)
					{
						printf("choose vertex = %d\n", choosed_vertex);
						printf("K = %d and Find a min_k_truss of size: %d\n", K, S + K);
						ans = true;
						return S + K;
						break;
					}
					
					if (!ret)
					{
						if (S == K - 1) return -1;
						//printf("Cannot find query vertex based ans!\n");
						break;
						//return -1;
					}

				}

			}
			


		   //Return to old graph

			n = temp_n;
			m = temp_m;


			for (ui i = 0; i < n; i++)
			{
				pstart[i] = temp_pstart[i];
				pend[i] = temp_pend[i];

				for (ept j = temp_pstart[i]; j < temp_pend[i]; j++)
				{
					ui e = temp_edgelist_pointer[j];
					edges[j] = temp_edges[j];
					edgelist_pointer[j] = temp_edgelist_pointer[j];
					//printf("j = %d and %d-%d %d\n", j, i, edges[j], e);
					//printf("%d %d\n", edge_list[e << 1], temp_edge_list[e << 1]);
					edge_list[e << 1] = temp_edge_list[e << 1];
					edge_list[(e << 1) + 1] = temp_edge_list[(e << 1) + 1];
				}
			}


			for (ui i = 0; i < (temp_m >> 1); i++)
			{
				tri_cnt[i] = temp_tri_cnt[i];
				//edge_list[i << 1] = temp_edge_list[i << 1];
				//edge_list[(i << 1) + 1] = temp_edge_list[(i << 1) + 1];
			}


		}

		if (!ans)
		{
			printf("There is no %d-truss in the graph.\n", K);
			return -1;
		}
	}

	static bool cmp(const pair<ui, ui> &a, const pair<ui, ui> &b)
	{
		return a.second > b.second;
	}


	int find_min_k_truss(clock_t start, bool& timeflag, string method)
	{
		read_graph();
		max_truss();
		ofstream outfile;

		//output maximum k-truss
		outfile.open("/home/zhangqifan/min_k_truss/dataset/ktruss1.txt");
		for (int i = 0; i < (m >> 1); i++)
		{
			int u = edge_list[i << 1], v = edge_list[(i << 1) + 1];
			if (start_from_zero)
				outfile << old_old_id[u] << " " << old_old_id[v] << endl;
			else
				outfile << old_old_id[u] + 1 << " " << old_old_id[v] + 1 << endl;
		}
		outfile.close();


		ui temp_n = n;
		ept temp_m = m;

		ept* temp_pstart = new ept[n];
		ept* temp_pend = new ept[n];
		ui* temp_edges = new ui[m];
		ui* temp_edgelist_pointer = new ui[m];
		ui* temp_tri_cnt = new ui[m >> 1];
		ui* temp_edge_list = new ui[m];

		memset(temp_pstart, 0, sizeof(ept) * n);
		memset(temp_pstart, 0, sizeof(ept) * n);
		memset(temp_edges, 0, sizeof(ui) * m);
		memset(temp_edgelist_pointer, 0, sizeof(ui) * m);
		memset(temp_tri_cnt, 0, sizeof(ui) * (m >> 1));
		memset(temp_edge_list, 0, sizeof(ui) * m);


		//store graph

		for (ui i = 0; i < n; i++)
		{
			temp_pstart[i] = pstart[i];
			temp_pend[i] = pend[i];
			for (ept j = pstart[i]; j < pend[i]; j++)
			{
				ui e = edgelist_pointer[j];
				temp_edges[j] = edges[j];
				temp_edgelist_pointer[j] = edgelist_pointer[j];
				temp_edge_list[e << 1] = edge_list[e << 1];
				temp_edge_list[(e << 1) + 1] = edge_list[(e << 1) + 1];
				//printf("j = %d and %d-%d %d\n",j, i, edges[j], temp_edgelist_pointer[j]);
			}
		}

		for (ui i = 0; i < (m >> 1); i++) temp_tri_cnt[i] = tri_cnt[i];


		bool ans = false;

		//Enter the iterative framework and start traversing from s = 0

		for (ui i = 0; i <= n - K && !ans; i++)
		{
			if ((int)(n - K) < 0) break;
			printf("n = %d, m = %d\n", n, m / 2);
			ui S = i;
			ui* active_edgelist = new ui[m >> 1];
			ui active_edgelist_n = m >> 1;
			for (ui i = 0; i < (m >> 1); i++) active_edgelist[i] = i;

			bool* deleted = new bool[m >> 1];
			memset(deleted, false, sizeof(bool) * (m >> 1));
			bool* exists = new bool[n];
			memset(exists, false, sizeof(bool) * n);

			ui* degree = new ui[n];
			for (ui i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];

			queue<ui> Qv;
			queue<ept> Qe;

			//The size of s-truss-bar is at least S + K
			ui lb = S + K;
			printf("S = %d and LB = %d\n", S, lb);
			ui deleted_nodes = 0;

			bool* inq = new bool[n];
			memset(inq, false, sizeof(bool) * n);


			if (method == "ENUM")
			{
				vector<pair<int, int>> vp;
				ui* rid = new ui[n];
				for (ui i = 0; i < n; i++) rid[i] = i;
				for (ui j = 0; j < (m >> 1); j++)
				{
					ui x = edge_list[j << 1], y = edge_list[(j << 1) + 1];
					vp.push_back(mp(x, y));
				}


				bool ret;
				ENUM* enum1 = new ENUM(S, lb);
				enum1->load_graph_GLOBAL(n, rid, vp); //enum1->load_graph(rid[min_id], ids_n, rid, vp, ids);

				
				ret = enum1->find_bounded_stb_Global(start, timeflag);
				enum1->~ENUM();
				if (ret)
				{
					printf("K = %d and ENUM_G11 finds a min_k_truss of size: %d\n", K, S + K);
					ans = true;
					return S + K;
					break;
				}
				else continue;
			}

			//preprocessing
			ept x = CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);//ֻҪ����lb�ͺ��ˣ�����Ҫ����
			//ept x = CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - 2 * S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
			m -= 2 * x;

			for (ui i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];

			
			//Choose vertex. Find 2-hop neighbours. Branch and Bound.
			ui* choosed = new ui[n];
			memset(choosed, 0, sizeof(ui) * n);
			for (int j = 0; j < n && m; j++) //&&ktruss_bar.size()<=UB
			{

				int min_key = n, min_id = -1;
				for (int k = 0; k < n; k++)
				{
					if (degree[k] < min_key && !choosed[k])
					{
						min_key = degree[k];
						min_id = k;
					}
				}
				choosed[min_id] = 1;


				//It is impossible to generate a solution larger than the current one��prune
				if (min_key < lb - (S + 1))
				{
					if (degree[min_key] != 0)  //"degree = 0" indicates that the vertex has been deleted
					{
						while (!Qv.empty()) Qv.pop();
						while (!Qe.empty()) Qe.pop();
						Qv.push(min_id);
						deleted_nodes = 0;
						m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
						//m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - 2*S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
					}
					continue;
				}

				if (m == 0) break;

				assert(degree[min_id] == min_key);
				ui* ids = new ui[n];
				ui ids_n = 0;
				ui* rid = new ui[n];
				vector<pair<int, int> > vp; vp.reserve(m / 2);
				ui* Q = new ui[n];
				ui* t_degree = new ui[n];
				memset(t_degree, 0, sizeof(ui) * n);
				ui* exist = new ui[n];
				memset(exist, 0, sizeof(ui) * n);

		

				//Find 2-hop neighbors
				extract_subgraph_and_prune(lb, S + 1, min_id, ids, ids_n, rid, vp, Q, t_degree, exist, pend, deleted, edgelist_pointer);

				if (ids_n != 0)
					printf("\nS = %d and ids_n = %d\n",S, ids_n);


				if (ids_n >= lb)
				{
					old_id = new ui[ids_n];
					for (int i = 0; i < ids_n; i++)
					{
						old_id[rid[ids[i]]] = old_old_id[ids[i]];
					}


					if (method == "IE_ENUM")
					{
						/*for (auto u : vp)
							printf("%d-%d\n", u.first, u.second);*/
						bool ret;
						ENUM* enum1 = new ENUM(S, lb);
						enum1->load_graph_IE(rid[min_id], ids_n, rid, vp, old_id); //enum1->load_graph(rid[min_id], ids_n, rid, vp, ids);

						ret = enum1->find_bounded_stb(start, timeflag);
						enum1->~ENUM();
						if (ret)
						{
							printf("choose vertex = %d\n", old_old_id[min_id]);
							printf("K = %d and Find a min_k_truss of size: %d\n", K, S + K);
							ans = true;
							return S + K;
							break;
						}
					}
					else if (method == "BNB")
					{
						BnB1* bnb = new BnB1(S, lb);
						bnb->load_graph(rid[min_id], ids_n, rid, vp, old_id);

						bnb->rearrange_graph();   ///��ԭ

						//ɾ
						/*bnb->test();
						exit(1);*/


						bool ret = bnb->bb_s_truss_bar(ids_n, start, timeflag);
						bnb->~BnB1();


						if (ret)
						{
							printf("choose vertex = %d\n", old_old_id[min_id]);
							printf("K = %d and Find a min_k_truss of size: %d\n", K, S + K);
							ans = true;
							return S + K;
							break;
						}
					}

					else if (method == "BKEPLEX")
					{
						BKEPLEX* bnbplex = new BKEPLEX(S, lb);
						bnbplex->load_graph(rid[min_id], ids_n, rid, vp, old_id);
						bnbplex->rearrange_graph();

						bool ret = bnbplex->bb_s_eplex(ids_n, start, timeflag, 1);
						bnbplex->~BKEPLEX();

						if (ret)
						{
							printf("choose vertex = %d\n", old_old_id[min_id]);
							printf("K = %d and Find a min_k_truss of size: %d\n", K, S + K);
							ans = true;
							return S + K;
							break;
						}
					}
				}

				while (!Qv.empty()) Qv.pop();
				while (!Qe.empty()) Qe.pop();
				Qv.push(min_id);
				deleted_nodes = 0;
				m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
				//m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - 2*S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);

			}



		   //Return to old graph

			n = temp_n;
			m = temp_m;


			for (ui i = 0; i < n; i++)
			{
				pstart[i] = temp_pstart[i];
				pend[i] = temp_pend[i];

				for (ept j = temp_pstart[i]; j < temp_pend[i]; j++)
				{
					ui e = temp_edgelist_pointer[j];
					edges[j] = temp_edges[j];
					edgelist_pointer[j] = temp_edgelist_pointer[j];
					//printf("j = %d and %d-%d %d\n", j, i, edges[j], e);
					//printf("%d %d\n", edge_list[e << 1], temp_edge_list[e << 1]);
					edge_list[e << 1] = temp_edge_list[e << 1];
					edge_list[(e << 1) + 1] = temp_edge_list[(e << 1) + 1];
				}
			}


			for (ui i = 0; i < (temp_m >> 1); i++)
			{
				tri_cnt[i] = temp_tri_cnt[i];
				//edge_list[i << 1] = temp_edge_list[i << 1];
				//edge_list[(i << 1) + 1] = temp_edge_list[(i << 1) + 1];
			}


		}

		if (!ans)
		{
			printf("There is no %d-truss in the graph.\n", K);
			return -1;
		}


	}

	int find_c_min_k_truss(string method, clock_t start, bool& timeflag, int c)
	{
		int C = 0;
		read_graph();
		max_truss();

		ui temp_n = n;
		ept temp_m = m;

		ept* temp_pstart = new ept[n];
		ept* temp_pend = new ept[n];
		ui* temp_edges = new ui[m];
		ui* temp_edgelist_pointer = new ui[m];
		ui* temp_tri_cnt = new ui[m >> 1];
		ui* temp_edge_list = new ui[m];

		memset(temp_pstart, 0, sizeof(ept) * n);
		memset(temp_pstart, 0, sizeof(ept) * n);
		memset(temp_edges, 0, sizeof(ui) * m);
		memset(temp_edgelist_pointer, 0, sizeof(ui) * m);
		memset(temp_tri_cnt, 0, sizeof(ui) * (m >> 1));
		memset(temp_edge_list, 0, sizeof(ui) * m);


		//store graph

		for (ui i = 0; i < n; i++)
		{
			temp_pstart[i] = pstart[i];
			temp_pend[i] = pend[i];
			for (ept j = pstart[i]; j < pend[i]; j++)
			{
				ui e = edgelist_pointer[j];
				temp_edges[j] = edges[j];
				temp_edgelist_pointer[j] = edgelist_pointer[j];
				temp_edge_list[e << 1] = edge_list[e << 1];
				temp_edge_list[(e << 1) + 1] = edge_list[(e << 1) + 1];
				//printf("j = %d and %d-%d %d\n",j, i, edges[j], temp_edgelist_pointer[j]);
			}
		}

		for (ui i = 0; i < (m >> 1); i++) temp_tri_cnt[i] = tri_cnt[i];


		bool ans = false;
		vector<vector<int>> topc_ans;

		for (ui i = 0; i <= n - K && !ans; i++)
		{
			//Preprocessing
			if ((int)(n - K) < 0) break;
			printf("n = %d, m = %d\n", n, m / 2);
			ui S = i;

			ui* active_edgelist = new ui[m >> 1];
			ui active_edgelist_n = m >> 1;
			for (ui i = 0; i < (m >> 1); i++) active_edgelist[i] = i;

			bool* deleted = new bool[m >> 1];
			memset(deleted, false, sizeof(bool) * (m >> 1));
			bool* exists = new bool[n];
			memset(exists, false, sizeof(bool) * n);

			ui* degree = new ui[n];
			for (ui i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];

			queue<ui> Qv;
			queue<ept> Qe;


			//ui lb = ktruss_bar.size(); ����,lbӦ����ΪS+K
			ui lb = S + K;
			printf("S = %d and LB = %d\n", S, lb);
			ui deleted_nodes = 0;

			bool* inq = new bool[n];
			memset(inq, false, sizeof(bool) * n);

			//ept x = CTCP(inq, deleted_nodes, Qv, lb - S, Qe, lb - 1 - S, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
			//ept x = CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);//ֻҪ����lb�ͺ��ˣ�����Ҫ����

			ept x = CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - 2 * S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);//ֻҪ����lb�ͺ��ˣ�����Ҫ����

			printf("deleted %lu edges\n", x);
			printf("deleted %d nodes\n", deleted_nodes);
			m -= 2 * x;
			printf("After CTCP, n = %u, m = %lu\n\n", n - deleted_nodes, m / 2);


			for (ui i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];


			if (method == "BNBGLOBAL" || S > K - 1)
			{
				ui* ids = new ui[n];
				ui ids_n = 0, rid_n = 0;
				ui* rid = new ui[n];
				vector<pair<int, int> > vp; vp.reserve(m / 2);
				for (ui i = 0; i < n; i++)
				{
					if (degree[i] > 0)
					{
						ids[ids_n] = i;
						rid[ids[ids_n]] = rid_n;
						ids_n++;
						rid_n++;
					}
				}

				for (ui i = 0; i < n; i++)
				{
					if (degree[i] > 0)
					{
						for (ept j = pstart[i]; j < pend[i]; j++)
						{
							if (edges[j] > i && !deleted[edgelist_pointer[j]]) vp.push_back(mp(rid[i], rid[edges[j]]));
						}
					}
				}

				bool ret;
				BnBG* bnbg = new BnBG(S, lb);
				bnbg->load_graph(ids_n, rid, vp, ids);
				bnbg->rearrange_graph();
				ret = bnbg->bb_s_truss_bar(ids_n, start, timeflag);
				bnbg->~BnBG();

				if (ret)
				{
					C++;
					printf("K = %d and Find the %d-th min_k_truss of size: %d\n", K, C, S + K);
					ans = true;
					if (C == c) return S + K;
				}
				continue;
			}

			//Choose vertex. Find 2-hop neighbours. Branch and Bound.
			ui* choosed = new ui[n];
			memset(choosed, 0, sizeof(ui) * n);


			for (int j = 0; j < n && m; j++) //&&ktruss_bar.size()<=UB
			{
				int min_key = n, min_id = -1;
				for (int k = 0; k < n; k++)
				{
					if (degree[k] < min_key && !choosed[k])
					{
						min_key = degree[k];
						min_id = k;
					}
				}
				choosed[min_id] = 1;

				if (min_key < lb - (S + 1))
				{
					if (degree[min_key] != 0) 
					{
						while (!Qv.empty()) Qv.pop();
						while (!Qe.empty()) Qe.pop();
						Qv.push(min_id);
						deleted_nodes = 0;
						m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
						//m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - 2*S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
					}
					continue;
				}

				if (m == 0) break;

				assert(degree[min_id] == min_key);
				ui* ids = new ui[n];
				ui ids_n = 0;
				ui* rid = new ui[n];
				vector<pair<int, int> > vp; vp.reserve(m / 2);
				ui* Q = new ui[n];
				ui* t_degree = new ui[n];
				memset(t_degree, 0, sizeof(ui) * n);
				ui* exist = new ui[n];
				memset(exist, 0, sizeof(ui) * n);

				//�Ҷ����ھ�
				extract_subgraph_and_prune(lb, S + 1, min_id, ids, ids_n, rid, vp, Q, t_degree, exist, pend, deleted, edgelist_pointer);

				if (ids_n != 0)
					printf("\n S= %d and ids_n = %d\n",S, ids_n);


				if (ids_n >= lb)
				{
					old_id = new ui[ids_n];
					for (int i = 0; i < ids_n; i++)
					{
						old_id[rid[ids[i]]] = old_old_id[ids[i]];
					}

					if (method == "IE_ENUM")
					{
						/*for (auto u : vp)
							printf("%d-%d\n", u.first, u.second);*/
						bool ret;
						ENUMC* enum1 = new ENUMC(S, lb);
						enum1->load_graph_IE(rid[min_id], ids_n, rid, vp, old_id, c, topc_ans); //enum1->load_graph(rid[min_id], ids_n, rid, vp, ids);

						ret = enum1->find_bounded_stb(start, timeflag, C);
						enum1->update_topc(topc_ans);
						enum1->~ENUMC();
						if (C >= c)
						{
							ans = true;
							return S + K;
						}
					}
					else if (method == "BNB")
					{
						bool ret;
						BnBC* bnb = new BnBC(S, lb);
						bnb->load_graph(rid[min_id], ids_n, rid, vp, ids, c, old_id, topc_ans);

						bnb->rearrange_graph();   ///��ԭ

						ret = bnb->bb_s_truss_bar(ids_n, start, timeflag, C);
						bnb->update_topc(topc_ans);
						bnb->~BnBC();
						if (C >= c)
						{
							ans = true;
							return S + K;
						}
					}
				}

				while (!Qv.empty()) Qv.pop();
				while (!Qe.empty()) Qe.pop();
				Qv.push(min_id);
				deleted_nodes = 0;
				m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
				//m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - 2*S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);

			}



			//Return to old graph

			n = temp_n;
			m = temp_m;


			for (ui i = 0; i < n; i++)
			{
				pstart[i] = temp_pstart[i];
				pend[i] = temp_pend[i];

				for (ept j = temp_pstart[i]; j < temp_pend[i]; j++)
				{
					ui e = temp_edgelist_pointer[j];
					edges[j] = temp_edges[j];
					edgelist_pointer[j] = temp_edgelist_pointer[j];
					//printf("j = %d and %d-%d %d\n", j, i, edges[j], e);
					//printf("%d %d\n", edge_list[e << 1], temp_edge_list[e << 1]);
					edge_list[e << 1] = temp_edge_list[e << 1];
					edge_list[(e << 1) + 1] = temp_edge_list[(e << 1) + 1];
				}
			}


			for (ui i = 0; i < (temp_m >> 1); i++)
			{
				tri_cnt[i] = temp_tri_cnt[i];
			}


		}

		if (!ans)
		{
			printf("There is no %d-truss in the graph.\n", K);
		}


	}

	int find_min_k_truss_iterall(clock_t start, bool& timeflag, string method)
	{
		int minsize = 100000000;
		read_graph();
		max_truss();

		ui temp_n = n;
		ept temp_m = m;

		ept* temp_pstart = new ept[n];
		ept* temp_pend = new ept[n];
		ui* temp_edges = new ui[m];
		ui* temp_edgelist_pointer = new ui[m];
		ui* temp_tri_cnt = new ui[m >> 1];
		ui* temp_edge_list = new ui[m];

		memset(temp_pstart, 0, sizeof(ept) * n);
		memset(temp_pstart, 0, sizeof(ept) * n);
		memset(temp_edges, 0, sizeof(ui) * m);
		memset(temp_edgelist_pointer, 0, sizeof(ui) * m);
		memset(temp_tri_cnt, 0, sizeof(ui) * (m >> 1));
		memset(temp_edge_list, 0, sizeof(ui) * m);


		//store graph

		for (ui i = 0; i < n; i++)
		{
			temp_pstart[i] = pstart[i];
			temp_pend[i] = pend[i];
			for (ept j = pstart[i]; j < pend[i]; j++)
			{
				ui e = edgelist_pointer[j];
				temp_edges[j] = edges[j];
				temp_edgelist_pointer[j] = edgelist_pointer[j];
				temp_edge_list[e << 1] = edge_list[e << 1];
				temp_edge_list[(e << 1) + 1] = edge_list[(e << 1) + 1];
				//printf("j = %d and %d-%d %d\n",j, i, edges[j], temp_edgelist_pointer[j]);
			}
		}

		for (ui i = 0; i < (m >> 1); i++) temp_tri_cnt[i] = tri_cnt[i];


		bool ans = false;

		for (ui i = 0; i <= n - K; i++)
		{
			clock_t finish;
			finish = clock();
			double time = (double)(finish - start) / CLOCKS_PER_SEC;
			if (time > 10000)
			{
				timeflag = false;
				return -1;
			}
			printf("**********current minsize = %d\n", minsize);
			//Preprocessing
			if ((int)(n - K) < 0) break;
			printf("n = %d, m = %d\n", n, m / 2);
			ui S = i;

			ui* active_edgelist = new ui[m >> 1];
			ui active_edgelist_n = m >> 1;
			for (ui i = 0; i < (m >> 1); i++) active_edgelist[i] = i;

			bool* deleted = new bool[m >> 1];
			memset(deleted, false, sizeof(bool) * (m >> 1));
			bool* exists = new bool[n];
			memset(exists, false, sizeof(bool) * n);

			ui* degree = new ui[n];
			for (ui i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];

			queue<ui> Qv;
			queue<ept> Qe;

			ui lb = S + K;
			printf("S = %d and LB = %d\n", S, lb);
			ui deleted_nodes = 0;

			bool* inq = new bool[n];
			memset(inq, false, sizeof(bool) * n);

			//ept x = CTCP(inq, deleted_nodes, Qv, lb - S, Qe, lb - 1 - S, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
			ept x = CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);//ֻҪ����lb�ͺ��ˣ�����Ҫ����

			//ept x = CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - 2*S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);//ֻҪ����lb�ͺ��ˣ�����Ҫ����

			printf("deleted %lu edges\n", x);
			printf("deleted %d nodes\n", deleted_nodes);
			m -= 2 * x;
			printf("After CTCP, n = %u, m = %lu\n\n", n - deleted_nodes, m / 2);




			for (ui i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];


			if (S <= K - 1)
			{

				//Choose vertex. Find 2-hop neighbours. Branch and Bound.
				ui* choosed = new ui[n];
				memset(choosed, 0, sizeof(ui) * n);

				for (int j = 0; j < n && m; j++) //&&ktruss_bar.size()<=UB
				{
					int min_key = n, min_id = -1;
					for (int k = 0; k < n; k++)
					{
						if (degree[k] < min_key && !choosed[k])
						{
							min_key = degree[k];
							min_id = k;
						}
					}
					choosed[min_id] = 1;

					if (min_key < lb - (S + 1))
					{
						if (degree[min_key] != 0)  
						{
							while (!Qv.empty()) Qv.pop();
							while (!Qe.empty()) Qe.pop();
							Qv.push(min_id);
							deleted_nodes = 0;
							m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
							//m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - 2*S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
						}
						continue;
					}

					if (m == 0) break;

					assert(degree[min_id] == min_key);
					ui* ids = new ui[n];
					ui ids_n = 0;
					ui* rid = new ui[n];
					vector<pair<int, int> > vp; vp.reserve(m / 2);
					ui* Q = new ui[n];
					ui* t_degree = new ui[n];
					memset(t_degree, 0, sizeof(ui) * n);
					ui* exist = new ui[n];
					memset(exist, 0, sizeof(ui) * n);

					extract_subgraph_and_prune(lb, S + 1, min_id, ids, ids_n, rid, vp, Q, t_degree, exist, pend, deleted, edgelist_pointer);


					if (ids_n != 0)
						printf("\nids_n = %d\n", ids_n);

					if (ids_n >= lb)
					{
				
						bool ret;
						BnBiterall* bnbiter = new BnBiterall(S, lb, start_from_zero);
						bnbiter->load_graph(rid[min_id], ids_n, rid, vp, ids);

						bnbiter->rearrange_graph();   ///��ԭ

						int size = n * 2;
						ret = bnbiter->bb_s_truss_bar(ids_n, size);
						bnbiter->~BnBiterall();

						if (ret)
						{
							printf("K = %d and Find a k_truss of size: %d\n", K, size);
							ans = true;
							if (size < minsize) minsize = size;
						}
					}

					while (!Qv.empty()) Qv.pop();
					while (!Qe.empty()) Qe.pop();
					Qv.push(min_id);
					deleted_nodes = 0;
					m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
					//m -= 2 * CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - 2*S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);

				}

			}
			else if (S > K - 1)
			{
				ui* ids = new ui[n];
				ui ids_n = 0, rid_n = 0;
				ui* rid = new ui[n];
				vector<pair<int, int> > vp; vp.reserve(m / 2);
				for (ui i = 0; i < n; i++)
				{
					if (degree[i] > 0)
					{
						ids[ids_n] = i;
						rid[ids[ids_n]] = rid_n;
						ids_n++;
						rid_n++;
					}
				}

				for (ui i = 0; i < n; i++)
				{
					if (degree[i] > 0)
					{
						for (ept j = pstart[i]; j < pend[i]; j++)
						{
							if (edges[j] > i && !deleted[edgelist_pointer[j]]) vp.push_back(mp(rid[i], rid[edges[j]]));
						}
					}
				}

				int size = 2 * n;
				bool ret;
				BnBG* bnbg = new BnBG(S, lb);
				bnbg->load_graph(ids_n, rid, vp, ids);
				bnbg->rearrange_graph();
				ret = bnbg->bb_s_truss_bar_iterall(ids_n, start, timeflag, size);
				bnbg->~BnBG();

				if (ret)
				{
					printf("K = %d and Find a k_truss of size: %d\n", K, size);
					ans = true;
					if (size < minsize) minsize = size;
				}
			}


			//Return to old graph

			n = temp_n;
			m = temp_m;


			for (ui i = 0; i < n; i++)
			{
				pstart[i] = temp_pstart[i];
				pend[i] = temp_pend[i];

				for (ept j = temp_pstart[i]; j < temp_pend[i]; j++)
				{
					ui e = temp_edgelist_pointer[j];
					edges[j] = temp_edges[j];
					edgelist_pointer[j] = temp_edgelist_pointer[j];
					//printf("j = %d and %d-%d %d\n", j, i, edges[j], e);
					//printf("%d %d\n", edge_list[e << 1], temp_edge_list[e << 1]);
					edge_list[e << 1] = temp_edge_list[e << 1];
					edge_list[(e << 1) + 1] = temp_edge_list[(e << 1) + 1];
				}
			}


			for (ui i = 0; i < (temp_m >> 1); i++)
			{
				tri_cnt[i] = temp_tri_cnt[i];
				//edge_list[i << 1] = temp_edge_list[i << 1];
				//edge_list[(i << 1) + 1] = temp_edge_list[(i << 1) + 1];
			}


		}

		if (!ans)
		{
			printf("There is no %d-truss in the graph.\n", K);
		}
		return minsize;

	}

	int find_SCKT_k_truss_query(clock_t start, bool& timeflag, string method, int choosed_vertex, int s)
	{
		read_graph();
		max_truss();

		if (deleted_oldid[choosed_vertex])
		{
			printf("Query vertex has been deleted by the maximum ktruss\n");
			return -1;
		}

		ofstream outfile;
		outfile.open("/home/zhangqifan/min_k_truss/dataset/ktruss.txt");
		for (int i = 0; i < (m >> 1); i++)
		{
			int u = edge_list[i << 1], v = edge_list[(i << 1) + 1];
			if (start_from_zero)
				outfile << old_old_id[u] << " " << old_old_id[v] << endl;
			else
				outfile << old_old_id[u] + 1 << " " << old_old_id[v] + 1 << endl;
		}
		outfile.close();


		ui temp_n = n;
		ept temp_m = m;

		ept* temp_pstart = new ept[n];
		ept* temp_pend = new ept[n];
		ui* temp_edges = new ui[m];
		ui* temp_edgelist_pointer = new ui[m];
		ui* temp_tri_cnt = new ui[m >> 1];
		ui* temp_edge_list = new ui[m];

		memset(temp_pstart, 0, sizeof(ept) * n);
		memset(temp_pstart, 0, sizeof(ept) * n);
		memset(temp_edges, 0, sizeof(ui) * m);
		memset(temp_edgelist_pointer, 0, sizeof(ui) * m);
		memset(temp_tri_cnt, 0, sizeof(ui) * (m >> 1));
		memset(temp_edge_list, 0, sizeof(ui) * m);


		//store graph

		for (ui i = 0; i < n; i++)
		{
			temp_pstart[i] = pstart[i];
			temp_pend[i] = pend[i];
			for (ept j = pstart[i]; j < pend[i]; j++)
			{
				ui e = edgelist_pointer[j];
				temp_edges[j] = edges[j];
				temp_edgelist_pointer[j] = edgelist_pointer[j];
				temp_edge_list[e << 1] = edge_list[e << 1];
				temp_edge_list[(e << 1) + 1] = edge_list[(e << 1) + 1];
				//printf("j = %d and %d-%d %d\n",j, i, edges[j], temp_edgelist_pointer[j]);
			}
		}

		for (ui i = 0; i < (m >> 1); i++) temp_tri_cnt[i] = tri_cnt[i];


		bool ans = false;

		for (ui i = s; i >= K && !ans; i--)
		{

			//Preprocessing
			if ((int)(n - K) < 0) break;
			printf("n = %d, m = %d\n", n, m / 2);
			ui S = i - K;

			ui* active_edgelist = new ui[m >> 1];
			ui active_edgelist_n = m >> 1;
			for (ui i = 0; i < (m >> 1); i++) active_edgelist[i] = i;

			bool* deleted = new bool[m >> 1];
			memset(deleted, false, sizeof(bool) * (m >> 1));
			bool* exists = new bool[n];
			memset(exists, false, sizeof(bool) * n);

			ui* degree = new ui[n];
			for (ui i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];

			queue<ui> Qv;
			queue<ept> Qe;


			ui lb = S + K;
			printf("S = %d and LB = %d\n", S, lb);
			ui deleted_nodes = 0;

			bool* inq = new bool[n];
			memset(inq, false, sizeof(bool) * n);

			ept x = CTCP(inq, deleted_nodes, Qv, lb - S - 1, Qe, lb - S - 2, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, degree, pstart, pend, edges, exists);
		
			printf("deleted %lu edges\n", x);
			printf("deleted %d nodes\n", deleted_nodes);
			m -= 2 * x;
			printf("After CTCP, n = %u, m = %lu\n\n", n - deleted_nodes, m / 2);

			for (ui i = 0; i < n; i++) degree[i] = pend[i] - pstart[i];


	

			if (method == "BNBGLOBAL")
			{
				ui* ids = new ui[n];
				ui ids_n = 0, rid_n = 0;
				ui* rid = new ui[n];
				vector<pair<int, int> > vp; vp.reserve(m / 2);
				for (ui i = 0; i < n; i++)
				{
					if (degree[i] > 0)
					{
						ids[ids_n] = i;
						rid[ids[ids_n]] = rid_n;
						ids_n++;
						rid_n++;
					}
				}

				for (ui i = 0; i < n; i++)
				{
					if (degree[i] > 0)
					{
						for (ept j = pstart[i]; j < pend[i]; j++)
						{
							if (edges[j] > i && !deleted[edgelist_pointer[j]]) vp.push_back(mp(rid[i], rid[edges[j]]));
						}
					}
				}

				bool ret;
				BnBG* bnbg = new BnBG(S, lb);
				bnbg->load_graph(ids_n, rid, vp, ids);
				bnbg->rearrange_graph();
				ret = bnbg->bb_s_truss_bar(ids_n, start, timeflag);
				bnbg->~BnBG();

				if (ret)
				{
					printf("K = %d and Find a min_k_truss of size: %d\n", K, S + K);
					ans = true;
					return S + K;
				}
				continue;
			}
			//Choose vertex. Find 2-hop neighbours. Branch and Bound.
			ui* choosed = new ui[n];
			memset(choosed, 0, sizeof(ui) * n);

			int query_v;
			for (int ii = 0; ii < n; ii++)
			{
				if (old_old_id[ii] == choosed_vertex)
				{
					query_v = ii;
					break;
				}
			}

			for (int j = 0; j < n && m; j++) //&&ktruss_bar.size()<=UB
			{

				int min_key = n, min_id = -1;
				min_id = query_v;
				min_key = degree[min_id];

				if (m == 0) break;

				assert(degree[min_id] == min_key);
				ui* ids = new ui[n];
				ui ids_n = 0;
				ui* rid = new ui[n];
				vector<pair<int, int> > vp; vp.reserve(m / 2);
				ui* Q = new ui[n];
				ui* t_degree = new ui[n];
				memset(t_degree, 0, sizeof(ui) * n);
				ui* exist = new ui[n];
				memset(exist, 0, sizeof(ui) * n);

				extract_subgraph_and_prune(lb, S + 1, min_id, ids, ids_n, rid, vp, Q, t_degree, exist, pend, deleted, edgelist_pointer);

				if (ids_n != 0)
					printf("\nids_n = %d\n", ids_n);


				if (ids_n >= lb)
				{
					old_id = new ui[ids_n];
					for (int i = 0; i < ids_n; i++)
					{
						old_id[rid[ids[i]]] = old_old_id[ids[i]];
					}
					
					BnB1* bnb = new BnB1(S, lb);
					bnb->load_graph(rid[min_id], ids_n, rid, vp, old_id);

					bnb->rearrange_graph();   ///��ԭ

					//ɾ
					/*bnb->test();
					exit(1);*/


					bool ret = bnb->bb_s_truss_bar(ids_n, start, timeflag);
					bnb->~BnB1();


					//test!!!!!!!!!!!!!!�ǵ�ɾ
					//ret = true;


					if (ret)
					{
						printf("choose vertex = %d\n", choosed_vertex);
						printf("K = %d and Find a min_k_truss of size: %d\n", K, S + K);
						ans = true;
						return S + K;
						break;
					}

					if (!ret)
					{
						if (S == K - 1) return -1;
						//printf("Cannot find query vertex based ans!\n");
						break;
						//return -1;
					}

				}

			}



			//Return to old graph

			n = temp_n;
			m = temp_m;


			for (ui i = 0; i < n; i++)
			{
				pstart[i] = temp_pstart[i];
				pend[i] = temp_pend[i];

				for (ept j = temp_pstart[i]; j < temp_pend[i]; j++)
				{
					ui e = temp_edgelist_pointer[j];
					edges[j] = temp_edges[j];
					edgelist_pointer[j] = temp_edgelist_pointer[j];
					//printf("j = %d and %d-%d %d\n", j, i, edges[j], e);
					//printf("%d %d\n", edge_list[e << 1], temp_edge_list[e << 1]);
					edge_list[e << 1] = temp_edge_list[e << 1];
					edge_list[(e << 1) + 1] = temp_edge_list[(e << 1) + 1];
				}
			}


			for (ui i = 0; i < (temp_m >> 1); i++)
			{
				tri_cnt[i] = temp_tri_cnt[i];
				//edge_list[i << 1] = temp_edge_list[i << 1];
				//edge_list[(i << 1) + 1] = temp_edge_list[(i << 1) + 1];
			}


		}

		if (!ans)
		{
			printf("There is no %d-truss in the graph.\n", K);
			return -1;
		}
	}


	void extract_subgraph_and_prune(ui lb, ui K, ui u, ui* ids, ui& ids_n, ui* rid, vector<pair<int, int> >& vp, ui* Q, ui* degree, ui* exists, ept* pend, bool* deleted, ui* edgelist_pointer) {
		vp.clear();
		ids_n = 0; ids[ids_n++] = u; exists[u] = 1;
		ui u_n = pstart[u];
		for (ept i = pstart[u]; i < pend[u]; i++) if (!deleted[edgelist_pointer[i]]) {
			edges[u_n] = edges[i]; edgelist_pointer[u_n++] = edgelist_pointer[i];
			ui v = edges[i];
			ids[ids_n++] = v; exists[v] = 2;
		}
		pend[u] = u_n;

		ui Q_n = 0;
		for (ui i = 1; i < ids_n; i++) {
			u = ids[i];
			u_n = pstart[u];
			degree[u] = 0;
			for (ept j = pstart[u]; j < pend[u]; j++) if (!deleted[edgelist_pointer[j]]) {
				edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (exists[edges[j]] == 2) ++degree[u];
			}
			pend[u] = u_n;
			//printf("degree + 2K = %d, lb = %d\n", degree[u] + 2 * K, lb);
			if (degree[u] + 2 * K < lb) Q[Q_n++] = u;
		}
		for (ui i = 0; i < Q_n; i++) exists[Q[i]] = 10;
		//printf("Qn = %d, ids_n = %d\n", Q_n, ids_n);

		for (ui i = 0; i < Q_n; i++) {
			u = Q[i];
			exists[u] = 10;
			for (ept j = pstart[u]; j < pend[u]; j++) if (exists[edges[j]] == 2) {
				if ((degree[edges[j]]--) + 2 * K == lb)
				{
					assert(Q_n < m / 2);
					Q[Q_n++] = edges[j];
				}

			}
		}
		assert(Q_n <= ids_n);
		if (ids_n - Q_n + K - 1 < lb) { //����u,u���ھӣ������ٽ���k-1��u�ķ��ھ�
			for (ui i = 0; i < ids_n; i++) exists[ids[i]] = 0;
			ids_n = 0;
			return;
		}

		ui nr_size = ids_n;//nr_size - 1��u��һ���ھӣ�prune֮ǰ��������

		for (ui i = 1; i < nr_size; i++) if (exists[ids[i]] == 2) {//��u��2���ھ�ֻ�ӻ�ûɾ��1���ھӳ���
			u = ids[i];
			for (ept j = pstart[u]; j < pend[u]; j++) {
				if (!exists[edges[j]]) {
					ids[ids_n++] = edges[j];
					exists[edges[j]] = 3;
					degree[edges[j]] = 1;
				}
				else if (exists[edges[j]] == 3) ++degree[edges[j]];//2���ھӵ�degreeָ����2���ھ���1���ھ��еĶ���
			}
		}

#ifndef NDEBUG
		//printf("Entire list: ");
		//for(ui i = 0;i < nr_size;i ++) printf(" %u", ids[i]);
		//printf("\n");
#endif

		//prune 1-hop neighbour
		ui new_size = 1;
		for (ui i = 1; i < nr_size; i++) {
			if (exists[ids[i]] == 10) exists[ids[i]] = 0;
			else ids[new_size++] = ids[i];
		}
#ifndef NDEBUG
		if (new_size + Q_n != nr_size) {
			printf("new_size: %u, Q_n: %u, nr_size: %u\n", new_size, Q_n, nr_size);
			printf("New list: ");
			for (ui i = 0; i < new_size; i++) printf(" %u", ids[i]);
			printf("\n");
			printf("Pruned list: ");
			for (ui i = 0; i < Q_n; i++) printf(" %u", Q[i]);
			printf("\n");
		}
#endif

		//prune 2-hop neighbour
		assert(new_size + Q_n == nr_size);
		ui old_nr_size = nr_size;
		nr_size = new_size;
		for (ui i = old_nr_size; i < ids_n; i++) {
			if (degree[ids[i]] + 2 * K < lb + 2) exists[ids[i]] = 0;
			else ids[new_size++] = ids[i];

		}
		ids_n = new_size;
#ifndef NDEBUG
		assert(exists[ids[0]] == 1);
		for (ui i = 1; i < nr_size; i++) assert(exists[ids[i]] == 2);
		for (ui i = nr_size; i < ids_n; i++) assert(exists[ids[i]] == 3);
#endif

		
		for (ui i = 0; i < ids_n; i++) {
			assert(exists[ids[i]]);
			rid[ids[i]] = i;
		}


		for (ui i = 0; i < nr_size; i++) {
			u = ids[i];
			for (ept j = pstart[u]; j < pend[u]; j++) if (exists[edges[j]] && edges[j] > u) {
				assert(!deleted[edgelist_pointer[j]]);
				vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
		}

		for (ui i = nr_size; i < ids_n; i++) {
			u = ids[i];
			u_n = pstart[u];
			for (ept j = pstart[u]; j < pend[u]; j++) if (!deleted[edgelist_pointer[j]]) {
				edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (edges[j] > u && exists[edges[j]]) vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
			pend[u] = u_n;
		}
		for (ui i = 0; i < ids_n; i++) exists[ids[i]] = 0;
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


#ifndef  NDEBUG
		long long sum = 0;
		for (ui i = 0; i < n; i++) sum += pend[i] - pstart[i];
		//cout << sum * 2 << " " << m << endl;
		assert(sum * 2 == m); 
#endif // ! NDEBUG

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
		printf("Total number of triangles:%d\n", cnt);
#endif // ! NDEBUG
	}


	void reorganize_oriented_graph(ui n, ui* tri_cnt, ui* edge_list, ept* pstart, ept* pend, ept* pend2, ui* edges, ui* edgelist_pointer)
	{
		memset(tri_cnt, 0, (m / 2) * sizeof(ui));
		for (ui i = 0; i < n; i++) pend2[i] = pend[i];


		//
		ept pos = 0;

		for (ui i = 0; i < n; i++)
		{
			for (ept j = pstart[i]; j < pend[i]; j++)
			{
				tri_cnt[pos >> 1] = edgelist_pointer[j];
				edge_list[pos++] = i, edge_list[pos++] = edges[j];

				ept& k = pend2[edges[j]];
				edgelist_pointer[j] = edgelist_pointer[k] = (pos >> 1) - 1; //edgelist_pointer��¼�ߵı��,j��k�Ǳߵĵ�ַ
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
};



#endif