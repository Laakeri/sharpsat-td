#include "decomposition.h"

#include <chrono>
#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <cassert>
#include <queue>
#include <algorithm>
#include <sstream>

using Graph = vector<vector<int>>;

namespace decomp {
namespace {
void bfs(int x, const Graph& graph, vector<int>& u, vector<int>& vs) {
	int n = graph.size()-1;
	assert(!u[x]);
	u[x] = x;
	vs.push_back(x);
	for (int i = 0; i < (int)vs.size(); i++) {
		int t = vs[i];
		assert(u[t]);
		for (int nt : graph[t]) {
			assert(1 <= nt && nt <= n);
			if (!u[nt]) {
				u[nt] = x;
				vs.push_back(nt);
			} else {
				assert(u[nt] == x);
			}
		}
	}
}

vector<pair<Graph, vector<int>>> Components(const Graph& graph) {
	int n = graph.size()-1;
	assert(n >= 1);
	vector<int> u(n+1);
	vector<pair<Graph, vector<int>>> ret;
	int tot = 0;
	for (int i = 1; i <= n; i++) {
		if (u[i]) continue;
		vector<int> vs;
		bfs(i, graph, u, vs);
		sort(vs.begin(), vs.end());
		int nn = vs.size();
		tot += nn;
		Graph ng(nn+1);
		vector<int> n_map(nn+1);
		for (int j = 0; j < (int)vs.size(); j++) {
			int t = vs[j];
			assert(u[t] == i);
			n_map[j+1] = t;
			for (int nt : graph[t]) {
				assert(u[nt] == i);
				int nj = lower_bound(vs.begin(), vs.end(), nt) - vs.begin();
				assert(vs[nj] == nt);
				assert(j != nj);
				ng[j+1].push_back(nj+1);
			}
		}
		ret.push_back({ng, n_map});
	}
	assert(tot == n);
	return ret;
}

string TmpInstance(int a, int b, int c) {
  chrono::time_point<chrono::system_clock> now = std::chrono::system_clock::now();
  auto duration = now.time_since_epoch();
  uint64_t micros = chrono::duration_cast<chrono::microseconds>(duration).count();
  string dir = "td_tmp";
  system(("mkdir -p " + dir).c_str());
  return dir+"/instance"+std::to_string(micros)+"_"+std::to_string(a)+"_"+std::to_string(b)+"_"+std::to_string(c)+".tmp";
}

int TWCen(int x, int p, int n, const vector<vector<int>>& bags, const vector<vector<int>>& tree, int& cen) {
	assert(x >= 1 && x < (int)bags.size());
	assert(bags.size() == tree.size());
	assert(cen == 0);
	int intro = 0;
	for (int nx : tree[x]) {
		if (nx == p) continue;
		int cintro = TWCen(nx, x, n, bags, tree, cen);
		intro += cintro;
		if (cintro >= n/2) {
			assert(cen);
			return intro;
		}
	}
	for (int v : bags[x]) {
		if (!binary_search(bags[p].begin(), bags[p].end(), v)) {
			intro++;
		}
	}
	if (intro >= n/2) {
		cen = x;
	}
	return intro;
}

void TWDes(int x, int p, int d, const vector<vector<int>>& bags, const vector<vector<int>>& tree, vector<int>& ret) {
	assert(x >= 1 && x < (int)bags.size());
	assert(bags.size() == tree.size());
	assert(d >= 1);
	bool new_vs = false;
	for (int v : bags[x]) {
		if (ret[v] == 0) {
			new_vs = true;
		} else {
			assert(ret[v] <= d);
			assert(binary_search(bags[p].begin(), bags[p].end(), v));
		}
	}
	if (new_vs) {
		d++;
		for (int v : bags[x]) {
			if (ret[v] == 0) {
				ret[v] = d;
			}
		}
	}
	for (int nx : tree[x]) {
		if (nx != p) {
			TWDes(nx, x, d, bags, tree, ret);
		}
	}
}

pair<int, vector<int>> Treewidth(const Graph& graph, double time) {
	int n = graph.size()-1;
	assert(n >= 1);
	if (n == 1) {
		return {1, {0, 1}};
	}
	if (n == 2) {
		return {2, {0, 1, 1}};
	}
	vector<int> ret(n + 1);
	if (time < 0.1) {
		cout << "WARNING: Low time budjet" << endl;
		for (int i = 1; i <= n; i++) {
			ret[i] = 1;
		}
		return {n, ret};
	}
	int m = 0;
	for (int i = 1; i <= n; i++) {
		m += (int)graph[i].size();
	}
	assert(m%2 == 0);
	m /= 2;
	string tmp1 = TmpInstance(n, m, 1);
	string tmp2 = TmpInstance(n, m, 2);
	ofstream out(tmp1);
	out<<"p tw "<<n<<" "<<m<<endl;
	for (int i = 1; i <= n; i++) {
		for (int ni : graph[i]) {
			assert(ni != i && ni >= 1 && ni <= n);
			if (i < ni) {
				out << i << " " << ni << endl;
			}
		}
	}
	out.close();
	string tw_binary = "../../../flow-cutter-pace17/flow_cutter_pace17";
	string cmd = "timeout " + to_string(time) + "s " + tw_binary + " <" + tmp1 + " >" + tmp2 + " 2>/dev/null";
	cout << "CMD:" << endl;
	cout << cmd << endl;
	int status = system(cmd.c_str());
	assert(status >= 0);
	assert(WIFEXITED(status));
	assert(WEXITSTATUS(status) == 124); // TIMEOUT timed out
	cout << "tw finish ok" << endl;
	int width = 0;
	vector<vector<int>> bags;
	vector<vector<int>> tree;
	ifstream in(tmp2);
	string tmp;
	int real_width = 0;
	while (getline(in, tmp)) {
		stringstream ss(tmp);
		ss>>tmp;
		if (tmp == "c") continue;
		if (tmp == "s") {
			ss>>tmp;
			assert(tmp == "td");
			assert(bags.empty());
			int bs,nn;
			ss>>bs>>width>>nn;
			assert(nn == n);
			bags.resize(bs+1);
			tree.resize(bs+1);
			cout << "claim width " << width << endl;
		} else if (tmp == "b") {
			int bid;
			ss>>bid;
			assert(bid >= 1 && bid < (int)bags.size());
			assert(bags[bid].empty());
			int v;
			while (ss>>v) {
				assert(v >= 1 && v <= n);
				bags[bid].push_back(v);
			}
			assert((int)bags[bid].size() <= width);
			sort(bags[bid].begin(), bags[bid].end());
			bags[bid].erase(unique(bags[bid].begin(), bags[bid].end()), bags[bid].end());
			real_width = max(real_width, (int)bags[bid].size());
		} else {
			int a = stoi(tmp);
			assert(1 <= a && a < (int)bags.size());
			int b;
			ss>>b;
			assert(a != b);
			tree[a].push_back(b);
			tree[b].push_back(a);
		}
	}
	in.close();
	system(("rm -f " + tmp1).c_str());
	system(("rm -f " + tmp2).c_str());
	assert(bags[0].empty());
	assert(tree[0].empty());
	int centroid = 0;
	TWCen(1, 0, n, bags, tree, centroid);
	assert(centroid >= 1 && centroid < (int)bags.size());
	TWDes(centroid, 0, 1, bags, tree, ret);
	for (int i = 1; i <= n; i++) {
		assert(ret[i] >= 1);
	}
	assert(real_width <= width);
	return {real_width, ret};
}
} // namespace

pair<int, vector<int>> ComputeTreewidth(const Graph& graph, double time) {
	int n = graph.size()-1;
	if (n == 0) {
		return {1, {}};
	}
	assert(n >= 1);
	assert(time > 0.45);
	double time_per_var = time / (double)n;
	vector<int> ret(n + 1);
	int w = 0;
	for (const pair<Graph, vector<int>>& comp : Components(graph)) {
		int tn = comp.first.size()-1;
		double ttime = (double)tn * time_per_var;
		auto ww = Treewidth(comp.first, ttime);
		assert((int)ww.second.size() == tn+1);
		assert(ww.second[0] == 0);
		assert(ww.first >= 1);
		assert(ww.first <= tn);
		w = max(w, ww.first);
		for (int i = 1; i <= tn; i++) {
			assert(ww.second[i] >= 1 && ww.second[i] <= tn);
			ret[comp.second[i]] = ww.second[i];
		}
	}
	assert(w >= 1);
	for (int i = 1; i <= n; i++) {
		assert(ret[i] >= 1 && ret[i] <= n);
	}
	cout << "TREEWIDTH " << w << endl;
	return {w, ret};
}

} // namespace decomp