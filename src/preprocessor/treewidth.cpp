#include "treewidth.hpp"

#include <chrono>
#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <cassert>
#include <queue>
#include <algorithm>
#include <sstream>

#include "utils.hpp"

namespace sspp {

namespace decomp {
namespace {
string TmpInstance(int a, int b, int c, string tmp_dir) {
  std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
  auto duration = now.time_since_epoch();
  uint64_t micros = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
  return tmp_dir+"/instance"+std::to_string(micros)+"_"+std::to_string(a)+"_"+std::to_string(b)+"_"+std::to_string(c)+".tmp";
}
} // namespace

TreeDecomposition Treedecomp(const Graph& graph, double time, string tmp_dir) {
	int n = graph.n();
	if (n == 0) {
		TreeDecomposition dec(0, 0);
		return dec;
	}
	if (n == 1) {
		TreeDecomposition dec(1, 1);
		dec.SetBag(1, {0});
		return dec;
	}
	if (time < 0.099) {
		TreeDecomposition dec(1, n);
		vector<int> all;
		for (int i = 0; i < n; i++) {
			all.push_back(i);
		}
		dec.SetBag(1, all);
		return dec;
	}
	assert(n >= 2);
	auto es = graph.Edges();
	int m = es.size();
	string tmp1 = TmpInstance(n, m, 1, tmp_dir);
	string tmp2 = TmpInstance(n, m, 2, tmp_dir);
	std::ofstream out(tmp1);
	out<<"p tw "<<n<<" "<<m<<'\n';
	for (auto e : es) {
		out << e.F+1 << " " << e.S+1 << '\n';
	}
	out << std::flush;
	out.close();
	cout<<"c o Primal edges "<<es.size()<<endl;
	string tw_binary = "./flow_cutter_pace17";
	string cmd = "timeout " + to_string(time) + "s " + tw_binary + " <" + tmp1 + " >" + tmp2 + " 2>/dev/null";
	cout << "c o CMD: " << cmd << endl;
	int status = system(cmd.c_str());
	assert(status >= 0);
	assert(WIFEXITED(status));
	assert(WEXITSTATUS(status) == 124); // TIMEOUT timed out
	cout << "c o tw finish ok" << endl;
	TreeDecomposition dec(0, 0);
	std::ifstream in(tmp2);
	string tmp;
	int claim_width = 0;
	while (getline(in, tmp)) {
		std::stringstream ss(tmp);
		ss>>tmp;
		if (tmp == "c") continue;
		if (tmp == "s") {
			ss>>tmp;
			assert(tmp == "td");
			int bs,nn;
			ss>>bs>>claim_width>>nn;
			assert(nn == n);
			claim_width--;
			dec = TreeDecomposition(bs, nn);
		} else if (tmp == "b") {
			int bid;
			ss>>bid;
			vector<int> bag;
			int v;
			while (ss>>v) {
				bag.push_back(v-1);
			}
			dec.SetBag(bid, bag);
		} else {
			int a = stoi(tmp);
			int b;
			ss>>b;
			dec.AddEdge(a, b);
		}
	}
	in.close();
	system(("rm -f " + tmp1).c_str());
	system(("rm -f " + tmp2).c_str());
	assert(dec.Width() <= claim_width);
	cout << "c o width " << dec.Width() << endl;
	assert(dec.Verify(graph));
	return dec;
}
} // namespace decomp
} // namespace sspp