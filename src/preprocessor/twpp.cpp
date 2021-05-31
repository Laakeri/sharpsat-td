#include "twpp.hpp"

#include <vector>
#include <cassert>
#include <map>
#include <iostream>
#include <set>

namespace sspp {
namespace {
class UniqQue {
 private:
  std::queue<int> q_;
  std::vector<int> inq_;
  int n_;
 public:
  explicit UniqQue(int n);
  void Add(int x);
  void Add(const std::vector<int>& xs);
  int Pop();
  bool Empty() const;
};
UniqQue::UniqQue(int n) : inq_(n), n_(n) {}
void UniqQue::Add(int x) {
  assert(x>=0 && x<n_);
  if (!inq_[x]) {
    q_.push(x);
    inq_[x] = true;
  }
}
void UniqQue::Add(const std::vector<int>& xs) {
  for (int x : xs) {
    Add(x);
  }
}
int UniqQue::Pop() {
  assert(q_.size()>0);
  int r = q_.front();
  assert(inq_[r]);
  inq_[r] = false;
  q_.pop();
  return r;
}
bool UniqQue::Empty() const {
  return q_.empty();
}
} // namespace

namespace mcs {
std::vector<int> Mcs(const Graph& graph) {
  std::vector<int> order(graph.n());
  std::vector<int> label(graph.n());
  std::vector<char> rm(graph.n());
	std::vector<std::vector<int> > labels(graph.n());
  for (int i = 0; i < graph.n(); i++) labels[i].clear();
  for (int i = 0; i < graph.n(); i++) labels[0].push_back(i);
  int max_label = 0;
  for (int it = graph.n() - 1; it >= 0; it--) {
    if (labels[max_label].size() == 0) {
      max_label--;
      it++;
      continue;
    }
    int x = labels[max_label].back();
    labels[max_label].pop_back();
    if (rm[x]) {
      it++;
      continue;
    }
    order[it] = x;
    for (int nx : graph.Neighbors(x)) {
      if (!rm[nx]) {
        label[nx]++;
        labels[label[nx]].push_back(nx);
        max_label = std::max(max_label, label[nx]);
      }
    }
    rm[x] = true;
  }
  return order;
}

struct McsMOutput {
  std::vector<Edge> fill_edges;
  std::vector<int> elimination_order;
  std::vector<char> is_maximal_clique_point;
};

McsMOutput McsM(const Graph& graph) {
  Timer mcst;
  mcst.start();
  std::vector<int> label(graph.n());
  std::vector<std::vector<int> > reach(graph.n());
  Bitset rm(graph.n());
  Bitset rc(graph.n());
  for (int i = 0; i < graph.n(); i++) reach[i].clear();
  std::vector<Edge> fill;
  std::vector<int> order(graph.n());
  std::vector<char> is_maximal_point(graph.n());
  // TODO: maybe better variable names?
  int chunks = rm.chunks_;
  int prev_label = -1;
  for (int it = graph.n() - 1; it >= 0; it--) {
    int x = 0;
    int max_label = 0;
    for (int i = 0; i < graph.n(); i++) {
      if (!rm.Get(i) && label[i] >= max_label) {
        x = i;
        max_label = label[x];
      }
    }
    assert(!rm.Get(x) && label[x] < graph.n());
    order[it] = x;
    is_maximal_point[it] = (label[x] <= prev_label);
    prev_label = label[x];
    rc.Clear();
    rc.SetTrue(x);
    rm.SetTrue(x);
    for (int y : graph.Neighbors(x)) {
      if (!rm.Get(y)) {
        rc.SetTrue(y);
        reach[label[y]].push_back(y);
      }
    }
    for (int i = 0; i < graph.n(); i++) {
      while (!reach[i].empty()) {
        int y = reach[i].back();
        reach[i].pop_back();
        for (int j=0;j<chunks;j++){
          uint64_t td = graph.adj_mat2_[y].data_[j] & (~rm.data_[j]) & (~rc.data_[j]);
          while (td) {
            int z = j*BITS + __builtin_ctzll(td);
            td &= ~-td;
            rc.SetTrue(z);
            if (label[z] > i) {
              reach[label[z]].push_back(z);
              label[z]++;
              fill.push_back({x, z});
            } else {
              reach[i].push_back(z);
            }
          }
        }
      }
    }
    for (int y : graph.Neighbors(x)) {
      if (!rm.Get(y)) label[y]++;
    }
  }
  return McsMOutput({fill, order, is_maximal_point});
}

int Treewidth(const Graph& graph) {
  if (graph.m() == 0) return 0;
  std::vector<int> order = Mcs(graph);
  std::vector<int> inv_order = PermInverse(order);
  int treewidth = 0;
  for (int i = 0; i < graph.n(); i++) {
    int x = order[i];
    int nb = 0;
    std::vector<int> clq;
    for (int nx : graph.Neighbors(x)) {
      if (inv_order[nx] > i) {
        nb++;
        clq.push_back(nx);
      }
    }
    treewidth = std::max(treewidth, nb);
  }
  return treewidth;
}

std::vector<Graph> Atoms(const Graph& graph, const McsMOutput& mcs_m_output) {
  if (graph.n() == 0 || graph.m() == 0) {
    return {};
  }
  Graph filled_graph(graph);
  filled_graph.AddEdges(mcs_m_output.fill_edges);
  std::vector<char> rm(graph.n());
  std::vector<char> block(graph.n());
  std::vector<Graph> atoms;
  for (int it = 0; it < graph.n(); it++) {
    int x = mcs_m_output.elimination_order[it];
    if (mcs_m_output.is_maximal_clique_point[it]) {
      std::vector<int> cand_clique;
      for (int nx : filled_graph.Neighbors(x)) {
        if (!rm[nx]) cand_clique.push_back(nx);
      }
      if (graph.IsClique(cand_clique)) {
        for (int y : cand_clique) block[y] = 1;
        std::vector<int> component = graph.FindComponentAndMark(x, block);
        for (int y : cand_clique) block[y] = 0;
        component.insert(component.end(), cand_clique.begin(), cand_clique.end());
        atoms.push_back(Graph(graph.EdgesIn(component)));
      }
    }
    rm[x] = true;
  }
  int found = 0; // TODO
  for (int x = 0; x < graph.n(); x++) {
    if (!block[x]) {
      std::vector<int> component = graph.FindComponentAndMark(x, block);
      atoms.push_back(Graph(graph.EdgesIn(component)));
      found++;
    }
  }
  assert(found == 1); // TODO
  return atoms;
}

int Heur(const Graph& graph, int v) {
  std::set<Edge> fes;
  for (const Bitset& cn : graph.CompNeighsBit(graph.adj_mat2_[v])) {
    for (auto fe : graph.FillEdges(cn.Elements())) {
      fes.insert(fe);
    }
  }
  return fes.size();
}

void LbTriang(Graph& graph) {
  Timer lbt;
  lbt.start();
  std::priority_queue<std::pair<int, int>> q;
  std::vector<int> hs(graph.n());
  for (int i=0;i<graph.n();i++) {
    hs[i] = Heur(graph, i);
    q.push({-hs[i], i});
  }
  while (!q.empty()) {
    int fi = -q.top().first;
    int v = q.top().second;
    q.pop();
    if (hs[v] != fi) continue;
    hs[v] = -1;
    std::set<int> upd1;
    for (const Bitset& cn : graph.CompNeighsBit(graph.adj_mat2_[v])) {
      for (auto fe : graph.FillEdges(cn.Elements())) {
        graph.AddEdge(fe);
        upd1.insert(fe.first);
        upd1.insert(fe.second);
      }
    }
    std::set<int> upd2;
    for (int u : upd1) {
      for (int nb : graph.Neighbors(u)) {
        if (hs[nb] != -1) {
          upd2.insert(nb);
        }
      }
    }
    for (int u : upd2) {
      assert(hs[u] != -1);
      hs[u] = Heur(graph, u);
      q.push({-hs[u], u});
    }
  }
}
} // namespace

void TWPP::UpdLB(const Graph& graph) {
	lb_ = std::max(lb_, graph.Degeneracy());
}

// THIS SHOULD BE CALLED WITH GRAPH THAT HAS BEEN BROKE INTO ATOMS
bool TWPP::GreedyElim(Graph graph) {
	UniqQue q(graph.n());
	q.Add(graph.Vertices());
	bool reduced = false;
	while (!q.Empty()) {
		int v = q.Pop();
		if (graph.IsClique(graph.Neighbors(v))) {
			reduced = true;
			auto nbs = graph.Neighbors(v);
			graph.RemoveEdgesBetween(v, nbs);
			q.Add(nbs);
			lb_ = std::max(lb_, (int)nbs.size());
		} else {
			for (auto cn : graph.CompNeighsBit(graph.adj_mat2_[v])) {
				assert(graph.IsMinsep(cn));
				if (graph.IsAlmostClique(cn.Elements())) {
					auto fes = graph.FillEdges(cn);
					graph.AddEdges(fes);
					Append(fill_, graph.MapBack(fes));
					lb_ = std::max(lb_, cn.Popcount());
					reduced = true;
					for (auto e : fes) {
						q.Add(graph.Neighbors(e.first));
						q.Add(graph.Neighbors(e.second));
					}
					break;
				}
			}
		}
	}
	if (reduced) {
		Graph t_graph(graph.Edges());
		if (t_graph.m() > 0) {
			t_graph.InheritMap(graph);
			kernel_.push(t_graph);
		}
		return true;
	} else {
		return false;
	}
}

bool TWPP::AlmostClique(Graph graph) {
	assert(graph.IsConnected());
	std::vector<int> found;
	int maxc = graph.n() + 1;
	for (int i=0;i<graph.n();i++) {
		Graph t_graph = graph;
		auto nbs = t_graph.Neighbors(i);
		assert((int)nbs.size() > 2);
		t_graph.RemoveEdgesBetween(i, nbs);
		mcs::McsMOutput minimal_triangulation = mcs::McsM(t_graph);
		t_graph.AddEdges(minimal_triangulation.fill_edges);
  	Bitset not_rm(graph.n());
  	not_rm.FillTrue();
		for (int it = 0; it < graph.n(); it++) {
	    int x = minimal_triangulation.elimination_order[it];
	    if (minimal_triangulation.is_maximal_clique_point[it]) {
	      Bitset cand_clique_b = t_graph.adj_mat2_[x] & not_rm;
	      cand_clique_b.SetFalse(x);
	      if (cand_clique_b.Popcount()>0 && graph.IsClique(cand_clique_b)) {
	      	auto cand_clique = cand_clique_b.Elements();
	      	assert(graph.IsClique(cand_clique));
	      	assert(t_graph.IsMinsep(cand_clique));
	      	assert(!graph.IsMinsep(cand_clique));
	      	cand_clique.push_back(i);
		  		SortAndDedup(cand_clique);
		  		assert(graph.IsMinsep(cand_clique));
		  		assert(!graph.IsClique(cand_clique));
		  		assert(graph.IsAlmostClique(cand_clique));
		  		auto comps = graph.NComponents(cand_clique);
		  		lb_ = std::max(lb_, (int)cand_clique.size());
		  		int tmaxc = 0;
		  		for (const auto& comp : graph.NComponents(cand_clique)) {
		  			tmaxc = std::max(tmaxc, (int)comp.size());
		  		}
		  		if (tmaxc < maxc) {
		  			maxc = tmaxc;
		  			found = cand_clique;
		  		}
	      }
	    }
	    not_rm.SetFalse(x);
	  }
	}
	if (!found.empty()) {
		assert(maxc < graph.n());
		assert(graph.IsMinsep(found));
		assert(!graph.IsClique(found));
		assert(graph.IsAlmostClique(found));
		auto fes = graph.FillEdges(found);
		assert(fes.size()>0);
		graph.AddEdges(fes);
		Append(fill_, graph.MapBack(fes));
		kernel_.push(graph);
		return true;
	} else {
		assert(maxc == graph.n() + 1);
		return false;
	}
}

bool TWPP::Reduce(Graph graph) {
	if (graph.m() == 0) return true;
  UpdLB(graph);
  mcs::McsMOutput minimal_triangulation = mcs::McsM(graph);
  Graph fill_graph = graph;
  fill_graph.AddEdges(minimal_triangulation.fill_edges);
  if (minimal_triangulation.fill_edges.size() <= 1 || mcs::Treewidth(fill_graph) <= lb_) {
  	Append(fill_, graph.MapBack(minimal_triangulation.fill_edges));
  	lb_ = std::max(lb_, mcs::Treewidth(fill_graph));
  	return true;
  }
  std::vector<Graph> atoms = mcs::Atoms(graph, minimal_triangulation);
  if (atoms.size() == 0) return true;
  if (atoms.size() > 1) {
  	for (auto& atom : atoms) {
  		atom.InheritMap(graph);
  		if (atom.n() >= 4) {
  			kernel_.push(atom);
  		}
  	}
  	return true;
  } else {
  	atoms[0].InheritMap(graph);
  	graph = atoms[0];
  }
  if (GreedyElim(graph)) {
  	return true;
  }
  if (graph.m() <= 10000) {
    Graph lb_graph = graph;
    mcs::LbTriang(lb_graph);
    int ub = mcs::Treewidth(lb_graph);
    if (ub <= lb_) {
      Append(fill_, graph.MapBack(graph.FillEdges(lb_graph)));
      return true;
    }
  }
  if (graph.m() <= 500) {
  	 if (AlmostClique(graph)) {
	  	return true;
	  }
  }
  kernel_.push(graph);
  return false;
}

Graph TWPP::PP(Graph graph) {
  assert(!bad_);
  bad_ = true;
  assert(kernel_.size() == 0);
  assert(lb_ == 0);
  kernel_.push(graph);
  int failed = 0;
  while (failed < (int)kernel_.size()) {
  	if (Reduce(kernel_.front())) {
  		failed = 0;
  	} else {
  		failed++;
  	}
  	kernel_.pop();
  }
  graph.AddEdges(fill_);
  return graph;
}

} // namespace sspp