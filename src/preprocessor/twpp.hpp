#pragma once

#include "graph.hpp"

#include "utils.hpp"

#include <queue>

namespace sspp {
class TWPP {
 public:
  Graph PP(Graph graph);
 private:
  void UpdLB(const Graph& graph);
  bool GreedyElim(Graph graph);
  bool Reduce(Graph graph);
  bool AlmostClique(Graph graph);
  bool bad_ = false;
  std::vector<Edge> fill_;
  int lb_ = 0;
  std::queue<Graph> kernel_;
};
} // namespace sspp