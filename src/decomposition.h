#ifndef DECOMPOSITION_H_
#define DECOMPOSITION_H_

#include <vector>

// SharpSAT uses this anyway globally so...
using namespace std;

namespace decomp {
pair<int, vector<int>> ComputeTreewidth(const vector<vector<int>>& graph, double time);
} // namespace decomp

#endif /* DECOMPOSITION_H_ */
