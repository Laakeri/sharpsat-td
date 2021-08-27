#pragma once

#include "utils.hpp"
#include "instance.hpp"

#include <map>

namespace sspp {
class Preprocessor {
 public:
 	Instance Preprocess(Instance ins, const string& techniques);
 	Instance Preprocess(int vars_, vector<vector<Lit>> clauses_, string techniques);
 	int FreeVars() const;
 	void SetMaxGTime(double time);
 	void SetMaxSparsTime(double time);
 private:
 	bool DoTechniques(const string& techniques, int l, int r);
 	bool EliminateDefSimplicial();
 	void MergeAdjEquivs();
  void Sparsify();
 	Instance MapBack();
 	Instance UnsatInst();
 	void BackBone();
 	void PropStren();
 	void Subsume();
 	void FailedLiterals();
 	void TakeUnits(vector<vector<Lit>>& cls);
 	void MapClauses(vector<vector<Lit>>& cls, Var& nvars, vector<Var>& map_to);
 	void Tighten(bool loop);
 	int vars = 0;
 	vector<vector<Lit>> clauses, learned_clauses;
 	Timer timer;

 	int orig_vars = 0;
 	vector<Var> var_map;
 	// 1 positive 2 negative
 	vector<char> assign;
 	bool unsat = false;

 	bool weighted = false;
 	vector<double> weights;
 	int free_vars = 0;

 	double max_g_time = 1e9;
 	double max_s_time = 1e9;

 	Timer g_timer;
 	Timer s_timer;

 	void eqdfs(Lit lit, Lit e, const vector<vector<Lit>>& eq, vector<Lit>& eqc);
};
inline bool ValidTechniques(const string& techniques, bool weighted) {
	int d = 0;
	for (char c : techniques) {
		if (c == '[') {
			d++;
		} else if (c == ']') {
			d--;
			if (d < 0) return false;
		} else if (c == 'F' // failed literals by propagation
			      || c == 'P' // vivification by propagation
			      || c == 'V' // complete vivification
			      || c == 'S' // sparsification
			      || c == 'E') { // equivalence merging
			// ok
		} else if (c == 'G' && !weighted) { // B+E
			// ok
		} else if (c == 'I' && !weighted) { // Inc B+E
			// ok
		} else {
			return false;
		}
	}
	return d == 0;
}
} // namespace sspp