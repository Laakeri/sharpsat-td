#pragma once

#include <vector>
#include <ostream>

#include "utils.hpp"
#include "graph.hpp"
#include "../mpfr/mpreal.h"

namespace sspp {
struct Instance {
	Instance(string input_file, bool weighted_);
	Instance(int vars_);
	Instance(int vars_, vector<vector<Lit>> clauses_);

	Var AddVar();
	void AddClause(vector<Lit> clause);
	void AddLearnedClause(vector<Lit> clause);
	void UpdClauseInfo();

	void Eliminate(Var v);

	//Graph PrimalGraph() const;
	void PrintInfo() const;

	void Print(std::ostream& out) const;

	int vars = 0;
	vector<vector<Lit>> clauses;
	vector<vector<Lit>> learned_clauses;

	bool weighted = false;
	mpfr::mpreal weight_factor = 1;

  vector<double> weights;
  int total_lits = 0;
 private:
 	Lit RecConstruct(vector<Lit> clause);
 	int format = 0;
 	int unit_clauses = 0;
	int binary_clauses = 0;
};
} // namespace sspp