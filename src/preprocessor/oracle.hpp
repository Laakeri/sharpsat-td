#pragma once

#include "instance.hpp"
#include "utils.hpp"
#include "graph.hpp"

namespace sspp {
namespace oracle {
struct Stats {
 	int64_t propagations = 0;
 	int64_t decisions = 0;
 	int64_t learned_clauses = 0;
 	int64_t learned_bin_clauses = 0;
 	int64_t learned_units = 0;
 	int64_t conflicts = 0;
 	int64_t nontriv_redu = 0;
 	int64_t forgot_clauses = 0;
 	int64_t restarts = 0;
 	Timer setup_timer, solve_timer, prop_timer, learn_timer, maint_time, cache_timer;
  void Print() const;
};

struct Watch {
	// should align to 8+4+4=16 bytes
	size_t cls;
	Lit blit;
	int size;
};

struct VarC {
	size_t reason = 0;
	int level = 0;
	char phase = 0;
};

struct CInfo {
	size_t pt;
	int glue;
	int used;
	bool Keep() const;
};

class Oracle {
 public:
 	Oracle(const Instance& inst);
 	Oracle(int vars_, const vector<vector<Lit>>& clauses_);
  Oracle(int vars_, const vector<vector<Lit>>& clauses_, const vector<vector<Lit>>& learned_clauses_);

  void SetAssumpLit(Lit lit, bool freeze);
  bool Solve(const vector<Lit>& assumps, bool usecache=true);
  bool FalseByProp(const vector<Lit>& assumps);
  bool FreezeUnit(Lit unit);
  bool FreezeUnits(const vector<Lit>& units);
  bool AddClauseIfNeeded(vector<Lit> clause, bool entailed);
  void AddClause(const vector<Lit>& clause, bool entailed);
  void PrintStats() const;
  double ConflictRate(int samples);
  vector<vector<Lit>> AllClauses() const;
  vector<vector<Lit>> LearnedClauses() const;
  bool RemoveIfPossible(const vector<Lit>& clause);

  vector<Lit> InferUnits(const vector<Lit>& assumps);
  int PropDg(const vector<Lit>& assumps);

  int CurLevel() const;
  int LitVal(Lit lit) const;
 private:
 	void ForgetLearned();
 	void AddOrigClause(vector<Lit> clause, bool entailed);
 	vector<Lit> clauses;
 	vector<vector<Watch>> watches;
 	vector<signed char> lit_val;
 	vector<VarC> vs;

 	bool unsat = false;

 	const int vars;
 	size_t orig_clauses_size = 0;
 	Stats stats;
 	vector<Lit> prop_q;
 	vector<Var> decided;
 	vector<char> in_cc;

 	std::mt19937 rand_gen;

 	int64_t og_bump = 0;

 	long double glue_long_ema = 0;
 	long double glue_short_ema = 0;
 	long double long_a = 1;
 	long double short_a = 1;
 	void UpdGlueEma(int glue);

 	long double var_ass_ema = 0;
 	long double var_ass_a = 1;
 	void UpdVarAssEma();

 	int64_t redu_it = 1;
 	vector<char> seen;
 	vector<int64_t> redu_seen;
 	vector<Lit> redu_s;

 	int64_t lvl_it = 1;
 	vector<int64_t> lvl_seen;

 	vector<Lit> learned_units;

 	int restart_factor = 0;
 	vector<int> luby;
 	void InitLuby();
 	int NextLuby();

 	size_t ideal_clause_db_size = 0;
 	void ResizeClauseDb();
 	void BumpClause(size_t cls);
 	vector<CInfo> cla_info;

 	double var_inc = 1;
 	double var_fact = 0;
 	size_t heap_N;
 	vector<double> var_act_heap;
 	void BumpVar(Var v);
 	void ActivateActivity(Var v);
 	Var PopVarHeap();

 	bool LitSat(Lit lit) const;
 	bool LitAssigned(Lit lit) const;
 	void Assign(Lit dec, size_t reason_clause, int level);
 	void Decide(Lit dec, int level);
 	void UnDecide(int level);

 	bool HardSolve();
 	// True if conflict
 	size_t Propagate(int level);

 	size_t AddLearnedClause(const vector<Lit>& clause);
 	bool LitReduntant(Lit lit);
 	void PopLit(vector<Lit>& clause, int& confl_levels, int& impl_lits, int level);
 	vector<Lit> LearnUip(size_t conflict_clause);
 	int CDCLBT(size_t confl_clause, int min_level=0);

 	vector<vector<char>> sol_cache;
 	void AddSolToCache();
 	bool SatByCache(const vector<Lit>& assumps);
 	void ClearSolCache();
};
} // namespace oracle
} // namespace sspp