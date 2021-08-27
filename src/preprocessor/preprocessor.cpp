#include "preprocessor.hpp"

#include "subsumer.hpp"
#include "oracle.hpp"
#include "treewidth.hpp"
#include "twpp.hpp"

namespace sspp {
using oracle::Oracle;

Instance Preprocessor::UnsatInst() {
	cout<<"c o UNSAT"<<endl;
	Instance ins(1);
	ins.AddClause({NegLit(1)});
	ins.AddClause({PosLit(1)});
	if (weighted) {
		ins.weighted = true;
		ins.weights = {0, 0, 0.5, 0.5};
	}
	return ins;
}

void Preprocessor::SetMaxGTime(double time) {
	max_g_time = time;
}

void Preprocessor::SetMaxSparsTime(double time) {
	max_s_time = time;
}

void Preprocessor::MapClauses(vector<vector<Lit>>& cls, Var& nvars, vector<Var>& map_to) {
	for (int i = 0; i < (int)cls.size(); i++) {
		auto& clause = cls[i];
		bool sat = false;
		for (int j = 0; j < (int)clause.size(); j++) {
			Var v = var_map[VarOf(clause[j])];
			if (assign[v] == 2 && IsPos(clause[j])) {
				SwapDel(clause, j);
				j--;
				continue;
			} else if (assign[v] == 1 && IsNeg(clause[j])) {
				SwapDel(clause, j);
				j--;
				continue;
			} else if (assign[v] != 0) {
				sat = true;
				break;
			}
		}
		if (sat) {
			SwapDel(cls, i);
			i--;
			continue;
		}
		if (clause.size() == 0) {
			unsat = true;
		}
		for (Lit& l : clause) {
			Var v = VarOf(l);
			assert(assign[var_map[v]] == 0);
			if (map_to[v] == 0) {
				map_to[v] = ++nvars;
			}
			l = MkLit(map_to[v], IsPos(l));
		}
		sort(clause.begin(), clause.end());
	}
}

void Preprocessor::TakeUnits(vector<vector<Lit>>& cls) {
	for (int i = 0; i < (int)cls.size(); i++) {
		const auto& clause = cls[i];
		for (int j = 1; j < (int)clause.size(); j++) {
			assert(VarOf(clause[j]) > VarOf(clause[j-1]));
		}
		if (clause.size() == 0) {
			unsat = true;
			return;
		}
		if (clause.size() == 1) {
			Var v = var_map[VarOf(clause[0])];
			if (assign[v] == 1) {
				if (IsNeg(clause[0])) {
					unsat = true;
					return;
				}
			} else if (assign[v] == 2) {
				if (IsPos(clause[0])) {
					unsat = true;
					return;
				}
			} else {
				if (IsPos(clause[0])) {
					assign[v] = 1;
				} else {
					assign[v] = 2;
				}
			}
			SwapDel(cls, i);
			i--;
			continue;
		}
	}
}

void Preprocessor::Tighten(bool loop) {
	SortAndDedup(clauses);
	SortAndDedup(learned_clauses);
	TakeUnits(clauses);
	TakeUnits(learned_clauses);
	if (unsat) return;
	Var nvars = 0;
	vector<Var> map_to(vars+1);
	MapClauses(clauses, nvars, map_to);
	MapClauses(learned_clauses, nvars, map_to);
	if (unsat) return;
	assert(nvars <= vars);
	vector<Var> new_var_map(nvars+1);
	for (Var v = 1; v <= vars; v++){
		if (map_to[v]) {
			assert(map_to[v] <= nvars);
			new_var_map[map_to[v]] = var_map[v];
		}
	}
	var_map = new_var_map;
	vars = nvars;
	bool unit = false;
	for (const auto& clause : clauses) {
		assert(IsClause(clause));
		if (clause.size() == 1) {
			unit = true;
		}
		for (Lit lit : clause) {
			assert(lit >= 2 && VarOf(lit) <= vars);
		}
	}
	for (const auto& clause : learned_clauses) {
		assert(IsClause(clause));
		if (clause.size() == 1) {
			unit = true;
		}
		for (Lit lit : clause) {
			assert(lit >= 2 && VarOf(lit) <= vars);
		}
	}
	if (unit && loop) {
		Tighten(loop);
	}
}

Instance Preprocessor::Preprocess(Instance inst, const string& techniques) {
	weighted = inst.weighted;
	weights = inst.weights;
	return Preprocess(inst.vars, inst.clauses, techniques);
}

void Preprocessor::FailedLiterals() {
	Oracle oracle(vars, clauses, learned_clauses);
	for (Lit lit = 2; lit <= vars*2+1; lit++) {
		if (oracle.FalseByProp({lit})) {
			oracle.FreezeUnit(Neg(lit));
			clauses.push_back({Neg(lit)});
		}
	}
	Tighten(false);
	Subsume();
}

void Preprocessor::Subsume() {
	Subsumer subsumer;
	clauses = subsumer.Subsume(clauses);
	SortAndDedup(clauses);
	if (learned_clauses.empty()) {
		return;
	}
	for (const auto& clause : clauses) {
		learned_clauses.push_back(clause);
	}
	Subsumer subsumer2;
	learned_clauses = subsumer2.Subsume(learned_clauses);
	SortAndDedup(learned_clauses);
	assert(IsSorted(clauses));
	assert(IsSorted(learned_clauses));
	for (int i = 0; i < (int)learned_clauses.size(); i++) {
		if (BS(clauses, learned_clauses[i])) {
			SwapDel(learned_clauses, i);
			i--;
		}
	}
}

int Preprocessor::FreeVars() const {
	return free_vars;
}

void Preprocessor::PropStren() {
	Oracle oracle(vars, clauses, learned_clauses);
	bool found = true;
	while (found) {
		found = false;
		for (int i = 0; i < (int)clauses.size(); i++) {
			for (int j = 0; j < (int)clauses[i].size(); j++) {
				auto assump = Negate(clauses[i]);
				SwapDel(assump, j);
				if (oracle.FalseByProp(assump)) {
					found = true;
					sort(assump.begin(), assump.end());
					auto clause = Negate(assump);
					oracle.AddClauseIfNeeded(clause, true);
					clauses[i] = clause;
					j = -1;
				}
			}
		}
	}
	Tighten(false);
	Subsume();
}

void Preprocessor::BackBone() {
	Oracle oracle(vars, clauses, learned_clauses);
	bool sat = false;
	for (int i = 0; i < (int)clauses.size(); i++) {
		for (int j = 0; j < (int)clauses[i].size(); j++) {
			auto assump = Negate(clauses[i]);
			SwapDel(assump, j);
			if (!oracle.Solve(assump)) {
				sort(assump.begin(), assump.end());
				auto clause = Negate(assump);
				oracle.AddClauseIfNeeded(clause, true);
				clauses[i] = clause;
				j = -1;
				if (clause.empty()) {
					unsat = true;
					cout<<"c o UNSAT"<<endl;
					return;
				}
			} else if(!sat) {
				sat = true;
				cout<<"c o SAT"<<endl;
			}
		}
	}
	if (!sat) {
		if (!oracle.Solve({})) {
			unsat = true;
			cout<<"c o UNSAT"<<endl;
			return;
		} else {
			sat = true;
			cout<<"c o SAT"<<endl;
		}
	}
	for (const auto& clause : oracle.LearnedClauses()) {
		learned_clauses.push_back(clause);
	}
	Tighten(false);
	Subsume();
}

Instance Preprocessor::MapBack() {
	if (unsat) return UnsatInst();
	assert(var_map.size() == (size_t)vars+1);
	Instance ret(vars);
	for (const auto& clause : clauses) {
		assert(clause.size() >= 2);
		ret.AddClause(clause);
	}
	for (const auto& clause : learned_clauses) {
		assert(clause.size() >= 2);
		ret.AddLearnedClause(clause);
	}
	free_vars = orig_vars - vars;
	for (Var v = 1; v <= orig_vars; v++) {
		if (assign[v]) {
			free_vars--;
		}
	}
	if (weighted) {
		ret.weighted = true;
		ret.weights.resize(vars*2+2);
		ret.weight_factor = 1;
		for (Var v = 1; v <= vars; v++) {
			Var ov = var_map[v];
			assert(ov <= orig_vars && ov >= 1);
			ret.weights[PosLit(v)] = weights[PosLit(ov)];
			ret.weights[NegLit(v)] = weights[NegLit(ov)];
		}
		vector<char> used(orig_vars+1);
		for (const auto& clause : clauses) {
			for (Lit lit : clause) {
				Var v = VarOf(lit);
				assert(v <= vars && v >= 1);
				v = var_map[v];
				assert(v <= orig_vars && v >= 1);
				assert(assign[v] == 0);
				used[v] = true;
			}
		}
		for (Var v = 1; v <= orig_vars; v++) {
			if (assign[v]) {
				assert(!used[v]);
				if (assign[v] == 1) {
					ret.weight_factor *= weights[PosLit(v)];
				} else if (assign[v] == 2) {
					ret.weight_factor *= weights[NegLit(v)];
				} else {
					assert(0);
				}
			} else if (!used[v]) {
				ret.weight_factor *= (weights[PosLit(v)] + weights[NegLit(v)]);
			}
		}
	}
	ret.UpdClauseInfo();
	return ret;
}

void Preprocessor::Sparsify() {
	s_timer.start();
	vector<vector<int>> edgew(vars+1);
	for (Var i = 1; i <= vars; i++) {
		edgew[i].resize(vars+1);
	}
	for (const auto& clause : clauses) {
		for (Lit l1 : clause) {
			for (Lit l2 : clause) {
				Var v1 = VarOf(l1);
				Var v2 = VarOf(l2);
				if (v1 < v2) {
					edgew[v1][v2]++;
					assert(edgew[v1][v2] <= (int)clauses.size());
				}
			}
		}
	}
	vector<pair<vector<int>, vector<Lit>>> cs;
	for (const auto& clause : clauses) {
		vector<int> ww;
		for (Lit l1 : clause) {
			for (Lit l2 : clause) {
				Var v1 = VarOf(l1);
				Var v2 = VarOf(l2);
				if (v1 < v2) {
					assert(edgew[v1][v2] >= 1);
					if ((int)ww.size() < edgew[v1][v2]) {
						ww.resize(edgew[v1][v2]);
					}
					ww[edgew[v1][v2]-1]--;
				}
			}
		}
		cs.push_back({ww, clause});
	}
	sort(cs.begin(), cs.end());
	assert(cs.size() == clauses.size());
	for (int i = 0; i < (int)clauses.size(); i++) {
		clauses[i] = cs[i].S;
	}
	cs.clear();
	// TODO: optimize
	Oracle oracle(vars+(int)clauses.size(), {});
	for (int i = 0; i < (int)clauses.size(); i++) {
		vector<Lit> clause = clauses[i];
		clause.push_back(PosLit(vars+1+i));
		oracle.AddClause(clause, false);
	}
	for (int i = 0; i < (int)clauses.size(); i++) {
		oracle.SetAssumpLit(NegLit(vars+1+i), false);
	}
	vector<int> in;
	for (int i = 0; i < (int)clauses.size(); i++) {
		oracle.SetAssumpLit(PosLit(vars+1+i), false);
		if (s_timer.get() > max_s_time || oracle.Solve(Negate(clauses[i]), false)) {
			oracle.SetAssumpLit(NegLit(vars+1+i), true);
			in.push_back(i);
		} else {
			oracle.SetAssumpLit(PosLit(vars+1+i), true);
			learned_clauses.push_back(clauses[i]);
		}
	}
	int j = 0;
	for (int i : in) {
		clauses[j++] = clauses[i];
	}
	clauses.resize(j);
	Subsume();
	s_timer.stop();
}

void Preprocessor::eqdfs(Lit lit, Lit e, const vector<vector<Lit>>& eq, vector<Lit>& eqc) {
	if (eqc[lit]) {
		assert(eqc[lit] == e);
		return;
	}
	eqc[lit] = e;
	for (Lit nl : eq[lit]) {
		eqdfs(nl, e, eq, eqc);
	}
}

void Preprocessor::MergeAdjEquivs() {
	vector<vector<char>> pg(vars+1);
	for (Var v = 1; v <= vars; v++) {
		pg[v].resize(vars+1);
	}
	for (const auto& clause : clauses) {
		for (Lit l1 : clause) {
			for (Lit l2 : clause) {
				Var v1 = VarOf(l1);
				Var v2 = VarOf(l2);
				if (v1 < v2) {
					pg[v1][v2] = true;
				}
			}
		}
	}
	Oracle oracle(vars, clauses, learned_clauses);
	vector<vector<Lit>> eq(vars*2+2);
	for (Var v1 = 1; v1 <= vars; v1++) {
		for (Var v2 = v1+1; v2 <= vars; v2++) {
			if (!pg[v1][v2]) continue;
			if (!oracle.Solve({PosLit(v1), PosLit(v2)}) && !oracle.Solve({NegLit(v1), NegLit(v2)})) {
				eq[PosLit(v2)].push_back(NegLit(v1));
				eq[NegLit(v2)].push_back(PosLit(v1));
				eq[NegLit(v1)].push_back(PosLit(v2));
				eq[PosLit(v1)].push_back(NegLit(v2));
			} else if (!oracle.Solve({PosLit(v1), NegLit(v2)}) && !oracle.Solve({NegLit(v1), PosLit(v2)})) {
				eq[PosLit(v2)].push_back(PosLit(v1));
				eq[NegLit(v2)].push_back(NegLit(v1));
				eq[PosLit(v1)].push_back(PosLit(v2));
				eq[NegLit(v1)].push_back(NegLit(v2));
			}
		}
	}
	int eq_classes = 0;
	vector<Lit> eqc(vars*2+2);
	for (Var v = 1; v <= vars; v++) {
		if (eqc[PosLit(v)] == 0) {
			eq_classes++;
			assert(eqc[NegLit(v)] == 0);
			eqdfs(PosLit(v), PosLit(v), eq, eqc);
			eqdfs(NegLit(v), NegLit(v), eq, eqc);
		}
		assert(eqc[NegLit(v)] == Neg(eqc[PosLit(v)]));
	}
	for (int i = 0; i < (int)clauses.size(); i++) {
		for (Lit& lit : clauses[i]) {
			lit = eqc[lit];
		}
		SortAndDedup(clauses[i]);
		bool taut = false;
		for (int j = 1; j < (int)clauses[i].size(); j++) {
			if (VarOf(clauses[i][j]) == VarOf(clauses[i][j-1])) {
				taut = true;
			}
		}
		if (taut) {
			SwapDel(clauses, i);
			i--;
			continue;
		}
	}
	for (Lit lit = 2; lit <= vars*2+1; lit++) {
		assert(Neg(eqc[lit]) == eqc[Neg(lit)]);
		if (eqc[lit] != lit) {
			clauses.push_back({lit, Neg(eqc[lit])});
			SortAndDedup(clauses.back());
		}
	}
	for (const auto& clause : oracle.LearnedClauses()) {
		learned_clauses.push_back(clause);
	}
	Subsume();
}

bool Preprocessor::EliminateDefSimplicial() {
	g_timer.start();
	Graph graph(vars, clauses);
	TWPP twpp;
	graph = twpp.PP(graph);
	bool found = false;
	while (true) {
		int simps = 0;
		vector<Var> extra(vars+1);
		for (Var v = 1; v <= vars; v++) {
			if (!graph.Neighbors(v).empty() && graph.IsSimp(v)) {
				int poss = 0;
				int negs = 0;
				for (const auto& clause : clauses) {
					if (BS(clause, PosLit(v))) {
						poss++;
					}
					if (BS(clause, NegLit(v))) {
						negs++;
					}
				}
				// TODO magic constant 4
				if (min(poss, negs) <= 4) {
					simps++;
					extra[v] = vars + simps;
				}
			}
		}
		if (simps == 0) {
			Tighten(true);
			Subsume();
			if (found && g_timer.get() < max_g_time) {
				EliminateDefSimplicial();
				return true;
			} else {
				g_timer.stop();
				return false;
			}
		}
		Oracle oracle(vars + 2*simps, {});
		for (const auto& cls : {clauses, learned_clauses}) {
			for (const auto& clause : cls) {
				oracle.AddClause(clause, false);
				bool hs = false;
				for (Lit lit : clause) {
					if (extra[VarOf(lit)]) {
						hs = true;
					}
				}
				if (hs) {
					vector<Lit> ac;
					for (Lit lit : clause) {
						if (extra[VarOf(lit)]) {
							ac.push_back(MkLit(extra[VarOf(lit)], IsPos(lit)));
						} else {
							ac.push_back(lit);
						}
					}
					SortAndDedup(ac);
					oracle.AddClause(ac, false);
				}
			}
		}
		for (Var v = 1; v <= vars; v++) {
			if (extra[v]) {
				oracle.AddClause({PosLit(v), NegLit(extra[v]), NegLit(extra[v]+simps)}, false);
				oracle.AddClause({NegLit(v), PosLit(extra[v]), NegLit(extra[v]+simps)}, false);
			}
		}
		vector<char> def(vars+1);
		int defs = 0;
		for (Var v = 1; v <= vars; v++) {
			if (extra[v] && g_timer.get() < max_g_time) {
				vector<Lit> assumps = {PosLit(v), NegLit(extra[v])};
				for (Var tv = 1; tv <= vars; tv++) {
					if (extra[tv] && tv != v && !def[tv]) {
						assumps.push_back(PosLit(extra[tv]+simps));
					}
				}
				if (!oracle.Solve(assumps)) {
					def[v] = 1;
					defs++;
				}
			}
		}
		for (Var v = 1; v <= vars; v++) {
			if (def[v]) {
				auto nbs = graph.Neighbors(v);
				graph.RemoveEdgesBetween(v, nbs);
				vector<vector<Lit>> pos;
				vector<vector<Lit>> neg;
				for (int i = 0; i < (int)clauses.size(); i++) {
					for (int j = 0; j < (int)clauses[i].size(); j++) {
						if (VarOf(clauses[i][j]) == v) {
							bool ph = IsPos(clauses[i][j]);
							ShiftDel(clauses[i], j);
							assert(IsClause(clauses[i]));
							if (ph) {
								pos.push_back(clauses[i]);
							} else {
								neg.push_back(clauses[i]);
							}
							SwapDel(clauses, i);
							i--;
							break;
						}
					}
				}
				assert(min(pos.size(), neg.size()) >= 1);
				for (const auto& c1 : pos) {
					for (const auto& c2 : neg) {
						vector<Lit> res;
						bool taut = false;
						int j = 0;
						for (int i = 0; i < (int)c1.size(); i++) {
							while (j < (int)c2.size() && VarOf(c2[j]) < VarOf(c1[i])) {
								res.push_back(c2[j]);
								j++;
							}
							if (j < (int)c2.size() && VarOf(c2[j]) == VarOf(c1[i])) {
								if (c2[j] == c1[i]) {
									res.push_back(c1[i]);
									j++;
								} else {
									taut = true;
									break;
								}
							} else {
								res.push_back(c1[i]);
							}
						}
						while (j < (int)c2.size() && !taut) {
							res.push_back(c2[j]);
							j++;
						}
						if (!taut) {
							assert(IsClause(res));
							clauses.push_back(res);
						}
					}
				}
				clauses.push_back({PosLit(v)});
			}
		}
		for (int i = 0; i < (int)learned_clauses.size(); i++) {
			bool fo = false;
			for (Lit lit : learned_clauses[i]) {
				if (def[VarOf(lit)]) {
					fo = true;
					break;
				}
			}
			if (fo) {
				SwapDel(learned_clauses, i);
				i--;
			}
		}
		Subsume();
		if (!defs) {
			Tighten(true);
			if (found && g_timer.get() < max_g_time) {
				EliminateDefSimplicial();
				return true;
			} else {
				g_timer.stop();
				return false;
			}
		}
		found = true;
	}
	g_timer.stop();
}

bool Preprocessor::DoTechniques(const string& techniques, int l, int r) {
	if (unsat || vars == 0) return true;
	if (l > r) return true;
	bool fixpoint = true;
	if (l == r) {
		if (techniques[l] == 'F') {
			FailedLiterals();
		} else if (techniques[l] == 'P') {
			PropStren();
		} else if (techniques[l] == 'V') {
			BackBone();
		} else if (techniques[l] == 'S') {
			Sparsify();
		} else if (techniques[l] == 'E') {
			MergeAdjEquivs();
		} else if (techniques[l] == 'G') {
			bool s = EliminateDefSimplicial();
			if (s) {
				fixpoint = false;
			}
		} else {
			assert(false);
		}
		cout<<"c o PP: "<<techniques[l]<<" "<<vars<<" "<<clauses.size()<<" "<<learned_clauses.size()<<" "<<timer.get()<<endl;
	} else if (techniques[l] == '[' && techniques[r] == ']') {
		fixpoint = false;
		while (!fixpoint) {
			fixpoint = DoTechniques(techniques, l+1, r-1);
		}
	} else if (techniques[l] == '[') {
		int d = 1;
		int f = -1;
		for (int j = l+1; j <= r; j++) {
			if (techniques[j] == '[') {
				d++;
			} else if (techniques[j] == ']') {
				d--;
				if (d == 0) {
					f = j;
					break;
				}
			}
		}
		assert(f >= l && f < r && techniques[f] == ']');
		fixpoint = DoTechniques(techniques, l, f);
		fixpoint &= DoTechniques(techniques, f+1, r);
	} else {
		fixpoint = DoTechniques(techniques, l, l);
		fixpoint &= DoTechniques(techniques, l+1, r);
	}
	return fixpoint;
}

Instance Preprocessor::Preprocess(int vars_, vector<vector<Lit>> clauses_, string techniques) {
	if (techniques.empty() || techniques[0] != 'F') {
		techniques = "F" + techniques;
	}
	assert(ValidTechniques(techniques, weighted));
	assert(orig_vars == 0);
	vars = vars_;
	clauses = clauses_;
	learned_clauses.clear();
	timer.start();
	orig_vars = vars;
	var_map.resize(vars+1);
	assign.resize(vars+1);
	for (Var v = 1; v <= vars; v++) {
		var_map[v] = v;
	}
	Tighten(false);
	cout<<"c o Init PP "<<vars<<" "<<clauses.size()<<" "<<timer.get()<<endl;
	DoTechniques(techniques, 0, (int)techniques.size()-1);
	return MapBack();
}
} // namespace sspp