#include <istream>
#include <fstream>
#include <cctype>
#include <iomanip>

#include "instance.hpp"
#include "utils.hpp"

namespace sspp {
namespace {
vector<string> Tokens(string t) {
	vector<string> ret;
	ret.push_back("");
	for (char c : t) {
		if (std::isspace(c)) {
			if (!ret.back().empty()) {
				ret.push_back("");
			}
		} else {
			ret.back() += c;
		}
	}
	if (ret.back().empty()) ret.pop_back();
	return ret;
}
} // namespace

Var Instance::AddVar() {
	return ++vars;
}

void Instance::AddClause(vector<Lit> clause) {
	assert(!clause.empty());
	for (Lit l : clause) {
		assert(1 <= VarOf(l) && VarOf(l) <= vars);
	}
	SortAndDedup(clause);
	for (size_t i = 1; i < clause.size(); i++) {
		if (VarOf(clause[i-1]) == VarOf(clause[i])) {
			// Tautology
			return;
		}
	}
	clauses.push_back(clause);
}

void Instance::AddLearnedClause(vector<Lit> clause) {
	assert(!clause.empty());
	for (Lit l : clause) {
		assert(1 <= VarOf(l) && VarOf(l) <= vars);
	}
	SortAndDedup(clause);
	for (size_t i = 1; i < clause.size(); i++) {
		if (VarOf(clause[i-1]) == VarOf(clause[i])) {
			// Tautology
			return;
		}
	}
	learned_clauses.push_back(clause);
}

void Instance::UpdClauseInfo() {
	total_lits = 0;
	unit_clauses = 0;
	binary_clauses = 0;
	for (const auto& clause : clauses) {
		assert(clause.size() > 0);
		if (clause.size() == 1) {
			unit_clauses++;
		} else if (clause.size() == 2) {
			binary_clauses++;
		}
		total_lits += clause.size();
		if (total_lits >= (1<<30)) {
			// Do not support so large formulas
			assert(0);
		}
	}
}

Instance::Instance(int vars_, vector<vector<Lit>> clauses_) : vars(vars_), clauses(clauses_) {
	UpdClauseInfo();
	for (const auto& clause : clauses) {
		for (int i = 0; i < (int)clause.size(); i++) {
			Var v = VarOf(clause[i]);
			assert(v >= 1 && v <= vars);
			if (i > 0) {
				assert(v > VarOf(clause[i-1]));
			}
		}
	}
}

Instance::Instance(int vars_) : vars(vars_) {}

Instance::Instance(string input_file, bool weighted_) {
	weighted = weighted_;
	std::ifstream in(input_file);
	string tmp;
	vector<Lit> cur_clause;
	int pline_clauses = 0;
	int read_clauses = 0;
	int format = 0;
	int read_weights = 0;
	while (std::getline(in, tmp)) {
		if (tmp.empty()) continue;
		auto tokens = Tokens(tmp);
		if (tokens.empty()) {
			continue;
		} else if (weighted && format == 1 && tokens.size() == 6 && tokens[0] == "c" && tokens[1] == "p" && tokens[2] == "weight") {
			assert(IsInt(tokens[3], -vars, vars));
			int dlit = stoi(tokens[3]);
			double w = stod(tokens[4]);
			assert(dlit != 0);
			Lit lit = FromDimacs(dlit);
			weights[lit] = w;
			weights[Neg(lit)] = (double)1-w;
			read_weights++;
		} else if (tokens[0] == "c") {
			continue;
		} else if (format == 0 && tokens.size() == 4 && tokens[0] == "p" && tokens[1] == "cnf") {
			format = 1;
			vars = stoi(tokens[2]);
			pline_clauses = stoi(tokens[3]);
			if (weighted) {
				weights.resize(vars*2+2);
				for (int i = 0; i < (int)vars*2+2; i++) {
					weights[i] = 1;
				}
			}
		} else if (format == 1 && IsInt(tokens[0])) {
			for (string& t : tokens) {
				assert(IsInt(t, -vars, vars));
				int dlit = stoi(t);
				if (dlit == 0) {
					read_clauses++;
					AddClause(cur_clause);
					cur_clause.clear();
				} else {
					Lit lit = FromDimacs(dlit);
					cur_clause.push_back(lit);
				}
			}
		}
	}
	if (weighted) {
		for (int v = 1; v <= vars; v++) {
			if (weights[PosLit(v)] == 0) {
				AddClause({NegLit(v)});
			}
			if (weights[NegLit(v)] == 0) {
				AddClause({PosLit(v)});
			}
		}
	}
	assert(format == 1);
	assert(cur_clause.empty());
	if (pline_clauses != read_clauses) {
		cout<<"c o Warning: p line mismatch. Claimed clauses: "<<pline_clauses<<" actual clauses: "<<read_clauses<<endl;
	}
	cout<<"c o Parsed "<<vars<<" vars, "<<clauses.size()<<" clauses, and "<<read_weights<<" weights."<<endl;
	UpdClauseInfo();
}

void Instance::PrintInfo() const {
	cerr<<"Vars: "<<vars<<endl;
	cerr<<"Clauses: "<<clauses.size()<<endl;
	cerr<<"Total literals: "<<total_lits<<endl;
	cerr<<"Unit clauses: "<<unit_clauses<<endl;
	cerr<<"Binary clauses: "<<binary_clauses<<endl;
	cerr<<endl;
}

void Instance::Print(std::ostream& out) const {
	out<<"p cnf "<<vars<<" "<<clauses.size()<<endl;
	for (const auto& clause : clauses) {
		for (Lit lit : clause) {
			out<<ToDimacs(lit)<<" ";
		}
		out<<0<<endl;
	}
}

void Instance::Eliminate(Var var) {
	vector<vector<Lit>> pos;
	vector<vector<Lit>> neg;
	for (int i = 0; i < (int)clauses.size(); i++) {
		bool found = false;
		for (Lit lit : clauses[i]) {
			if (VarOf(lit) == var) {
				vector<Lit> others;
				for (Lit l : clauses[i]) {
					if (VarOf(l) != var) {
						others.push_back(l);
					}
				}
				if (IsPos(lit)) {
					pos.push_back(others);
				} else {
					neg.push_back(others);
				}
				found = true;
				break;
			}
		}
		if (found) {
			swap(clauses.back(), clauses[i]);
			clauses.pop_back();
			i--;
		}
	}
	cerr<<"elim "<<pos.size()<<" "<<neg.size()<<endl;
	vector<vector<Lit>> cadd;
	for (const auto& c1 : pos) {
		for (const auto& c2 : neg) {
			vector<Lit> clause;
			for (Lit lit : c1) {
				clause.push_back(lit);
			}
			for (Lit lit : c2) {
				clause.push_back(lit);
			}
			SortAndDedup(clause);
			bool taut = false;
			for (int i=1;i<(int)clause.size();i++){
				if (VarOf(clause[i]) == VarOf(clause[i-1])){
					taut = true;
					break;
				}
			}
			if (!taut) {
				cadd.push_back(clause);
			}
		}
	}
	// Remove subsumed clauses
	std::sort(cadd.begin(), cadd.end(), [](const vector<Lit>& a, const vector<Lit>& b){
		if (a.size() == b.size()) {
			return a < b;
		}
		return a.size() < b.size();
	});
	{
		int subsumed = 0;
		vector<vector<Lit>> tclauses;
		for (int i = 0; i < (int)cadd.size(); i++) {
			bool sfo = false;
			for (int j = 0; j < i; j++) {
				if (Subsumes(cadd[j], cadd[i])) {
					sfo = true;
					break;
				}
			}
			if (!sfo) {
				tclauses.push_back(cadd[i]);
			} else {
				subsumed++;
			}
		}
		cerr<<"subsumed "<<subsumed<<endl;
		for (auto cl : tclauses) {
			AddClause(cl);
		}
	}
}
} // namespace sspp