/*
 * solver.h
 *
 *  Created on: Aug 23, 2012
 *      Author: marc
 */

#ifndef SOLVER_H_
#define SOLVER_H_


#include "statistics.h"
#include "instance.h"
#include "component_management.h"

#include "preprocessor/instance.hpp"
#include "preprocessor/graph.hpp"

#include "solver_config.h"

#include <sys/time.h>

#include "decomposition.h"
#include <deque>

enum retStateT {
	EXIT, RESOLVED, PROCESS_COMPONENT, BACKTRACK
};

class StopWatch {
public:

  StopWatch() {
	  interval_length_.tv_sec = 60;
	  gettimeofday(&last_interval_start_, NULL);
	  start_time_ = stop_time_ = last_interval_start_;
  }

  bool start() {
    bool ret = gettimeofday(&last_interval_start_, NULL);
    start_time_ = stop_time_ = last_interval_start_;
    return !ret;
  }

  bool stop() {
    return gettimeofday(&stop_time_, NULL) == 0;
  }

  double getElapsedSeconds() {
    timeval r = getElapsedTime();
    return r.tv_sec + (double) r.tv_usec / 1000000;
  }

  bool interval_tick() {
    timeval actual_time;
    gettimeofday(&actual_time, NULL);
    if (actual_time.tv_sec - last_interval_start_.tv_sec
        > interval_length_.tv_sec) {
      gettimeofday(&last_interval_start_, NULL);
      return true;
    }
    return false;
  }

private:
  timeval start_time_;
  timeval stop_time_;

  timeval interval_length_;
  timeval last_interval_start_;

  // if we have started and then stopped the watch, this returns
  // the elapsed time
  // otherwise, time elapsed from start_time_ till now is returned
  timeval getElapsedTime() {
	  timeval r;
	  timeval other_time = stop_time_;
	  if (stop_time_.tv_sec == start_time_.tv_sec
	      && stop_time_.tv_usec == start_time_.tv_usec)
	    gettimeofday(&other_time, NULL);
	  long int ad = 0;
	  long int bd = 0;

	  if (other_time.tv_usec < start_time_.tv_usec) {
	    ad = 1;
	    bd = 1000000;
	  }
	  r.tv_sec = other_time.tv_sec - ad - start_time_.tv_sec;
	  r.tv_usec = other_time.tv_usec + bd - start_time_.tv_usec;
	  return r;
  }
};

template<class T_num>
class Solver: public Instance<T_num> {
public:
	Solver(std::mt19937_64& gen) : hasher_(gen) {}

	T_num solve(const sspp::Instance& pp_ins, const sspp::TreeDecomposition& tdec);

	SolverConfiguration &config() {
		return config_;
	}

	DataAndStatistics<T_num> &statistics() {
	        return Instance<T_num>::statistics_;
	}

private:
	Hasher hasher_;

	SolverConfiguration config_;

	DecisionStack<T_num> stack_; // decision stack
	vector<LiteralID> literal_stack_;

	StopWatch stopwatch_;

	ComponentManager<T_num> comp_manager_ = ComponentManager<T_num>(config_,
			Instance<T_num>::statistics_, Instance<T_num>::literal_values_, Instance<T_num>::lit_weights_);

	int clause_budjet = 10000;
	unsigned last_ccl_cleanup_time_ = 0;

	unsigned last_del_decs_ = 0;
	unsigned last_del_confls_ = 0;

	bool simplePreProcess();
	bool prepFailedLiteralTest();
	// we assert that the formula is consistent
	// and has not been found UNSAT yet
	// hard wires all assertions in the literal stack into the formula
	// removes all set variables and essentially reinitiallizes all
	// further data
	void HardWireAndCompact();

	SOLVER_StateT countSAT();

	void decideLiteral();
	bool bcp();


	void decayActivitiesOf(Component & comp) {
		for (auto it = comp.varsBegin(); *it != varsSENTINEL; it++) {
	  	Instance<T_num>::literal(LiteralID(*it,true)).activity_score_ *=0.5;
	    Instance<T_num>::literal(LiteralID(*it,false)).activity_score_ *=0.5;
	  }
	}
	///  this method performs Failed literal tests online
	bool implicitBCP();

	// this is the actual BCP algorithm
	// starts propagating all literal in literal_stack_
	// beginingg at offset start_at_stack_ofs
	bool BCP(unsigned start_at_stack_ofs);

	retStateT backtrack();

	// if on the current decision level
	// a second branch can be visited, RESOLVED is returned
	// otherwise returns BACKTRACK
	retStateT resolveConflict();

	/////////////////////////////////////////////
	//  BEGIN small helper functions
	/////////////////////////////////////////////

	float scoreOf(VariableIndex v) {
		float score = 0;
		if (config_.vsads_freq) {
			score += comp_manager_.scoreOf(v);
		}
		if (config_.vsads_act) {
			score += 10.0 * Instance<T_num>::literal(LiteralID(v, true)).activity_score_;
			score += 10.0 * Instance<T_num>::literal(LiteralID(v, false)).activity_score_;
		}
		assert(v >= 0 && v < Instance<T_num>::extra_score.size());
		score += Instance<T_num>::extra_score[v];
		return score;
	}

	bool setLiteralIfFreeBase(LiteralID lit,
			Antecedent ant) {
		if (Instance<T_num>::literal_values_[lit] != X_TRI)
			return false;
		Instance<T_num>::var(lit).decision_level = stack_.get_decision_level();
		Instance<T_num>::var(lit).ante = ant;
		literal_stack_.push_back(lit);
		if (ant.isAClause() && ant.asCl() != NOT_A_CLAUSE) {
			Instance<T_num>::setGlue(ant.asCl());
		}
		Instance<T_num>::literal_values_[lit] = T_TRI;
		Instance<T_num>::literal_values_[lit.neg()] = F_TRI;
		return true;
	}

	bool setLiteralIfFree(LiteralID lit,
			Antecedent ant = Antecedent(NOT_A_CLAUSE));

	void unsetLiteral(LiteralID lit);

	void printOnlineStats();

	void print(vector<LiteralID> &vec);
	void print(vector<unsigned> &vec);


	void setConflictState(LiteralID litA, LiteralID litB) {
		violated_clause.clear();
		violated_clause.push_back(litA);
		violated_clause.push_back(litB);
	}
	void setConflictState(ClauseOfs cl_ofs) {
		Instance<T_num>::setGlue(cl_ofs);
		violated_clause.clear();
		for (auto it = Instance<T_num>::beginOf(cl_ofs); *it != SENTINEL_LIT; it++)
			violated_clause.push_back(*it);
	}

	vector<LiteralID>::const_iterator TOSLiteralsBegin() {
		return literal_stack_.begin() + stack_.top().literal_stack_ofs();
	}

	void initStack(unsigned int resSize) {
		stack_.clear();
		stack_.reserve(resSize);
		literal_stack_.clear();
		literal_stack_.reserve(resSize);
		// initialize the stack to contain at least level zero
		stack_.push_back(StackLevel<T_num>(1, 0, 2));
		stack_.back().changeBranch();
	}

	const LiteralID &TOS_decLit() {
		assert(stack_.top().literal_stack_ofs() < literal_stack_.size());
		return literal_stack_[stack_.top().literal_stack_ofs()];
	}

	void reactivateTOS() {
		for (auto it = TOSLiteralsBegin(); it != literal_stack_.end(); it++)
			unsetLiteral(*it);
		comp_manager_.cleanRemainingComponentsOf(stack_.top());
		literal_stack_.resize(stack_.top().literal_stack_ofs());
		stack_.top().resetRemainingComps();
	}

	bool fail_test(LiteralID lit) {
		unsigned sz = literal_stack_.size();
		// we increase the decLev artificially
		// s.t. after the tentative BCP call, we can learn a conflict clause
		// relative to the assignment of *jt
		stack_.startFailedLitTest();
		setLiteralIfFree(lit);

		assert(!Instance<T_num>::hasAntecedent(lit));

		bool bSucceeded = BCP(sz);
		if (!bSucceeded)
			recordAllUIPCauses();

		stack_.stopFailedLitTest();

		while (literal_stack_.size() > sz) {
			unsetLiteral(literal_stack_.back());
			literal_stack_.pop_back();
		}
		return bSucceeded;
	}
	/////////////////////////////////////////////
	//  BEGIN conflict analysis
	/////////////////////////////////////////////

	// if the state name is CONFLICT,
	// then violated_clause contains the clause determining the conflict;
	vector<LiteralID> violated_clause;
	// this is an array of all the clauses found
	// during the most recent conflict analysis
	// it might contain more than 2 clauses
	// but always will have:
	//      uip_clauses_.front() the 1UIP clause found
	//      uip_clauses_.back() the lastUIP clause found
	//  possible clauses in between will be other UIP clauses
	vector<vector<LiteralID> > uip_clauses_;

	// the assertion level of uip_clauses_.back()
	// or (if the decision variable did not have an antecedent
	// before) then assertionLevel_ == DL;
	int assertion_level_ = 0;

	// build conflict clauses from most recent conflict
	// as stored in state_.violated_clause
	// solver state must be CONFLICT to work;
	// this first method record only the last UIP clause
	// so as to create clause that asserts the current decision
	// literal
	void recordLastUIPCauses();
	void recordAllUIPCauses();

	void minimizeAndStoreUIPClause(LiteralID uipLit,
			vector<LiteralID> & tmp_clause, bool seen[]);
	void storeUIPClause(LiteralID uipLit, vector<LiteralID> & tmp_clause);
	int getAssertionLevel() const {
		return assertion_level_;
	}

	bool Weighted() const {
		return !Instance<T_num>::lit_weights_.empty();
	}

	/////////////////////////////////////////////
	//  END conflict analysis
	/////////////////////////////////////////////
};

template<typename T_num>
inline bool Solver<T_num>::setLiteralIfFree(LiteralID lit,
		Antecedent ant) {
	if (Instance<T_num>::literal_values_[lit] != X_TRI)
		return false;
	if (Weighted()) {
		if (stack_.size() >= 2 && Instance<T_num>::lit_weights_.size() > 0) {
			if (Instance<T_num>::dec_cands_[stack_.size()][lit.var()] == Instance<T_num>::dec_cands_[stack_.size()][0]) {
				Instance<T_num>::lit_mul_[lit.var()] = true;
				stack_.top().dec_weight *= Instance<T_num>::lit_weights_[lit];
			}
		}
	}
	return setLiteralIfFreeBase(lit, ant);
}

template<typename T_num>
inline void Solver<T_num>::unsetLiteral(LiteralID lit) {
	if (Weighted()) {
		if (stack_.size() >= 2 && Instance<T_num>::lit_mul_[lit.var()]) {
			stack_.top().dec_weight /= Instance<T_num>::lit_weights_[lit];
		}
		Instance<T_num>::lit_mul_[lit.var()] = false;
	}
	Instance<T_num>::unSet(lit);
}

template <class T_num>
void Solver<T_num>::print(vector<LiteralID> &vec) {
	for (auto l : vec)
		cout << l.toInt() << " ";
	cout << endl;
}

template <class T_num>
void Solver<T_num>::print(vector<unsigned> &vec) {
	for (auto l : vec)
		cout << l << " ";
	cout << endl;
}

template <class T_num>
bool Solver<T_num>::simplePreProcess() {
	assert(literal_stack_.size() == 0);
	unsigned start_ofs = 0;
//BEGIN process unit clauses
	for (auto lit : Instance<T_num>::unit_clauses_)
		setLiteralIfFree(lit);
//END process unit clauses
	bool succeeded = BCP(start_ofs);

	if (succeeded)
		succeeded &= prepFailedLiteralTest();

	if (succeeded)
		HardWireAndCompact();
	return succeeded;
}

template <class T_num>
bool Solver<T_num>::prepFailedLiteralTest() {
	unsigned last_size;
	do {
		last_size = literal_stack_.size();
		for (unsigned v = 1; v < Instance<T_num>::variables_.size(); v++)
			if (Instance<T_num>::isActive(v)) {
				unsigned sz = literal_stack_.size();
				setLiteralIfFree(LiteralID(v, true));
				bool res = BCP(sz);
				while (literal_stack_.size() > sz) {
					unsetLiteral(literal_stack_.back());
					literal_stack_.pop_back();
				}

				if (!res) {
					sz = literal_stack_.size();
					setLiteralIfFree(LiteralID(v, false));
					if (!BCP(sz))
						return false;
				} else {

					sz = literal_stack_.size();
					setLiteralIfFree(LiteralID(v, false));
					bool resb = BCP(sz);
					while (literal_stack_.size() > sz) {
						unsetLiteral(literal_stack_.back());
						literal_stack_.pop_back();
					}
					if (!resb) {
						sz = literal_stack_.size();
						setLiteralIfFree(LiteralID(v, true));
						if (!BCP(sz))
							return false;
					}
				}
			}
	} while (literal_stack_.size() > last_size);

	return true;
}

template <class T_num>
void Solver<T_num>::HardWireAndCompact() {
	Instance<T_num>::compactClauses();
	//Instance<T_num>::compactVariables();
	literal_stack_.clear();

	for (auto l = LiteralID(1, false); l != Instance<T_num>::literals_.end_lit(); l.inc()) {
		if (config_.vsads_act_init) {
			Instance<T_num>::literal(l).activity_score_ = Instance<T_num>::literal(l).binary_links_.size() - 1;
			Instance<T_num>::literal(l).activity_score_ += Instance<T_num>::occurrence_lists_[l].size();
		} else {
			Instance<T_num>::literal(l).activity_score_ = 0;
		}
	}

	Instance<T_num>::statistics_.num_unit_clauses_ = Instance<T_num>::unit_clauses_.size();

	Instance<T_num>::statistics_.num_original_binary_clauses_ = Instance<T_num>::statistics_.num_binary_clauses_;
	Instance<T_num>::statistics_.num_original_unit_clauses_ = Instance<T_num>::statistics_.num_unit_clauses_ =
			Instance<T_num>::unit_clauses_.size();
	initStack(Instance<T_num>::num_variables());
	Instance<T_num>::original_lit_pool_size_ = Instance<T_num>::literal_pool_.size();
}

template <class T_num>
T_num Solver<T_num>::solve(const sspp::Instance& pp_ins, const sspp::TreeDecomposition& tdec) {
	stopwatch_.start();
	Instance<T_num>::createfromPPIns(pp_ins);
	initStack(Instance<T_num>::num_variables());

	// Keep simplepreprocess here just to not break any invariant. Should be cheap.
	bool notfoundUNSAT = simplePreProcess();
	cout<<"c o SharpSAT loaded"<<endl;

	if (notfoundUNSAT) {

		if (!config_.quiet) {
			Instance<T_num>::statistics_.printShortFormulaInfo();
		}

		Instance<T_num>::PrepareTWScore(tdec, config_.decomp_weight, config_.weight_mode);

		violated_clause.reserve(Instance<T_num>::num_variables());

		comp_manager_.initialize(Instance<T_num>::literals_, Instance<T_num>::literal_pool_, hasher_);

		Instance<T_num>::statistics_.exit_state_ = countSAT();

		Instance<T_num>::statistics_.set_final_solution_count(stack_.top().getTotalModelCount());
		Instance<T_num>::statistics_.num_long_conflict_clauses_ = Instance<T_num>::num_conflict_clauses();

	} else {
		Instance<T_num>::statistics_.exit_state_ = SUCCESS;
		Instance<T_num>::statistics_.set_final_solution_count(T_num::Zero());
	}
	stopwatch_.stop();
	Instance<T_num>::statistics_.time_elapsed_ = stopwatch_.getElapsedSeconds();

	comp_manager_.gatherStatistics();
	Instance<T_num>::statistics_.printShort();
	return Instance<T_num>::statistics_.final_solution_count();
}

template <class T_num>
SOLVER_StateT Solver<T_num>::countSAT() {
	retStateT state = RESOLVED;

	while (true) {
		//debugStack();
		while (comp_manager_.findNextRemainingComponentOf(stack_.top(), hasher_)) {
			//debugStack();
			decideLiteral();
			if (stopwatch_.interval_tick())
				printOnlineStats();
			//debugStack();
			while (!bcp()) {
				//debugStack();
				state = resolveConflict();
				if (state == BACKTRACK)
					break;
			}
			if (state == BACKTRACK)
				break;
		}
		//debugStack();
		state = backtrack();
		//debugStack();
		if (state == EXIT)
			return SUCCESS;
		while (state != PROCESS_COMPONENT && !bcp()) {
			//debugStack();
			state = resolveConflict();
			//debugStack();
			if (state == BACKTRACK) {
				state = backtrack();
				if (state == EXIT)
					return SUCCESS;
			}
		}
	}
	return SUCCESS;
}


template <class T_num>
void Solver<T_num>::decideLiteral() {
	// establish another decision stack level
	stack_.push_back(
			StackLevel<T_num>(stack_.top().currentRemainingComponent(),
					literal_stack_.size(),
					comp_manager_.component_stack_size()));
	float max_score = -1;
	float score;
	unsigned max_score_var = 0;
	if (Instance<T_num>::lit_weights_.size() > 0) {
		while (stack_.size() + 1 >= Instance<T_num>::dec_cands_.size()) {
			Instance<T_num>::dec_cands_.push_back({});
			Instance<T_num>::dec_cands_.back().resize(Instance<T_num>::num_variables()+1);
		}
		Instance<T_num>::dec_cands_[stack_.size()][0]++;
	}
	for (auto it =
			comp_manager_.superComponentOf(stack_.top()).varsBegin();
			*it != varsSENTINEL; it++) {
		if (Instance<T_num>::lit_weights_.size() > 0) {
			Instance<T_num>::dec_cands_[stack_.size()][*it] = Instance<T_num>::dec_cands_[stack_.size()][0];
		}
		score = scoreOf(*it);
		if (score > max_score) {
			max_score = score;
			max_score_var = *it;
		}
	}
	// this assert should always hold,
	// if not then there is a bug in the logic of countSAT();
	assert(max_score_var != 0);

	LiteralID theLit(max_score_var,
			Instance<T_num>::literal(LiteralID(max_score_var, true)).activity_score_
					> Instance<T_num>::literal(LiteralID(max_score_var, false)).activity_score_);


	setLiteralIfFree(theLit);
	Instance<T_num>::statistics_.num_decisions_++;

	if (Instance<T_num>::statistics_.num_decisions_ % 128 == 0)
//    if (statistics_.num_conflicts_ % 128 == 0)
     Instance<T_num>::decayActivities();
       // decayActivitiesOf(comp_manager_.superComponentOf(stack_.top()));
	assert(
			stack_.top().remaining_components_ofs() <= comp_manager_.component_stack_size());
}

template <class T_num>
retStateT Solver<T_num>::backtrack() {
	assert(
			stack_.top().remaining_components_ofs() <= comp_manager_.component_stack_size());
	do {
		if (stack_.top().branch_found_unsat())
			comp_manager_.removeAllCachePollutionsOf(stack_.top());
		else if (stack_.top().anotherCompProcessible())
			return PROCESS_COMPONENT;

		if (!stack_.top().isSecondBranch()) {
			LiteralID aLit = TOS_decLit();
			assert(stack_.get_decision_level() > 0);
			stack_.top().changeBranch();
			reactivateTOS();
			setLiteralIfFree(aLit.neg(), NOT_A_CLAUSE);
			return RESOLVED;
		}
		// OTHERWISE:  backtrack further
		comp_manager_.cacheModelCountOf(stack_.top().super_component(),
				stack_.top().getTotalModelCount());

		if (stack_.get_decision_level() <= 0)
			break;
		reactivateTOS();

		assert(stack_.size()>=2);
		(stack_.end() - 2)->includeSolution(stack_.top().getTotalModelCount());
		stack_.pop_back();
		// step to the next component not yet processed
		stack_.top().nextUnprocessedComponent();

		assert(
				stack_.top().remaining_components_ofs() < comp_manager_.component_stack_size()+1);

	} while (stack_.get_decision_level() >= 0);
	return EXIT;
}

template <class T_num>
retStateT Solver<T_num>::resolveConflict() {
	recordLastUIPCauses();

	if (Instance<T_num>::statistics_.num_decisions_ > 2000000
	 && Instance<T_num>::statistics_.num_conflicts_ < 1000) {
		clause_budjet = 1000;
	}
	if ((int)Instance<T_num>::conflict_clauses_.size() > clause_budjet) {
		cout << "c o cleaning clauses "<<Instance<T_num>::conflict_clauses_.size()<<" "<<clause_budjet<<endl;
		int cutoff = Instance<T_num>::deleteConflictClauses();
		assert(Instance<T_num>::statistics_.num_decisions_ > last_del_decs_);
		assert(Instance<T_num>::statistics_.num_conflicts_ > last_del_confls_);
		unsigned new_decs_ = Instance<T_num>::statistics_.num_decisions_-last_del_decs_;
		unsigned new_confls_ = Instance<T_num>::statistics_.num_conflicts_-last_del_confls_;
		int target = 10000;
		if (new_confls_*2 > new_decs_) {
			target = 200000;
		} else if (new_confls_*4 > new_decs_) {
			target = 80000;
		} else if (new_confls_*8 > new_decs_) {
			target = 35000;
		} else if (new_confls_*16 > new_decs_) {
			target = 15000;
		} else if (new_confls_*32 > new_decs_) {
			target = 8000;
		} else {
			target = 1000;
		}
		if (clause_budjet < target) {
			if (cutoff == 3) {
				clause_budjet += 1000;
			} else if (cutoff <= 5) {
				clause_budjet += 500;
			} else if (cutoff <= 7) {
				clause_budjet += 250;
			} else {
				clause_budjet += 100;
			}
		} else if (cutoff > 5) {
			clause_budjet -= 500;
		}
		if (clause_budjet <= 1000) {
			clause_budjet = 1000;
		}
		if (clause_budjet >= 200000) {
			clause_budjet = 200000;
		}
		cout << "c o cleaned clauses "<<Instance<T_num>::conflict_clauses_.size()<<" "<<clause_budjet<<" "<<target<<endl;
	}

	if (Instance<T_num>::statistics_.num_clauses_learned_ - last_ccl_cleanup_time_ > 100000) {
		Instance<T_num>::compactConflictLiteralPool();
		last_ccl_cleanup_time_ = Instance<T_num>::statistics_.num_clauses_learned_;
	}

	Instance<T_num>::statistics_.num_conflicts_++;

	assert(
			stack_.top().remaining_components_ofs() <= comp_manager_.component_stack_size());

	assert(uip_clauses_.size() == 1);

	// DEBUG
	if (uip_clauses_.back().size() == 0)
		cout << " EMPTY CLAUSE FOUND" << endl;
	// END DEBUG

	stack_.top().mark_branch_unsat();
	//BEGIN Backtracking
	// maybe the other branch had some solutions
	if (stack_.top().isSecondBranch()) {
		return BACKTRACK;
	}

	Antecedent ant(NOT_A_CLAUSE);
	// this has to be checked since using implicit BCP
	// and checking literals there not exhaustively
	// we cannot guarantee that uip_clauses_.back().front() == TOS_decLit().neg()
	// this is because we might have checked a literal
	// during implict BCP which has been a failed literal
	// due only to assignments made at lower decision levels
	if (uip_clauses_.back().front() == TOS_decLit().neg()) {
		assert(TOS_decLit().neg() == uip_clauses_.back()[0]);
		Instance<T_num>::var(TOS_decLit().neg()).ante = Instance<T_num>::addUIPConflictClause(
				uip_clauses_.back());
		ant = Instance<T_num>::var(TOS_decLit()).ante;
	}
//	// RRR
//	else if(var(uip_clauses_.back().front()).decision_level
//			< stack_.get_decision_level()
//			&& assertion_level_ <  stack_.get_decision_level()){
//         stack_.top().set_both_branches_unsat();
//         return BACKTRACK;
//	}
//
//
//	// RRR
	assert(stack_.get_decision_level() > 0);
	assert(stack_.top().branch_found_unsat());

	// we do not have to remove pollutions here,
	// since conflicts only arise directly before
	// remaining components are stored
	// hence
	assert(
			stack_.top().remaining_components_ofs() == comp_manager_.component_stack_size());

	stack_.top().changeBranch();
	LiteralID lit = TOS_decLit();
	reactivateTOS();
	setLiteralIfFree(lit.neg(), ant);
//END Backtracking
	return RESOLVED;
}

template <class T_num>
bool Solver<T_num>::bcp() {
// the asserted literal has been set, so we start
// bcp on that literal
	unsigned start_ofs = literal_stack_.size() - 1;

//BEGIN process unit clauses
	for (auto lit : Instance<T_num>::unit_clauses_)
		setLiteralIfFree(lit);
//END process unit clauses

	bool bSucceeded = BCP(start_ofs);

	return bSucceeded;
}

template <class T_num>
bool Solver<T_num>::BCP(unsigned start_at_stack_ofs) {
	for (unsigned int i = start_at_stack_ofs; i < literal_stack_.size(); i++) {
		LiteralID unLit = literal_stack_[i].neg();
		//BEGIN Propagate Bin Clauses
		for (auto bt = Instance<T_num>::literal(unLit).binary_links_.begin();
				*bt != SENTINEL_LIT; bt++) {
			if (Instance<T_num>::isResolved(*bt)) {
				setConflictState(unLit, *bt);
				return false;
			}
			setLiteralIfFree(*bt, Antecedent(unLit));
		}
		//END Propagate Bin Clauses
		for (auto itcl = Instance<T_num>::literal(unLit).watch_list_.rbegin();
				*itcl != SENTINEL_CL; itcl++) {
			bool isLitA = (*Instance<T_num>::beginOf(*itcl) == unLit);
			auto p_watchLit = Instance<T_num>::beginOf(*itcl) + 1 - isLitA;
			auto p_otherLit = Instance<T_num>::beginOf(*itcl) + isLitA;

			if (Instance<T_num>::isSatisfied(*p_otherLit))
				continue;
			auto itL = Instance<T_num>::beginOf(*itcl) + 2;
			while (Instance<T_num>::isResolved(*itL))
				itL++;
			// either we found a free or satisfied lit
			if (*itL != SENTINEL_LIT) {
				Instance<T_num>::literal(*itL).addWatchLinkTo(*itcl);
				swap(*itL, *p_watchLit);
				*itcl = Instance<T_num>::literal(unLit).watch_list_.back();
				Instance<T_num>::literal(unLit).watch_list_.pop_back();
			} else {
				// or p_unLit stays resolved
				// and we have hence no free literal left
				// for p_otherLit remain poss: Active or Resolved
				if (setLiteralIfFree(*p_otherLit, Antecedent(*itcl))) { // implication
					if (isLitA)
						swap(*p_otherLit, *p_watchLit);
				} else {
					setConflictState(*itcl);
					return false;
				}
			}
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// BEGIN module conflictAnalyzer
///////////////////////////////////////////////////////////////////////////////////////////////

template <class T_num>
void Solver<T_num>::minimizeAndStoreUIPClause(LiteralID uipLit,
		vector<LiteralID> & tmp_clause, bool seen[]) {
	static deque<LiteralID> clause;
	clause.clear();
	assertion_level_ = 0;
	for (auto lit : tmp_clause) {
		if (Instance<T_num>::existsUnitClauseOf(lit.var()))
			continue;
		bool resolve_out = false;
		if (Instance<T_num>::hasAntecedent(lit)) {
			resolve_out = true;
			if (Instance<T_num>::getAntecedent(lit).isAClause()) {
				for (auto it = Instance<T_num>::beginOf(Instance<T_num>::getAntecedent(lit).asCl()) + 1;
						*it != SENTINEL_CL; it++)
					if (!seen[it->var()]) {
						resolve_out = false;
						break;
					}
			} else if (!seen[Instance<T_num>::getAntecedent(lit).asLit().var()]) {
				resolve_out = false;
			}
		}

		if (!resolve_out) {
			// uipLit should be the sole literal of this Decision Level
			if (Instance<T_num>::var(lit).decision_level >= assertion_level_) {
				assertion_level_ = Instance<T_num>::var(lit).decision_level;
				clause.push_front(lit);
			} else
				clause.push_back(lit);
		}
	}

	if(uipLit.var())
	 assert(Instance<T_num>::var(uipLit).decision_level == stack_.get_decision_level());

	//assert(uipLit.var() != 0);
	if (uipLit.var() != 0)
		clause.push_front(uipLit);
	uip_clauses_.push_back(vector<LiteralID>(clause.begin(), clause.end()));
}

template <class T_num>
void Solver<T_num>::recordLastUIPCauses() {
// note:
// variables of lower dl: if seen we dont work with them anymore
// variables of this dl: if seen we incorporate their
// antecedent and set to unseen
	bool seen[Instance<T_num>::num_variables() + 1];
	memset(seen, false, sizeof(bool) * (Instance<T_num>::num_variables() + 1));

	static vector<LiteralID> tmp_clause;
	tmp_clause.clear();

	assertion_level_ = 0;
	uip_clauses_.clear();

	unsigned lit_stack_ofs = literal_stack_.size();
	int DL = stack_.get_decision_level();
	unsigned lits_at_current_dl = 0;

	for (auto l : violated_clause) {
		if (Instance<T_num>::var(l).decision_level == 0 || Instance<T_num>::existsUnitClauseOf(l.var()))
			continue;
		if (Instance<T_num>::var(l).decision_level < DL)
			tmp_clause.push_back(l);
		else
			lits_at_current_dl++;
		Instance<T_num>::literal(l).increaseActivity();
		seen[l.var()] = true;
	}

	LiteralID curr_lit;
	while (lits_at_current_dl) {
		assert(lit_stack_ofs != 0);
		curr_lit = literal_stack_[--lit_stack_ofs];

		if (!seen[curr_lit.var()])
			continue;

		seen[curr_lit.var()] = false;

		if (lits_at_current_dl-- == 1) {
			// perform UIP stuff
			if (!Instance<T_num>::hasAntecedent(curr_lit)) {
				// this should be the decision literal when in first branch
				// or it is a literal decided to explore in failed literal testing
				//assert(stack_.TOS_decLit() == curr_lit);
//				cout << "R" << curr_lit.toInt() << "S"
//				     << var(curr_lit).ante.isAnt() << " "  << endl;
				break;
			}
		}

		assert(Instance<T_num>::hasAntecedent(curr_lit));

		//cout << "{" << curr_lit.toInt() << "}";
		if (Instance<T_num>::getAntecedent(curr_lit).isAClause()) {
			Instance<T_num>::updateActivities(Instance<T_num>::getAntecedent(curr_lit).asCl());
			assert(curr_lit == *Instance<T_num>::beginOf(Instance<T_num>::getAntecedent(curr_lit).asCl()));

			for (auto it = Instance<T_num>::beginOf(Instance<T_num>::getAntecedent(curr_lit).asCl()) + 1;
					*it != SENTINEL_CL; it++) {
				if (seen[it->var()] || (Instance<T_num>::var(*it).decision_level == 0)
						|| Instance<T_num>::existsUnitClauseOf(it->var()))
					continue;
				if (Instance<T_num>::var(*it).decision_level < DL)
					tmp_clause.push_back(*it);
				else
					lits_at_current_dl++;
				seen[it->var()] = true;
			}
		} else {
			LiteralID alit = Instance<T_num>::getAntecedent(curr_lit).asLit();
			Instance<T_num>::literal(alit).increaseActivity();
			Instance<T_num>::literal(curr_lit).increaseActivity();
			if (!seen[alit.var()] && !(Instance<T_num>::var(alit).decision_level == 0)
					&& !Instance<T_num>::existsUnitClauseOf(alit.var())) {
				if (Instance<T_num>::var(alit).decision_level < DL)
					tmp_clause.push_back(alit);
				else
					lits_at_current_dl++;
				seen[alit.var()] = true;
			}
		}
		curr_lit = NOT_A_LIT;
	}

//	cout << "T" << curr_lit.toInt() << "U "
//     << var(curr_lit).decision_level << ", " << stack_.get_decision_level() << endl;
//	cout << "V"  << var(curr_lit).ante.isAnt() << " "  << endl;
	minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, seen);

//	if (var(curr_lit).decision_level > assertion_level_)
//		assertion_level_ = var(curr_lit).decision_level;
}

template <class T_num>
void Solver<T_num>::recordAllUIPCauses() {
// note:
// variables of lower dl: if seen we dont work with them anymore
// variables of this dl: if seen we incorporate their
// antecedent and set to unseen
	bool seen[Instance<T_num>::num_variables() + 1];
	memset(seen, false, sizeof(bool) * (Instance<T_num>::num_variables() + 1));

	static vector<LiteralID> tmp_clause;
	tmp_clause.clear();

	assertion_level_ = 0;
	uip_clauses_.clear();

	unsigned lit_stack_ofs = literal_stack_.size();
	int DL = stack_.get_decision_level();
	unsigned lits_at_current_dl = 0;

	for (auto l : violated_clause) {
		if (Instance<T_num>::var(l).decision_level == 0 || Instance<T_num>::existsUnitClauseOf(l.var()))
			continue;
		if (Instance<T_num>::var(l).decision_level < DL)
			tmp_clause.push_back(l);
		else
			lits_at_current_dl++;
		Instance<T_num>::literal(l).increaseActivity();
		seen[l.var()] = true;
	}
	unsigned n = 0;
	LiteralID curr_lit;
	while (lits_at_current_dl) {
		assert(lit_stack_ofs != 0);
		curr_lit = literal_stack_[--lit_stack_ofs];

		if (!seen[curr_lit.var()])
			continue;

		seen[curr_lit.var()] = false;

		if (lits_at_current_dl-- == 1) {
			n++;
			if (!Instance<T_num>::hasAntecedent(curr_lit)) {
				// this should be the decision literal when in first branch
				// or it is a literal decided to explore in failed literal testing
				//assert(stack_.TOS_decLit() == curr_lit);
				break;
			}
			// perform UIP stuff
			minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, seen);
		}

		assert(Instance<T_num>::hasAntecedent(curr_lit));

		if (Instance<T_num>::getAntecedent(curr_lit).isAClause()) {
			Instance<T_num>::updateActivities(Instance<T_num>::getAntecedent(curr_lit).asCl());
			assert(curr_lit == *Instance<T_num>::beginOf(Instance<T_num>::getAntecedent(curr_lit).asCl()));

			for (auto it = Instance<T_num>::beginOf(Instance<T_num>::getAntecedent(curr_lit).asCl()) + 1;
					*it != SENTINEL_CL; it++) {
				if (seen[it->var()] || (Instance<T_num>::var(*it).decision_level == 0)
						|| Instance<T_num>::existsUnitClauseOf(it->var()))
					continue;
				if (Instance<T_num>::var(*it).decision_level < DL)
					tmp_clause.push_back(*it);
				else
					lits_at_current_dl++;
				seen[it->var()] = true;
			}
		} else {
			LiteralID alit = Instance<T_num>::getAntecedent(curr_lit).asLit();
			Instance<T_num>::literal(alit).increaseActivity();
			Instance<T_num>::literal(curr_lit).increaseActivity();
			if (!seen[alit.var()] && !(Instance<T_num>::var(alit).decision_level == 0)
					&& !Instance<T_num>::existsUnitClauseOf(alit.var())) {
				if (Instance<T_num>::var(alit).decision_level < DL)
					tmp_clause.push_back(alit);
				else
					lits_at_current_dl++;
				seen[alit.var()] = true;
			}
		}
	}
	if (!Instance<T_num>::hasAntecedent(curr_lit)) {
		minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, seen);
	}
}

template <class T_num>
void Solver<T_num>::printOnlineStats() {
	cout << "c o time elapsed: " << stopwatch_.getElapsedSeconds() << "s" << endl;
	cout << "c o decisions "<<Instance<T_num>::statistics_.num_decisions_<<endl;
	cout << "c o conflicts "<<Instance<T_num>::statistics_.num_conflicts_<<endl;
  cout << "c o conflict clauses (all / bin / unit) \t";
  cout << Instance<T_num>::num_conflict_clauses();
  cout << "/" << Instance<T_num>::statistics_.num_binary_conflict_clauses_ << "/"
      << Instance<T_num>::unit_clauses_.size() << endl;
}

#endif /* SOLVER_H_ */
