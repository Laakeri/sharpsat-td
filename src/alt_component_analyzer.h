/*
 * alt_component_analyzer.h
 *
 *  Created on: Mar 5, 2013
 *      Author: mthurley
 */

#ifndef ALT_COMPONENT_ANALYZER_H_
#define ALT_COMPONENT_ANALYZER_H_



#include "statistics.h"
#include "component_types/component.h"
#include "component_types/component_archetype.h"

#include <vector>
#include <cmath>
#include <gmpxx.h>
#include "containers.h"
#include "stack.h"

using namespace std;

template<class T_num>
class AltComponentAnalyzer {
public:
	AltComponentAnalyzer(DataAndStatistics<T_num> &statistics,
        LiteralIndexedVector<TriValue> & lit_values,
        const LiteralIndexedVector<T_num> & lit_weights) :
        statistics_(statistics), literal_values_(lit_values), lit_weights_(lit_weights) {
  }

  unsigned scoreOf(VariableIndex v) {
    return var_frequency_scores_[v];
  }

  ComponentArchetype<T_num> &current_archetype(){
    return archetype_;
  }

  void initialize(LiteralIndexedVector<Literal> & literals,
      vector<LiteralID> &lit_pool);

  int GetSeen(VariableIndex v) {
    assert(v <= max_variable_id_);
    return archetype_.var_seen(v);
  }

  bool isUnseenAndActive(VariableIndex v){
    assert(v <= max_variable_id_);
    return archetype_.var_unseen_in_sup_comp(v);
  }

  // manages the literal whenever it occurs in component analysis
  // returns true iff the underlying variable was unseen before
  //
  bool manageSearchOccurrenceOf(LiteralID lit){
      if(archetype_.var_unseen_in_sup_comp(lit.var())){
        search_stack_.push_back(lit.var());
        archetype_.setVar_seen(lit.var());
        return true;
      }
      return false;
    }
  bool manageSearchOccurrenceAndScoreOf(LiteralID lit){
    var_frequency_scores_[lit.var()]+= isActive(lit);
    return manageSearchOccurrenceOf(lit);
  }

  void setSeenAndStoreInSearchStack(VariableIndex v){
    assert(isActive(v));
    search_stack_.push_back(v);
        archetype_.setVar_seen(v);
  }


  void setupAnalysisContext(StackLevel<T_num> &top, Component & super_comp){
     archetype_.reInitialize(top,super_comp);

     for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++)
       if (isActive(*vt)) {
         archetype_.setVar_in_sup_comp_unseen(*vt);
         var_frequency_scores_[*vt] = 0;
       }

     for (auto itCl = super_comp.clsBegin(); *itCl != clsSENTINEL; itCl++)
         archetype_.setClause_in_sup_comp_unseen(*itCl);
  }

  // returns true, iff the component found is non-trivial
  bool exploreRemainingCompOf(VariableIndex v);


  inline Component *makeComponentFromArcheType(){
    return archetype_.makeComponentFromState(search_stack_.size());
  }

  unsigned max_clause_id(){
     return max_clause_id_;
  }
  unsigned max_variable_id(){
    return max_variable_id_;
  }

  ComponentArchetype<T_num> &getArchetype(){
    return archetype_;
  }

private:
  DataAndStatistics<T_num> &statistics_;

  // the id of the last clause
  // note that clause ID is the clause number,
  // different from the offset of the clause in the literal pool
  unsigned max_clause_id_ = 0;
  unsigned max_variable_id_ = 0;

  // this contains clause offsets of the clauses
  // where each variable occurs in;
  vector<ClauseOfs> variable_occurrence_lists_pool_;

  // this is a new idea,
  // for every variable we have a list
  // 0 binarylinks 0 occurrences 0
  // this should give better cache behaviour,
  // because all links of one variable (binary and nonbinray) are found
  // in one contiguous chunk of memory
  vector<unsigned> unified_variable_links_lists_pool_;

  vector<unsigned> variable_link_list_offsets_;

  LiteralIndexedVector<TriValue> & literal_values_;

  const LiteralIndexedVector<T_num> & lit_weights_;

  inline T_num LitWeight(const LiteralID lit) {
    if (lit_weights_.empty()) {
      return T_num::One();
    } else {
      return lit_weights_[lit];
    }
  }

  vector<unsigned> var_frequency_scores_;

  ComponentArchetype<T_num>  archetype_;

  vector<VariableIndex> search_stack_;

  bool isResolved(const LiteralID lit) {
    return literal_values_[lit] == F_TRI;
  }

  bool isSatisfied(const LiteralID lit) {
    return literal_values_[lit] == T_TRI;
  }
  bool isActive(const LiteralID lit) {
      return literal_values_[lit] == X_TRI;
  }

  bool isActive(const VariableIndex v) {
    return literal_values_[LiteralID(v, true)] == X_TRI;
  }

  unsigned *beginOfLinkList(VariableIndex v) {
    return &unified_variable_links_lists_pool_[variable_link_list_offsets_[v]];
  }

  // stores all information about the component of var
  // in variables_seen_, clauses_seen_ and
  // component_search_stack
  // we have an isolated variable iff
  // after execution component_search_stack.size()==1
  void recordComponentOf(const VariableIndex var);


  void getClause(vector<unsigned> &tmp,
   		       vector<LiteralID>::iterator & it_start_of_cl,
   		       LiteralID & omitLit){
  	  tmp.clear();
  	  for (auto it_lit = it_start_of_cl;*it_lit != SENTINEL_LIT; it_lit++) {
   		  if(it_lit->var() != omitLit.var())
   			 tmp.push_back(it_lit->raw());
   	  }
     }


  void searchClause(VariableIndex vt, ClauseIndex clID, LiteralID * pstart_cls){
    auto itVEnd = search_stack_.end();
    bool all_lits_active = true;
    for (auto itL = pstart_cls; *itL != SENTINEL_LIT; itL++) {
      assert(itL->var() <= max_variable_id_);
      if(!archetype_.var_nil(itL->var()))
        manageSearchOccurrenceAndScoreOf(*itL);
      else {
        assert(!isActive(*itL));
        all_lits_active = false;
        if (isResolved(*itL))
          continue;
        //BEGIN accidentally entered a satisfied clause: undo the search process
        while (search_stack_.end() != itVEnd) {
          assert(search_stack_.back() <= max_variable_id_);
          archetype_.setVar_in_sup_comp_unseen(search_stack_.back());
          search_stack_.pop_back();
        }
        archetype_.setClause_nil(clID);
        while(*itL != SENTINEL_LIT)
          if(isActive(*(--itL)))
            var_frequency_scores_[itL->var()]--;
        //END accidentally entered a satisfied clause: undo the search process
        break;
      }
    }

    if (!archetype_.clause_nil(clID)){
      var_frequency_scores_[vt]++;
      archetype_.setClause_seen(clID,all_lits_active);
    }
  }
};


template<class T_num>
void AltComponentAnalyzer<T_num>::initialize(LiteralIndexedVector<Literal> & literals,
    vector<LiteralID> &lit_pool) {

  max_variable_id_ = literals.end_lit().var() - 1;

  search_stack_.reserve(max_variable_id_ + 1);
  var_frequency_scores_.resize(max_variable_id_ + 1, 0);
  variable_occurrence_lists_pool_.clear();
  variable_link_list_offsets_.clear();
  variable_link_list_offsets_.resize(max_variable_id_ + 1, 0);

  vector<vector<ClauseOfs> > occs(max_variable_id_ + 1);
  vector<vector<unsigned> > occ_long_clauses(max_variable_id_ + 1);
  vector<vector<unsigned> > occ_ternary_clauses(max_variable_id_ + 1);

  vector<unsigned> tmp;
  max_clause_id_ = 0;
  unsigned curr_clause_length = 0;
  auto it_curr_cl_st = lit_pool.begin();

  for (auto it_lit = lit_pool.begin(); it_lit < lit_pool.end(); it_lit++) {
    if (*it_lit == SENTINEL_LIT) {

      if (it_lit + 1 == lit_pool.end())
        break;

      max_clause_id_++;
      it_lit += ClauseHeader::overheadInLits();
      it_curr_cl_st = it_lit + 1;
      curr_clause_length = 0;

    } else {
      assert(it_lit->var() <= max_variable_id_);
      curr_clause_length++;

      getClause(tmp,it_curr_cl_st, *it_lit);

      assert(tmp.size() > 1);

      if(tmp.size() == 2) {
      //if(false){
        occ_ternary_clauses[it_lit->var()].push_back(max_clause_id_);
        occ_ternary_clauses[it_lit->var()].insert(occ_ternary_clauses[it_lit->var()].end(),
            tmp.begin(), tmp.end());
      } else {
        occs[it_lit->var()].push_back(max_clause_id_);
        occs[it_lit->var()].push_back(occ_long_clauses[it_lit->var()].size());
        occ_long_clauses[it_lit->var()].insert(occ_long_clauses[it_lit->var()].end(),
            tmp.begin(), tmp.end());
        occ_long_clauses[it_lit->var()].push_back(SENTINEL_LIT.raw());
      }

    }
  }

  ComponentArchetype<T_num>::initArrays(max_variable_id_, max_clause_id_);
  // the unified link list
  unified_variable_links_lists_pool_.clear();
  unified_variable_links_lists_pool_.push_back(0);
  unified_variable_links_lists_pool_.push_back(0);
  for (unsigned v = 1; v < occs.size(); v++) {
    // BEGIN data for binary clauses
    variable_link_list_offsets_[v] = unified_variable_links_lists_pool_.size();
    for (auto l : literals[LiteralID(v, false)].binary_links_)
      if (l != SENTINEL_LIT)
        unified_variable_links_lists_pool_.push_back(l.var());

    for (auto l : literals[LiteralID(v, true)].binary_links_)
      if (l != SENTINEL_LIT)
        unified_variable_links_lists_pool_.push_back(l.var());

    unified_variable_links_lists_pool_.push_back(0);

    // BEGIN data for ternary clauses
    unified_variable_links_lists_pool_.insert(
        unified_variable_links_lists_pool_.end(),
        occ_ternary_clauses[v].begin(),
        occ_ternary_clauses[v].end());

    unified_variable_links_lists_pool_.push_back(0);

    // BEGIN data for long clauses
    for(auto it = occs[v].begin(); it != occs[v].end(); it+=2){
      unified_variable_links_lists_pool_.push_back(*it);
      unified_variable_links_lists_pool_.push_back(*(it + 1) +(occs[v].end() - it));
    }

    unified_variable_links_lists_pool_.push_back(0);

    unified_variable_links_lists_pool_.insert(
        unified_variable_links_lists_pool_.end(),
        occ_long_clauses[v].begin(),
        occ_long_clauses[v].end());
  }



}

template<typename T_num>
inline bool AltComponentAnalyzer<T_num>::exploreRemainingCompOf(VariableIndex v) {
  assert(archetype_.var_unseen_in_sup_comp(v));
  recordComponentOf(v);

  if (search_stack_.size() == 1) {
    T_num iw = LitWeight(LiteralID(v, true)) + LitWeight(LiteralID(v, false));
    archetype_.stack_level().includeSolution(iw);
    archetype_.setVar_in_other_comp(v);
    return false;
  }
  return true;
}

template <class T_num>
void AltComponentAnalyzer<T_num>::recordComponentOf(const VariableIndex var) {

  search_stack_.clear();
  setSeenAndStoreInSearchStack(var);

  for (auto vt = search_stack_.begin(); vt != search_stack_.end(); vt++) {
    //BEGIN traverse binary clauses
    assert(isActive(*vt));
    unsigned *p = beginOfLinkList(*vt);
    for (; *p; p++) {
      if(manageSearchOccurrenceOf(LiteralID(*p,true))){
        var_frequency_scores_[*p]++;
        var_frequency_scores_[*vt]++;
      }
    }
    //END traverse binary clauses

    for ( p++; *p ; p+=3) {
      if(archetype_.clause_unseen_in_sup_comp(*p)){
        LiteralID litA = *reinterpret_cast<const LiteralID *>(p + 1);
        LiteralID litB = *(reinterpret_cast<const LiteralID *>(p + 1) + 1);
        if(isSatisfied(litA)|| isSatisfied(litB))
          archetype_.setClause_nil(*p);
        else {
          var_frequency_scores_[*vt]++;
          manageSearchOccurrenceAndScoreOf(litA);
          manageSearchOccurrenceAndScoreOf(litB);
          archetype_.setClause_seen(*p,isActive(litA) &
              isActive(litB));
        }
      }
    }
    //END traverse ternary clauses

    for (p++; *p ; p +=2)
      if(archetype_.clause_unseen_in_sup_comp(*p))
        searchClause(*vt,*p, reinterpret_cast<LiteralID *>(p + 1 + *(p+1)));
  }
}



#endif /* ALT_COMPONENT_ANALYZER_H_ */
