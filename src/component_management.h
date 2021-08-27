/*
 * component_management.h
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#ifndef COMPONENT_MANAGEMENT_H_
#define COMPONENT_MANAGEMENT_H_



#include "component_types/component.h"
#include "component_cache.h"
#include "alt_component_analyzer.h"

#include <vector>
#include <gmpxx.h>
#include "containers.h"
#include "stack.h"
#include "hasher.h"

#include "solver_config.h"
using namespace std;

template <class T_num>
class ComponentManager {
public:
  ComponentManager(SolverConfiguration &config, DataAndStatistics<T_num> &statistics,
        LiteralIndexedVector<TriValue> & lit_values,
        LiteralIndexedVector<T_num> & lit_weights) :
        config_(config), statistics_(statistics), cache_(statistics),
        ana_(statistics,lit_values,lit_weights) {
  }

  void initialize(LiteralIndexedVector<Literal> & literals,
        vector<LiteralID> &lit_pool, const Hasher& hasher);

  unsigned scoreOf(VariableIndex v) {
      return ana_.scoreOf(v);
  }

  void cacheModelCountOf(unsigned stack_comp_id, const T_num &value) {
    cache_.storeValueOf(component_stack_[stack_comp_id]->id(), value);
  }

  Component & superComponentOf(StackLevel<T_num> &lev) {
    assert(component_stack_.size() > lev.super_component());
    return *component_stack_[lev.super_component()];
  }

  unsigned component_stack_size() {
    return component_stack_.size();
  }

  void cleanRemainingComponentsOf(StackLevel<T_num> &top) {
    while (component_stack_.size() > top.remaining_components_ofs()) {
      if (cache_.hasEntry(component_stack_.back()->id()))
        cache_.entry(component_stack_.back()->id()).set_deletable();
      delete component_stack_.back();
      component_stack_.pop_back();
    }
    assert(top.remaining_components_ofs() <= component_stack_.size());
  }

  Component & currentRemainingComponentOf(StackLevel<T_num> &top) {
    assert(component_stack_.size() > top.currentRemainingComponent());
    return *component_stack_[top.currentRemainingComponent()];
  }

  // checks for the next yet to explore remaining component of top
  // returns true if a non-trivial non-cached component
  // has been found and is now stack_.TOS_NextComp()
  // returns false if all components have been processed;
  inline bool findNextRemainingComponentOf(StackLevel<T_num> &top, const Hasher& hasher);

  inline void recordRemainingCompsFor(StackLevel<T_num> &top, const Hasher& hasher);

  inline void sortComponentStackRange(unsigned start, unsigned end);

  void gatherStatistics(){
//     statistics_.cache_bytes_memory_usage_ =
//	     cache_.recompute_bytes_memory_usage();
    cache_.compute_byte_size_infrasture();
  }

  void removeAllCachePollutionsOf(StackLevel<T_num> &top);

  AltComponentAnalyzer<T_num> ana_;
private:

  SolverConfiguration &config_;
  DataAndStatistics<T_num> &statistics_;

  vector<Component *> component_stack_;
  ComponentCache<T_num> cache_;
};


template <class T_num>
void ComponentManager<T_num>::sortComponentStackRange(unsigned start, unsigned end){
    assert(start <= end);
    // sort the remaining components for processing
    for (unsigned i = start; i < end; i++)
      for (unsigned j = i + 1; j < end; j++) {
        if (component_stack_[i]->num_variables()
            < component_stack_[j]->num_variables())
          swap(component_stack_[i], component_stack_[j]);
      }
  }

template <class T_num>
bool ComponentManager<T_num>::findNextRemainingComponentOf(StackLevel<T_num> &top, const Hasher& hasher) {
    // record Remaining Components if there are none!
    if (component_stack_.size() <= top.remaining_components_ofs())
      recordRemainingCompsFor(top, hasher);
    assert(!top.branch_found_unsat());
    if (top.hasUnprocessedComponents())
      return true;
    // if no component remains
    // make sure, at least that the current branch is considered SAT
    top.includeSolution(T_num::One());
    return false;
  }

template <class T_num>
inline void ComponentManager<T_num>::recordRemainingCompsFor(StackLevel<T_num> &top, const Hasher& hasher) {
   Component & super_comp = superComponentOf(top);
   unsigned new_comps_start_ofs = component_stack_.size();

   ana_.setupAnalysisContext(top, super_comp);

   for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++)
     if (ana_.isUnseenAndActive(*vt) &&
         ana_.exploreRemainingCompOf(*vt)){

       Component *p_new_comp = ana_.makeComponentFromArcheType();
       CacheableComponent<T_num> *packed_comp = new CacheableComponent<T_num>(ana_.getArchetype().current_comp_for_caching_, hasher);
         if (!cache_.manageNewComponent(top, *packed_comp)){
            component_stack_.push_back(p_new_comp);
            p_new_comp->set_id(cache_.storeAsEntry(*packed_comp, super_comp.id()));
         }
         else {
           delete packed_comp;
           delete p_new_comp;
         }
     }

   top.set_unprocessed_components_end(component_stack_.size());
   sortComponentStackRange(new_comps_start_ofs, component_stack_.size());
}

template<class T_num>
void ComponentManager<T_num>::initialize(LiteralIndexedVector<Literal> & literals,
    vector<LiteralID> &lit_pool, const Hasher& hasher) {

  ana_.initialize(literals, lit_pool);
  // BEGIN CACHE INIT
  //CacheableComponent<T_num>::adjustPackSize(ana_.max_variable_id(), ana_.max_clause_id());

  component_stack_.clear();
  component_stack_.reserve(ana_.max_variable_id() + 2);
  component_stack_.push_back(new Component());
  component_stack_.push_back(new Component());
  assert(component_stack_.size() == 2);
  component_stack_.back()->createAsDummyComponent(ana_.max_variable_id(),
      ana_.max_clause_id());


  cache_.init(*component_stack_.back(), hasher);
}

template<class T_num>
void ComponentManager<T_num>::removeAllCachePollutionsOf(StackLevel<T_num> &top) {
  // all processed components are found in
  // [top.currentRemainingComponent(), component_stack_.size())
  // first, remove the list of descendants from the father
  assert(top.remaining_components_ofs() <= component_stack_.size());
  assert(top.super_component() != 0);
  assert(cache_.hasEntry(superComponentOf(top).id()));

  if (top.remaining_components_ofs() == component_stack_.size())
    return;

  for (unsigned u = top.remaining_components_ofs(); u < component_stack_.size();
      u++) {
    assert(cache_.hasEntry(component_stack_[u]->id()));
    cache_.cleanPollutionsInvolving(component_stack_[u]->id());
  }

#ifdef DEBUG
  cache_.test_descendantstree_consistency();
#endif
}

#endif /* COMPONENT_MANAGEMENT_H_ */
