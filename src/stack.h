/*
 * stack.h
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#ifndef STACK_H_
#define STACK_H_

#include <gmpxx.h>
#include <vector>

template<class T_num>
class StackLevel {
public:
  /// active Component, once initialized, it should not change
  const unsigned super_component_ = 0;
  // branch
  bool active_branch_ = false;

  // offset in the literal stack where to store set lits
  const unsigned literal_stack_ofs_ = 0;

  //  Solutioncount
  T_num branch_model_count_[2] = {0,0};
  bool branch_found_unsat_[2] = {false,false};

  /// remaining Components

  // the start offset in the component stack for
  // the remaining components in this decision level
  // all remaining components can hence be found in
  // [remaining_components_ofs_, "nextLevel".remaining_components_begin_)
  const unsigned remaining_components_ofs_ = 0;

  // boundary of the stack marking which components still need to be processed
  // all components to be processed can be found in
  // [remaining_components_ofs_, unprocessed_components_end_)
  // also, all processed, can be found
  // in [unprocessed_components_end_, component_stack.size())
  unsigned unprocessed_components_end_ = 0;

  double dec_weight = 1;
  LogNum dec_weight_ln = (int)1;

  bool hasUnprocessedComponents() {
    assert(unprocessed_components_end_ >= remaining_components_ofs_);
    return unprocessed_components_end_ > remaining_components_ofs_;
  }
  void nextUnprocessedComponent() {
    assert(unprocessed_components_end_ > remaining_components_ofs_);
    unprocessed_components_end_--;
  }

  void resetRemainingComps() {
    unprocessed_components_end_ = remaining_components_ofs_;
  }
  unsigned super_component() {
    return super_component_;
  }
  unsigned remaining_components_ofs() {
    return remaining_components_ofs_;
  }
  void set_unprocessed_components_end(unsigned end) {
    unprocessed_components_end_ = end;
    assert(remaining_components_ofs_ <= unprocessed_components_end_);
  }

  StackLevel(unsigned super_comp, unsigned lit_stack_ofs,
      unsigned comp_stack_ofs) :
      super_component_(super_comp),
      literal_stack_ofs_(lit_stack_ofs),
      remaining_components_ofs_(comp_stack_ofs),
      unprocessed_components_end_(comp_stack_ofs) {
    assert(super_comp < comp_stack_ofs);
  }

  unsigned currentRemainingComponent() {
    assert(remaining_components_ofs_ <= unprocessed_components_end_ - 1);
    return unprocessed_components_end_ - 1;
  }
  bool isSecondBranch() {
    return active_branch_;
  }

  void changeBranch() {
    active_branch_ = true;
  }

  bool anotherCompProcessible() {
    return (!branch_found_unsat()) && hasUnprocessedComponents();
  }

  unsigned literal_stack_ofs() {
    return literal_stack_ofs_;
  }

  void includeSolution(const T_num &solutions);
  void includeSolution(unsigned solutions);

  bool branch_found_unsat() {
    return branch_found_unsat_[active_branch_];
  }
  void mark_branch_unsat() {
    branch_found_unsat_[active_branch_] = true;
  }

//  void set_both_branches_unsat(){
//	  branch_found_unsat_[0] =
//			  branch_found_unsat_[1] = true;
//	  branch_model_count_[0] = branch_model_count_[1] = 0;
//	  active_branch_ = 1;
//  }
  const T_num getTotalModelCount() const;
};

template <>
inline void StackLevel<LogNum>::includeSolution(unsigned solutions) {
  //cerr<<"includesolution_ln "<<solutions<<" "<<dec_weight.get()<<endl;
  if (branch_found_unsat_[active_branch_]) {
    assert(branch_model_count_[active_branch_] == 0);
    return;
  }
  if (solutions == 0)
    branch_found_unsat_[active_branch_] = true;
  if (branch_model_count_[active_branch_] == 0) {
    branch_model_count_[active_branch_] = dec_weight_ln * (int)solutions;
  }
  else {
    branch_model_count_[active_branch_] *= solutions;
  }
}
template <>
inline void StackLevel<LogNum>::includeSolution(const LogNum &solutions) {
  //cerr<<"includesolution_d "<<solutions<<" "<<dec_weight<<endl;
  if (branch_found_unsat_[active_branch_]) {
    assert(branch_model_count_[active_branch_] == 0);
    return;
  }
  if (solutions == 0)
    branch_found_unsat_[active_branch_] = true;
  if (branch_model_count_[active_branch_] == 0) {
    branch_model_count_[active_branch_] = solutions * dec_weight_ln;
  }
  else {
    branch_model_count_[active_branch_] *= solutions;
  }
}
template <>
inline void StackLevel<double>::includeSolution(const double &solutions) {
  //cerr<<"includesolution_d "<<solutions<<" "<<dec_weight<<endl;
  if (branch_found_unsat_[active_branch_]) {
    assert(branch_model_count_[active_branch_] == 0);
    return;
  }
  if (solutions == 0)
    branch_found_unsat_[active_branch_] = true;
  if (branch_model_count_[active_branch_] == 0) {
    branch_model_count_[active_branch_] = solutions * dec_weight;
  }
  else {
    branch_model_count_[active_branch_] *= solutions;
  }
}
template <>
inline void StackLevel<double>::includeSolution(unsigned solutions) {
  //cerr<<"includesolution_u "<<solutions<<" "<<dec_weight<<endl;
  if (branch_found_unsat_[active_branch_]) {
    assert(branch_model_count_[active_branch_] == 0);
    return;
  }
  if (solutions == 0)
    branch_found_unsat_[active_branch_] = true;
  if (branch_model_count_[active_branch_] == 0) {
    branch_model_count_[active_branch_] = (double)solutions * dec_weight;
  }
  else {
    branch_model_count_[active_branch_] *= (double)solutions;
  }
}

template <>
inline void StackLevel<mpz_class>::includeSolution(const mpz_class &solutions) {
  if (branch_found_unsat_[active_branch_]) {
    assert(branch_model_count_[active_branch_] == 0);
    return;
  }
  if (solutions == 0)
    branch_found_unsat_[active_branch_] = true;
  if (branch_model_count_[active_branch_] == 0)
    branch_model_count_[active_branch_] = solutions;
  else
    branch_model_count_[active_branch_] *= solutions;
}
template <>
inline void StackLevel<mpz_class>::includeSolution(unsigned solutions) {
  if (branch_found_unsat_[active_branch_]) {
    assert(branch_model_count_[active_branch_] == 0);
    return;
  }
  if (solutions == 0)
    branch_found_unsat_[active_branch_] = true;
  if (branch_model_count_[active_branch_] == 0)
    branch_model_count_[active_branch_] = solutions;
  else
    branch_model_count_[active_branch_] *= solutions;
}

template <>
const inline LogNum StackLevel<LogNum>::getTotalModelCount() const {
  return branch_model_count_[0] + branch_model_count_[1];
}

template <>
const inline double StackLevel<double>::getTotalModelCount() const {
  return branch_model_count_[0] + branch_model_count_[1];
}

template <>
const inline mpz_class StackLevel<mpz_class>::getTotalModelCount() const {
  return branch_model_count_[0] + branch_model_count_[1];
}

template<class T_num>
class DecisionStack: public vector<StackLevel<T_num>> {
public:

  StackLevel<T_num> &top() {
    assert(vector<StackLevel<T_num>>::size() > 0);
    return vector<StackLevel<T_num>>::back();
  }
  int get_decision_level() const {
    assert(vector<StackLevel<T_num>>::size() > 0);
    return vector<StackLevel<T_num>>::size() - 1;
  } // 0 means pre-1st-decision

};



#endif /* STACK_H_ */
