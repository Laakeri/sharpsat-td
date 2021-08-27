/*
 * instance.h
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#ifndef INSTANCE_H_
#define INSTANCE_H_

#include "statistics.h"
#include "structures.h"
#include "containers.h"
#include "preprocessor/instance.hpp"
#include "preprocessor/graph.hpp"

#include <assert.h>

template <class T_num>
class Instance {
protected:

  void unSet(LiteralID lit) {
    var(lit).ante = Antecedent(NOT_A_CLAUSE);
    var(lit).decision_level = INVALID_DL;
    literal_values_[lit] = X_TRI;
    literal_values_[lit.neg()] = X_TRI;
    //cerr<<"unset literal "<<lit.toInt()<<endl;
  }

  Antecedent & getAntecedent(LiteralID lit) {
    return variables_[lit.var()].ante;
  }

  bool hasAntecedent(LiteralID lit) {
    return variables_[lit.var()].ante.isAnt();
  }

  bool isAntecedentOf(ClauseOfs ante_cl, LiteralID lit) {
    return var(lit).ante.isAClause() && (var(lit).ante.asCl() == ante_cl);
  }

  bool isolated(VariableIndex v) {
    LiteralID lit(v, false);
    return (literal(lit).binary_links_.size() <= 1)
        & occurrence_lists_[lit].empty()
        & (literal(lit.neg()).binary_links_.size() <= 1)
        & occurrence_lists_[lit.neg()].empty();
  }

  bool free(VariableIndex v) {
    return isolated(v) & isActive(v);
  }

  int deleteConflictClauses();
  void deleteAllConflictClauses();
  bool markClauseDeleted(ClauseOfs cl_ofs);

  // Compact the literal pool erasing all the clause
  // information from deleted clauses
  void compactConflictLiteralPool();

  // we assert that the formula is consistent
  // and has not been found UNSAT yet
  // hard wires all assertions in the literal stack into the formula
  // removes all set variables and essentially reinitiallizes all
  // further data
  void compactClauses();
  //void compactVariables();
  void cleanClause(ClauseOfs cl_ofs);

  /////////////////////////////////////////////////////////
  // END access to variables and literals
  /////////////////////////////////////////////////////////


  unsigned int num_conflict_clauses() const {
    return conflict_clauses_.size();
  }

  unsigned int num_variables() {
    return variables_.size() - 1;
  }

  bool createfromFile(const string &file_name);
  bool createfromPPIns(const sspp::Instance& pp_ins);

  DataAndStatistics<T_num> statistics_;

  /** literal_pool_: the literals of all clauses are stored here
   *   INVARIANT: first and last entries of literal_pool_ are a SENTINEL_LIT
   *
   *   Clauses begin with a ClauseHeader structure followed by the literals
   *   terminated by SENTINEL_LIT
   */
  vector<LiteralID> literal_pool_;

  // this is to determine the starting offset of
  // conflict clauses
  unsigned original_lit_pool_size_;

  LiteralIndexedVector<Literal> literals_;

  LiteralIndexedVector<vector<ClauseOfs> > occurrence_lists_;

  LiteralIndexedVector<T_num> lit_weights_;
  vector<char> lit_mul_;
  vector<vector<unsigned long>> dec_cands_;

  vector<ClauseOfs> conflict_clauses_;
  vector<LiteralID> unit_clauses_;

  vector<Variable> variables_;
  LiteralIndexedVector<TriValue> literal_values_;

  vector<unsigned long> gluec_;
  unsigned long gluec_it_ = 1;

  void decayActivities() {
    for (auto l_it = literals_.begin(); l_it != literals_.end(); l_it++)
      l_it->activity_score_ *= 0.5;

  }
//  void decayActivities();

  void updateActivities(ClauseOfs clause_ofs) {
    setGlue(clause_ofs);
    for (auto it = beginOf(clause_ofs); *it != SENTINEL_LIT; it++) {
      literal(*it).increaseActivity();
    }
  }

  void setGlue(ClauseOfs clause_ofs, bool init=false) {
    if (!init && !getHeaderOf(clause_ofs).isLearned()) return;
    gluec_it_++;
    unsigned glue = 0;
    unsigned len = 0;
    for (auto it = beginOf(clause_ofs); *it != SENTINEL_LIT; it++) {
      len++;
      int dl = var(*it).decision_level;
      unsigned var = (*it).var();
      assert(dl != INVALID_DL);
      assert(dl < (int)gluec_.size());
      if (gluec_[dl] != gluec_it_) {
        gluec_[dl] = gluec_it_;
        glue++;
      }
    }
    assert(glue <= len);
    assert(glue >= 1);
    assert(len >= 3);
    getHeaderOf(clause_ofs).setGlue(glue);
  }

  bool isUnitClause(const LiteralID lit) {
    for (auto l : unit_clauses_)
      if (l == lit)
        return true;
    return false;
  }

  bool existsUnitClauseOf(VariableIndex v) {
    for (auto l : unit_clauses_)
      if (l.var() == v)
        return true;
    return false;
  }

  // addUnitClause checks whether lit or lit.neg() is already a
  // unit clause
  // a negative return value implied that the Instance is UNSAT
  bool addUnitClause(const LiteralID lit) {
    for (auto l : unit_clauses_) {
      if (l == lit)
        return true;
      if (l == lit.neg())
        return false;
    }
    unit_clauses_.push_back(lit);
    return true;
  }

  inline ClauseIndex addClause(vector<LiteralID> &literals);

  // adds a UIP Conflict Clause
  // and returns it as an Antecedent to the first
  // literal stored in literals
  inline Antecedent addUIPConflictClause(vector<LiteralID> &literals);

  inline bool addBinaryClause(LiteralID litA, LiteralID litB);

  /////////////////////////////////////////////////////////
  // BEGIN access to variables, literals, clauses
  /////////////////////////////////////////////////////////

  inline Variable &var(const LiteralID lit) {
    return variables_[lit.var()];
  }

  Literal & literal(LiteralID lit) {
    return literals_[lit];
  }

  inline bool isSatisfied(const LiteralID &lit) const {
    return literal_values_[lit] == T_TRI;
  }

  bool isResolved(LiteralID lit) {
    return literal_values_[lit] == F_TRI;
  }

  bool isActive(LiteralID lit) const {
    return literal_values_[lit] == X_TRI;
  }

  vector<LiteralID>::const_iterator beginOf(ClauseOfs cl_ofs) const {
    return literal_pool_.begin() + cl_ofs;
  }
  vector<LiteralID>::iterator beginOf(ClauseOfs cl_ofs) {
    return literal_pool_.begin() + cl_ofs;
  }

  decltype(literal_pool_.begin()) conflict_clauses_begin() {
     return literal_pool_.begin() + original_lit_pool_size_;
   }

  ClauseHeader &getHeaderOf(ClauseOfs cl_ofs) {
    return *reinterpret_cast<ClauseHeader *>(&literal_pool_[cl_ofs
        - ClauseHeader::overheadInLits()]);
  }

  bool isSatisfied(ClauseOfs cl_ofs) {
    for (auto lt = beginOf(cl_ofs); *lt != SENTINEL_LIT; lt++)
      if (isSatisfied(*lt))
        return true;
    return false;
  }

  void LoadWeights(const sspp::Instance& pp_ins, const unsigned int nVars);

protected:
  vector<float> extra_score;
  void PrepareTWScore(const sspp::TreeDecomposition& tdec, double weight, int weight_mode);
};

template <class T_num>
ClauseIndex Instance<T_num>::addClause(vector<LiteralID> &literals) {
  if (literals.size() == 1) {
    //TODO Deal properly with the situation that opposing unit clauses are learned
    assert(!isUnitClause(literals[0].neg()));
    unit_clauses_.push_back(literals[0]);
    return 0;
  }
  if (literals.size() == 2) {
    addBinaryClause(literals[0], literals[1]);
    return 0;
  }
  for (unsigned i = 0; i < ClauseHeader::overheadInLits(); i++)
    literal_pool_.push_back(0);
  ClauseOfs cl_ofs = literal_pool_.size();

  for (auto l : literals) {
    literal_pool_.push_back(l);
    literal(l).increaseActivity(1);
  }
  // make an end: SENTINEL_LIT
  literal_pool_.push_back(SENTINEL_LIT);
  literal(literals[0]).addWatchLinkTo(cl_ofs);
  literal(literals[1]).addWatchLinkTo(cl_ofs);
  getHeaderOf(cl_ofs).set_creation_time(statistics_.num_conflicts_);
  return cl_ofs;
}


template <class T_num>
Antecedent Instance<T_num>::addUIPConflictClause(vector<LiteralID> &literals) {
    Antecedent ante(NOT_A_CLAUSE);
    statistics_.num_clauses_learned_++;
    ClauseOfs cl_ofs = addClause(literals);
    if (cl_ofs != 0) {
      conflict_clauses_.push_back(cl_ofs);
      getHeaderOf(cl_ofs).set_length(literals.size());
      ante = Antecedent(cl_ofs);
      setGlue(cl_ofs, true);
    } else if (literals.size() == 2){
      ante = Antecedent(literals.back());
      statistics_.num_binary_conflict_clauses_++;
    } else if (literals.size() == 1)
      statistics_.num_unit_clauses_++;
    return ante;
  }

template <class T_num>
bool Instance<T_num>::addBinaryClause(LiteralID litA, LiteralID litB) {
   if (literal(litA).hasBinaryLinkTo(litB))
     return false;
   literal(litA).addBinLinkTo(litB);
   literal(litB).addBinLinkTo(litA);
   literal(litA).increaseActivity();
   literal(litB).increaseActivity();
   return true;
 }

 
#include <algorithm>
#include <fstream>
#include <sys/stat.h>

using namespace std;

template <class T_num>
void Instance<T_num>::cleanClause(ClauseOfs cl_ofs) {
  bool satisfied = false;
  for (auto it = beginOf(cl_ofs); *it != SENTINEL_LIT; it++)
    if (isSatisfied(*it)) {
      satisfied = true;
      break;
    }
  // mark the clause as empty if satisfied
  if (satisfied) {
    *beginOf(cl_ofs) = SENTINEL_LIT;
    return;
  }
  auto jt = beginOf(cl_ofs);
  auto it = beginOf(cl_ofs);
  // from now, all inactive literals are resolved
  for (; *it != SENTINEL_LIT; it++, jt++) {
    while (*jt != SENTINEL_LIT && !isActive(*jt))
      jt++;
    *it = *jt;
    if (*jt == SENTINEL_LIT)
      break;
  }
  unsigned length = it - beginOf(cl_ofs);
  // if it has become a unit clause, it should have already been asserted
  if (length == 1) {
    *beginOf(cl_ofs) = SENTINEL_LIT;
    // if it has become binary, transform it to binary and delete it
  } else if (length == 2) {
    addBinaryClause(*beginOf(cl_ofs), *(beginOf(cl_ofs) + 1));
    *beginOf(cl_ofs) = SENTINEL_LIT;
  }
}

template <class T_num>
void Instance<T_num>::compactClauses() {
  vector<ClauseOfs> clause_ofs;
  clause_ofs.reserve(statistics_.num_long_clauses_);

  // clear watch links and occurrence lists
  for (auto it_lit = literal_pool_.begin(); it_lit != literal_pool_.end();
      it_lit++) {
    if (*it_lit == SENTINEL_LIT) {
      if (it_lit + 1 == literal_pool_.end())
        break;
      it_lit += ClauseHeader::overheadInLits();
      clause_ofs.push_back(1 + it_lit - literal_pool_.begin());
    }
  }

  for (auto ofs : clause_ofs)
    cleanClause(ofs);

  for (auto &l : literals_)
    l.resetWatchList();

  occurrence_lists_.clear();
  occurrence_lists_.resize(variables_.size());

  vector<LiteralID> tmp_pool = literal_pool_;
  literal_pool_.clear();
  literal_pool_.push_back(SENTINEL_LIT);
  ClauseOfs new_ofs;
  unsigned num_clauses = 0;
  for (auto ofs : clause_ofs) {
    auto it = (tmp_pool.begin() + ofs);
    if (*it != SENTINEL_LIT) {
      for (unsigned i = 0; i < ClauseHeader::overheadInLits(); i++)
        literal_pool_.push_back(0);
      new_ofs = literal_pool_.size();
      literal(*it).addWatchLinkTo(new_ofs);
      literal(*(it + 1)).addWatchLinkTo(new_ofs);
      num_clauses++;
      for (; *it != SENTINEL_LIT; it++) {
        literal_pool_.push_back(*it);
        occurrence_lists_[*it].push_back(new_ofs);
      }
      literal_pool_.push_back(SENTINEL_LIT);
    }
  }

  vector<LiteralID> tmp_bin;
  unsigned bin_links = 0;
  for (auto &l : literals_) {
    tmp_bin.clear();
    for (auto it = l.binary_links_.begin(); *it != SENTINEL_LIT; it++)
      if (isActive(*it))
        tmp_bin.push_back(*it);
    bin_links += tmp_bin.size();
    tmp_bin.push_back(SENTINEL_LIT);
    l.binary_links_ = tmp_bin;
  }
  statistics_.num_long_clauses_ = num_clauses;
  statistics_.num_binary_clauses_ = bin_links >> 1;
}

/*
template <class T_num>
void Instance<T_num>::compactVariables() {
  vector<unsigned> var_map(variables_.size(), 0);
  unsigned last_ofs = 0;
  unsigned num_isolated = 0;
  LiteralIndexedVector<vector<LiteralID> > _tmp_bin_links(1);
  LiteralIndexedVector<TriValue> _tmp_values = literal_values_;

  for (auto l : literals_)
    _tmp_bin_links.push_back(l.binary_links_);

  assert(_tmp_bin_links.size() == literals_.size());
  for (unsigned v = 1; v < variables_.size(); v++)
    if (isActive(v)) {
      if (isolated(v)) {
        num_isolated++;
        continue;
      }
      last_ofs++;
      var_map[v] = last_ofs;
    }

  variables_.clear();
  variables_.resize(last_ofs + 1);
  occurrence_lists_.clear();
  occurrence_lists_.resize(variables_.size());
  literals_.clear();
  literals_.resize(variables_.size());
  literal_values_.clear();
  literal_values_.resize(variables_.size(), X_TRI);

  unsigned bin_links = 0;
  LiteralID newlit;
  for (auto l = LiteralID(0, false); l != _tmp_bin_links.end_lit(); l.inc()) {
    if (var_map[l.var()] != 0) {
      newlit = LiteralID(var_map[l.var()], l.sign());
      for (auto it = _tmp_bin_links[l].begin(); *it != SENTINEL_LIT; it++) {
        assert(var_map[it->var()] != 0);
        literals_[newlit].addBinLinkTo(
            LiteralID(var_map[it->var()], it->sign()));
      }
      bin_links += literals_[newlit].binary_links_.size() - 1;
    }
  }

  vector<ClauseOfs> clause_ofs;
  clause_ofs.reserve(statistics_.num_long_clauses_);
  // clear watch links and occurrence lists
  for (auto it_lit = literal_pool_.begin(); it_lit != literal_pool_.end();
      it_lit++) {
    if (*it_lit == SENTINEL_LIT) {
      if (it_lit + 1 == literal_pool_.end())
        break;
      it_lit += ClauseHeader::overheadInLits();
      clause_ofs.push_back(1 + it_lit - literal_pool_.begin());
    }
  }

  for (auto ofs : clause_ofs) {
    literal(LiteralID(var_map[beginOf(ofs)->var()], beginOf(ofs)->sign())).addWatchLinkTo(
        ofs);
    literal(LiteralID(var_map[(beginOf(ofs) + 1)->var()],
            (beginOf(ofs) + 1)->sign())).addWatchLinkTo(ofs);
    for (auto it_lit = beginOf(ofs); *it_lit != SENTINEL_LIT; it_lit++) {
      *it_lit = LiteralID(var_map[it_lit->var()], it_lit->sign());
      occurrence_lists_[*it_lit].push_back(ofs);
    }
  }

  literal_values_.clear();
  literal_values_.resize(variables_.size(), X_TRI);
  unit_clauses_.clear();

  statistics_.num_variables_ = variables_.size() - 1 + num_isolated;

  statistics_.num_used_variables_ = num_variables();
  statistics_.num_free_variables_ = num_isolated;
}*/

template <class T_num>
void Instance<T_num>::compactConflictLiteralPool(){
  auto write_pos = conflict_clauses_begin();
  vector<ClauseOfs> tmp_conflict_clauses = conflict_clauses_;
  conflict_clauses_.clear();
  for(auto clause_ofs: tmp_conflict_clauses){
    auto read_pos = beginOf(clause_ofs) - ClauseHeader::overheadInLits();
    for(unsigned i = 0; i < ClauseHeader::overheadInLits(); i++)
      *(write_pos++) = *(read_pos++);
    ClauseOfs new_ofs =  write_pos - literal_pool_.begin();
    conflict_clauses_.push_back(new_ofs);
    // first substitute antecedent if clause_ofs implied something
    if(isAntecedentOf(clause_ofs, *beginOf(clause_ofs)))
      var(*beginOf(clause_ofs)).ante = Antecedent(new_ofs);

    // now redo the watches
    literal(*beginOf(clause_ofs)).replaceWatchLinkTo(clause_ofs,new_ofs);
    literal(*(beginOf(clause_ofs)+1)).replaceWatchLinkTo(clause_ofs,new_ofs);
    // next, copy clause data
    assert(read_pos == beginOf(clause_ofs));
    while(*read_pos != SENTINEL_LIT)
      *(write_pos++) = *(read_pos++);
    *(write_pos++) = SENTINEL_LIT;
  }
  literal_pool_.erase(write_pos,literal_pool_.end());
}

template <class T_num>
void Instance<T_num>::deleteAllConflictClauses() {
  unsigned i2 = 0;
  for(unsigned i = 0; i < conflict_clauses_.size(); i++){
      if(!markClauseDeleted(conflict_clauses_[i])) {
        conflict_clauses_[i2++] = conflict_clauses_[i];
      }
  }
  assert(i2 <= conflict_clauses_.size());
  conflict_clauses_.resize(i2);
}

template <class T_num>
int Instance<T_num>::deleteConflictClauses() {
  statistics_.times_conflict_clauses_cleaned_++;
  vector<ClauseOfs> tmp_conflict_clauses = conflict_clauses_;
  conflict_clauses_.clear();
  vector<unsigned> tmp_glues;
  if (tmp_conflict_clauses.empty()) {
    return true;
  }
  for(auto clause_ofs: tmp_conflict_clauses){
    assert(getHeaderOf(clause_ofs).isLearned());
    tmp_glues.push_back(getHeaderOf(clause_ofs).glue());
    getHeaderOf(clause_ofs).see();
  }
  vector<unsigned> tmp_gluesB = tmp_glues;

  sort(tmp_gluesB.begin(), tmp_gluesB.end());

  unsigned cutoff = tmp_gluesB[tmp_gluesB.size()/2];

  cutoff = max(cutoff, (unsigned)3);

  cout<<"c o cutoff "<<cutoff<<endl;

  for(unsigned i = 0; i < tmp_conflict_clauses.size(); i++){
    if(tmp_glues[i] > cutoff || (tmp_glues[i] == cutoff && conflict_clauses_.size() > tmp_conflict_clauses.size()/2)){
      if(!markClauseDeleted(tmp_conflict_clauses[i])) {
        conflict_clauses_.push_back(tmp_conflict_clauses[i]);
      }
    } else
      conflict_clauses_.push_back(tmp_conflict_clauses[i]);
  }
  return cutoff;
}


template <class T_num>
bool Instance<T_num>::markClauseDeleted(ClauseOfs cl_ofs){
  // only first literal may possibly have cl_ofs as antecedent
  if(isAntecedentOf(cl_ofs, *beginOf(cl_ofs)))
    return false;

  literal(*beginOf(cl_ofs)).removeWatchLinkTo(cl_ofs);
  literal(*(beginOf(cl_ofs)+1)).removeWatchLinkTo(cl_ofs);
  return true;
}

template <class T_num>
bool Instance<T_num>::createfromFile(const string &file_name) {
  unsigned int nVars, nCls;
  int lit;
  unsigned max_ignore = 1000000;
  unsigned clauses_added = 0;
  LiteralID llit;
  vector<LiteralID> literals;
  string idstring;
  char c;

  // clear everything
  literal_pool_.clear();
  literal_pool_.push_back(SENTINEL_LIT);

  variables_.clear();
  variables_.push_back(Variable()); //initializing the Sentinel
  literal_values_.clear();
  unit_clauses_.clear();

  ///BEGIN File input
  ifstream input_file(file_name);
  if (!input_file) {
    cerr << "Cannot open file: " << file_name << endl;
    exit(0);
  }

  struct stat filestatus;
  stat(file_name.c_str(), &filestatus);

  literals.reserve(10000);
  while (input_file >> c && c != 'p')
    input_file.ignore(max_ignore, '\n');
  if (!(input_file >> idstring && idstring == "cnf" && input_file >> nVars
      && input_file >> nCls)) {
    cerr << "Invalid CNF file" << endl;
    exit(0);
  }

  variables_.resize(nVars + 1);
  literal_values_.resize(nVars + 1, X_TRI);
  literal_pool_.reserve(filestatus.st_size);
  conflict_clauses_.reserve(2*nCls);
  occurrence_lists_.clear();
  occurrence_lists_.resize(nVars + 1);

  literals_.clear();
  literals_.resize(nVars + 1);

  while ((input_file >> c) && clauses_added < nCls) {
    input_file.unget(); //extracted a nonspace character to determine if we have a clause, so put it back
    if ((c == '-') || isdigit(c)) {
      literals.clear();
      bool skip_clause = false;
      while ((input_file >> lit) && lit != 0) {
        bool duplicate_literal = false;
        for (auto i : literals) {
          if (i.toInt() == lit) {
            duplicate_literal = true;
            break;
          }
          if (i.toInt() == -lit) {
            skip_clause = true;
            break;
          }
        }
        if (!duplicate_literal) {
          literals.push_back(lit);
        }
      }
      if (!skip_clause) {
        assert(!literals.empty());
        clauses_added++;
        statistics_.incorporateClauseData(literals);
        ClauseOfs cl_ofs = addClause(literals);
        if (literals.size() >= 3)
          for (auto l : literals)
            occurrence_lists_[l].push_back(cl_ofs);
      }
    }
    input_file.ignore(max_ignore, '\n');
  }
  ///END NEW
  input_file.close();
  //  /// END FILE input

  statistics_.num_variables_ = statistics_.num_original_variables_ = nVars;
  statistics_.num_used_variables_ = num_variables();
  statistics_.num_free_variables_ = nVars - num_variables();

  statistics_.num_original_clauses_ = nCls;

  statistics_.num_original_binary_clauses_ = statistics_.num_binary_clauses_;
  statistics_.num_original_unit_clauses_ = statistics_.num_unit_clauses_ =
      unit_clauses_.size();

  original_lit_pool_size_ = literal_pool_.size();
  return true;
}

template <class T_num>
bool Instance<T_num>::createfromPPIns(const sspp::Instance& pp_ins) {
  unsigned int nVars, nCls;
  /*
  int lit;
  unsigned max_ignore = 1000000;
  */
  unsigned clauses_added = 0;
  LiteralID llit;
  vector<LiteralID> literals;
  string idstring;
  //char c;

  // clear everything
  literal_pool_.clear();
  literal_pool_.push_back(SENTINEL_LIT);

  variables_.clear();
  variables_.push_back(Variable()); //initializing the Sentinel
  literal_values_.clear();
  unit_clauses_.clear();


  /*
  ///BEGIN File input
  ifstream input_file(file_name);
  if (!input_file) {
    cerr << "Cannot open file: " << file_name << endl;
    exit(0);
  }

  struct stat filestatus;
  stat(file_name.c_str(), &filestatus);

  literals.reserve(10000);
  while (input_file >> c && c != 'p')
    input_file.ignore(max_ignore, '\n');
  if (!(input_file >> idstring && idstring == "cnf" && input_file >> nVars
      && input_file >> nCls)) {
    cerr << "Invalid CNF file" << endl;
    exit(0);
  }*/

  nVars = pp_ins.vars;
  nCls = pp_ins.clauses.size();

  variables_.resize(nVars + 1);
  literal_values_.resize(nVars + 1, X_TRI);
  literal_pool_.reserve(pp_ins.total_lits + pp_ins.clauses.size() + 10);
  conflict_clauses_.reserve(2*nCls);
  occurrence_lists_.clear();
  occurrence_lists_.resize(nVars + 1);
  lit_weights_.clear();

  literals_.clear();
  literals_.resize(nVars + 1);

  gluec_.resize(nVars*2 + 2);

  for (const auto& clause : pp_ins.clauses) {
    literals.clear();
    for (sspp::Lit lit : clause) {
      literals.push_back(sspp::ToDimacs(lit));
    }
    assert(!literals.empty());
    clauses_added++;
    statistics_.incorporateClauseData(literals);
    ClauseOfs cl_ofs = addClause(literals);
    if (literals.size() >= 3)
      for (auto l : literals)
        occurrence_lists_[l].push_back(cl_ofs);
  }
  if (pp_ins.weighted) {
    LoadWeights(pp_ins, nVars);
  }

  statistics_.num_variables_ = statistics_.num_original_variables_ = nVars;
  statistics_.num_used_variables_ = num_variables();
  statistics_.num_free_variables_ = nVars - num_variables();

  statistics_.num_original_clauses_ = nCls;

  statistics_.num_original_binary_clauses_ = statistics_.num_binary_clauses_;
  statistics_.num_original_unit_clauses_ = statistics_.num_unit_clauses_ =
      unit_clauses_.size();

  original_lit_pool_size_ = literal_pool_.size();
  return true;
}

template<typename T_num>
inline void Instance<T_num>::LoadWeights(const sspp::Instance& pp_ins, const unsigned int nVars) {
  assert(pp_ins.weighted);
  lit_weights_.resize(nVars + 1);
  lit_mul_.resize(nVars+1);
  for (int v = 1; v <= nVars; v++) {
    lit_weights_[LiteralID(-v)].Init(pp_ins.weights[sspp::NegLit(v)]);
    lit_weights_[LiteralID(v)].Init(pp_ins.weights[sspp::PosLit(v)]);
  }
}

template <class T_num>
void Instance<T_num>::PrepareTWScore(const sspp::TreeDecomposition& tdec, double weight, int weight_mode) {
  const int n = literals_.end_lit().var() - 1;
  assert(extra_score.empty());
  extra_score.resize(n+1);
  if (n <= 2) {
    return;
  }
  int width = tdec.Width();
  auto ord = tdec.GetOrd();
  // We use 1-indexing, ignore index 0
  assert(ord.size() == extra_score.size());
  int max_ord = 0;
  for (int i = 1; i <= n; i++) {
    assert(ord[i] >= 1);
    max_ord = max(max_ord, ord[i]);
  }
  assert(max_ord >= 1);
  // Normalize
  for (int i = 1; i <= n; i++) {
    extra_score[i] = max_ord - ord[i];
    extra_score[i] /= (double)max_ord;
    assert(extra_score[i] > -0.01 && extra_score[i] < 1.01);
  }
  // Now scores are between 0..1
  double coef = 1;
  if (weight_mode == 1) {
    if (weight > 0) {
      double rt = (double)n/(double)width;
      if (rt > 40) {
        coef = 1e7;
      } else {
        coef = weight*exp(rt)/(double)n;
      }
    } else {
      coef = 1e7;
    }
  } else if (weight_mode == 2) {
    coef = weight;
  } else if (weight_mode == 3) {
    double rt = (double)n/(double)width;
    if (rt > 40) {
      coef = 1e7;
    } else {
      coef = weight*exp(rt);
    }
  }
  coef = min(coef, 1e7);
  cout << "c o COEF: " << coef << " Width: " << width << endl;
  for (int i = 1; i <= n; i++) {
    extra_score[i] *= coef;
  }
}


#endif /* INSTANCE_H_ */
