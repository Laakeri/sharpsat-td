/*
 * statistics.h
 *
 *  Created on: Feb 13, 2013
 *      Author: mthurley
 */

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <string>
#include <cstdint>
#include <vector>

#include <gmpxx.h>

#include "structures.h"
#include "component_types/cacheable_component.h"

#include "primitive_types.h"

using namespace std;

template <class T_num>
class DataAndStatistics {
public:
  string input_file_;
  double time_elapsed_ = 0.0;
  uint64_t maximum_cache_size_bytes_ = 0;

  SOLVER_StateT exit_state_ = NO_STATE;
  // different variable counts
  // number of variables  and clauses before preprocessing
  unsigned long num_original_variables_ = 0;
  unsigned long num_original_clauses_ = 0;
  unsigned long num_original_binary_clauses_ = 0;
  unsigned long num_original_unit_clauses_ = 0;

  // number of variables remaining
  unsigned long num_variables_ = 0;
  // number of variables that actually occurs in clauses
  unsigned long num_used_variables_ = 0;
  unsigned long num_free_variables_ = 0;

  /// different clause counts

  // number of clauses after preprocessing
  unsigned long num_long_clauses_ = 0;
  unsigned long num_binary_clauses_ = 0;

  unsigned long num_long_conflict_clauses_ = 0;
  unsigned long num_binary_conflict_clauses_ = 0;

  unsigned long times_conflict_clauses_cleaned_ = 0;

  unsigned long num_unit_clauses_ = 0;

  /// number of all decisions made
  unsigned long num_decisions_ = 0;
  /// number of all implications derived
  unsigned long num_implications_ = 0;
  // number of all failed literal detections
  unsigned long num_failed_literals_detected_ = 0;
  unsigned long num_failed_literal_tests_ = 0;
  // number of all conflicts occurred
  unsigned long num_conflicts_ = 0;

  // number of clauses overall learned
  unsigned num_clauses_learned_ = 0;


  /* cache statistics */
  uint64_t num_cache_hits_ = 0;
  uint64_t num_cache_look_ups_ = 0;
  uint64_t sum_cache_hit_sizes_ = 0;

  uint64_t num_cached_components_ = 0;
  uint64_t sum_size_cached_components_ = 0;

  // the number of bytes occupied by all
  // components
  uint64_t sum_bytes_cached_components_ = 0;
  // the same number, summing over all components ever stored
  uint64_t overall_bytes_components_stored_ = 0;

  // the above numbers, but without any overhead,
  // counting only the pure data size of the components - without model counts
  uint64_t sum_bytes_pure_cached_component_data_ = 0;
  // the same number, summing over all components ever stored
  uint64_t overall_bytes_pure_stored_component_data_ = 0;


  uint64_t sys_overhead_sum_bytes_cached_components_ = 0;
    // the same number, summing over all components ever stored
  uint64_t sys_overhead_overall_bytes_components_stored_ = 0;

  uint64_t cache_infrastructure_bytes_memory_usage_ = 0;


  uint64_t overall_num_cache_stores_ = 0;
  /*end statistics */

  bool cache_full(){
    return cache_bytes_memory_usage() >= maximum_cache_size_bytes_;
  }

  uint64_t cache_bytes_memory_usage(){
    return cache_infrastructure_bytes_memory_usage_
           + sum_bytes_cached_components_;
  }

  uint64_t overall_cache_bytes_memory_stored(){
      return cache_infrastructure_bytes_memory_usage_
             + overall_bytes_components_stored_;
    }

  void incorporate_cache_store(CacheableComponent<T_num> &ccomp){
    sum_bytes_cached_components_ += ccomp.SizeInBytes();
    sum_size_cached_components_ += ccomp.num_variables();
    num_cached_components_++;
    overall_bytes_components_stored_ += ccomp.SizeInBytes();
    overall_num_cache_stores_ += ccomp.num_variables();
    sys_overhead_sum_bytes_cached_components_ += ccomp.sys_overhead_SizeInBytes();
    sys_overhead_overall_bytes_components_stored_ += ccomp.sys_overhead_SizeInBytes();

  }
  void incorporate_cache_erase(CacheableComponent<T_num> &ccomp){
      sum_bytes_cached_components_ -= ccomp.SizeInBytes();
      sum_size_cached_components_ -= ccomp.num_variables();
      num_cached_components_--;

      sys_overhead_sum_bytes_cached_components_ -= ccomp.sys_overhead_SizeInBytes();
  }

  void incorporate_cache_hit(CacheableComponent<T_num> &ccomp){
      num_cache_hits_++;
      sum_cache_hit_sizes_ += ccomp.num_variables();
  }
  unsigned long cache_MB_memory_usage() {
      return cache_bytes_memory_usage() / 1000000;
  }
  T_num final_solution_count_ = T_num::Zero();

  double implicitBCP_miss_rate() {
      if(num_failed_literal_tests_ == 0) return 0.0;
      return (num_failed_literal_tests_ - num_failed_literals_detected_) / (double) num_failed_literal_tests_;
  }
  unsigned long num_clauses() {
    return num_long_clauses_ + num_binary_clauses_ + num_unit_clauses_;
  }
  unsigned long num_conflict_clauses() {
    return num_long_conflict_clauses_ + num_binary_conflict_clauses_;
  }

  unsigned long clause_deletion_interval() {
    return 10000 + 10 * times_conflict_clauses_cleaned_;
  }

  void set_final_solution_count(const T_num &count);

  const T_num &final_solution_count() const {
    return final_solution_count_;
  }

  void incorporateConflictClauseData(const vector<LiteralID> &clause) {
    if (clause.size() == 1)
      num_unit_clauses_++;
    else if (clause.size() == 2)
      num_binary_conflict_clauses_++;
    num_long_conflict_clauses_++;
  }
  void incorporateClauseData(const vector<LiteralID> &clause) {
    if (clause.size() == 1)
      num_unit_clauses_++;
    else if (clause.size() == 2)
      num_binary_clauses_++;
    else
      num_long_clauses_++;
  }

  void printShort();

  void printShortFormulaInfo() {
    cout << "c o variables (all/used/free): \t";
    cout << num_variables_ << "/" << num_used_variables_ << "/";
    cout << num_variables_ - num_used_variables_ << endl;

    cout << "c o clauses (all/long/binary/unit): ";
    cout << num_clauses() << "/" << num_long_clauses_;
    cout << "/" << num_binary_clauses_ << "/" << num_unit_clauses_ << endl;
  }
  unsigned getTime() {
    return num_decisions_;
  }

  double avgCachedSize() {
    if (num_cache_hits_ == 0)
      return 0.0;
    return (double) sum_size_cached_components_
        / (double) num_cached_components_;
  }

  double avgCacheHitSize() {
    if (num_cache_hits_ == 0)
      return 0.0;
    return (double) sum_cache_hit_sizes_ / (double) num_cache_hits_;
  }

  long double getAvgComponentSize() {
    return sum_size_cached_components_ / (long double) num_cached_components_;
  }

  unsigned long cached_component_count() {
    return num_cached_components_;
  }

  unsigned long cache_hits() {
    return num_cache_hits_;
  }

  double cache_miss_rate() {
    if(num_cache_look_ups_ == 0) return 0.0;
    return (num_cache_look_ups_ - num_cache_hits_)
        / (double) num_cache_look_ups_;
  }

  long double getAvgCacheHitSize() {
    if(num_cache_hits_ == 0) return 0.0;
    return sum_cache_hit_sizes_ / (long double) num_cache_hits_;
  }

};

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

template<typename T_num>
inline void DataAndStatistics<T_num>::set_final_solution_count(const T_num &count) {
  assert(num_variables_ == num_used_variables_);
  final_solution_count_ = count;
}

template <class T_num>
void DataAndStatistics<T_num>::printShort() {
  /*
  if (exit_state_ == TIMEOUT) {
    cout << endl << " TIMEOUT !" << endl;
    return;
  }
  cout << endl << endl;
  cout << "variables (total / active / free)\t" << num_variables_ << "/"
      << num_used_variables_ << "/" << num_variables_ - num_used_variables_
      << endl;
  cout << "clauses (removed) \t\t\t" << num_original_clauses_ << " ("
      << num_original_clauses_ - num_clauses() << ")" << endl;
  */
  cout << "c o decisions \t\t\t\t" << num_decisions_ << endl;
  cout << "c o conflicts \t\t\t\t" << num_conflicts_ << endl;
  cout << "c o conflict clauses (all/bin/unit) \t";
  cout << num_conflict_clauses();
  cout << "/" << num_binary_conflict_clauses_ << "/" << num_unit_clauses_
      << endl;
  cout << "c o bytes cache size     \t" << cache_bytes_memory_usage()  << "\t"
      << endl;

  cout << "c o bytes cache (overall) \t" << overall_cache_bytes_memory_stored()
      << "" << endl;
  cout << "c o bytes cache (infra / comps) "
      << (cache_infrastructure_bytes_memory_usage_) << "/"
      << sum_bytes_cached_components_  << "\t" << endl;

  cout << "c o bytes pure comp data (curr)    " << sum_bytes_pure_cached_component_data_  << "" << endl;
  cout << "c o bytes pure comp data (overall) " <<overall_bytes_pure_stored_component_data_ << "" << endl;

  cout << "c o bytes cache with sysoverh (curr)    " << sys_overhead_sum_bytes_cached_components_  << "" << endl;
  cout << "c o bytes cache with sysoverh (overall) " << sys_overhead_overall_bytes_components_stored_ << "" << endl;


  cout << "c o cache (stores / hits) \t\t\t" << num_cached_components_ << "/"
      << num_cache_hits_ << endl;
  cout << "c o cache miss rate " << cache_miss_rate() * 100 << "%" << endl;
  cout << "c o avg. variable count (stores / hits) \t" << getAvgComponentSize()
      << "/" << getAvgCacheHitSize() << endl;
  /*
  cout << "\n# solutions " << endl;

  print_final_solution_count();

  cout << "\n# END" << endl << endl;
  cout << "time: " << time_elapsed_ << "s\n\n";
  */
}


#endif /* STATISTICS_H_ */
