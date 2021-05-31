/*
 * base_packed_component.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

#ifndef BASE_PACKED_COMPONENT_H_
#define BASE_PACKED_COMPONENT_H_

#include <assert.h>
#include <gmpxx.h>
#include <iostream>
#include <type_traits>

#include "../hasher.h"
#include "../structures.h"

using namespace std;

template <class T_num>
class BasePackedComponent {
public:
  BasePackedComponent() {}
  BasePackedComponent(unsigned creation_time): creation_time_(creation_time) {}

  unsigned creation_time() {
    return creation_time_;
  }

  const T_num &model_count() const {
    return model_count_;
  }

  unsigned alloc_of_model_count() const{
    if (std::is_same<T_num, mpz_class>::value) {
      return sizeof(T_num) + model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
    } else if (std::is_same<T_num, double>::value) {
      return sizeof(T_num);
    } else if (std::is_same<T_num, LogNum>::value) {
      return sizeof(T_num);
    } else {
      static_assert(std::is_same<T_num, mpz_class>::value || std::is_same<T_num, double>::value, "Type must be mpz_class or double");
      assert(false);
    }
  }

  void set_creation_time(unsigned time) {
    creation_time_ = time;
  }

  void set_model_count(const T_num &rn, unsigned time) {
    model_count_ = rn;
    length_solution_period_and_flags_ = (time - creation_time_) | (length_solution_period_and_flags_ & 1);
  }

  bool modelCountFound(){
    return (length_solution_period_and_flags_ >> 1);
  }

 // inline bool equals(const BasePackedComponent &comp) const;

  // a cache entry is deletable
  // only if it is not connected to an active
  // component in the component stack
  bool isDeletable() const {
    return length_solution_period_and_flags_ & 1;
  }
  void set_deletable() {
    length_solution_period_and_flags_ |= 1;
  }

  void clear() {
    // before deleting the contents of this component,
    // we should make sure that this component is not present in the component stack anymore!
    assert(isDeletable());
  }

  static unsigned _debug_static_val;

protected:
  // data_ contains in packed form the variable indices
  // and clause indices of the component ordered
  // structure is
  // var var ... clause clause ...
  // clauses begin at clauses_ofs_

  //unsigned hashkey_ = 0;

  T_num model_count_;

  unsigned creation_time_ = 1;


  // this is:  length_solution_period = length_solution_period_and_flags_ >> 1
  // length_solution_period == 0 means unsolved
  // and the first bit is "delete_permitted"
  unsigned length_solution_period_and_flags_ = 0;

  // deletion is permitted only after
  // the copy of this component in the stack
  // does not exist anymore

};

#endif /* BASE_PACKED_COMPONENT_H_ */
