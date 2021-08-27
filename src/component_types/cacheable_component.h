/*
 * cacheable_component.h
 *
 *  Created on: Feb 21, 2013
 *      Author: mthurley
 */

#ifndef CACHEABLE_COMPONENT_H_
#define CACHEABLE_COMPONENT_H_

#include <assert.h>
#include <vector>

#include "../primitive_types.h"
#include "../hasher.h"
#include "component.h"

using namespace std;

template<class T_num>
class CacheableComponent {
public:
  CacheableComponent() {
  }

  CacheableComponent(const Component &comp, const Hasher& hasher) {
    num_variables_ = comp.num_variables();
    clhashkey_ = hasher.Hash(comp.RawData());
  }

  unsigned num_variables() const{
    return num_variables_;
  }

  bool equals(const CacheableComponent &comp) const {
    return clhashkey_ == comp.clhashkey();
  }

  HashType clhashkey() const {
    return clhashkey_;
  }

  unsigned hashkey() {
    return clhashkey_[0];
  }

  unsigned creation_time() {
    return creation_time_;
  }
  
  const T_num &model_count() const {
    return model_count_;
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

  // a cache entry is deletable
  // only if it is not connected to an active
  // component in the component stack
  
  bool isDeletable() const {
    return length_solution_period_and_flags_ & 1;
  }

  void set_deletable() {
    length_solution_period_and_flags_ |= 1;
  }

  unsigned long SizeInBytes() const {
    return sizeof(CacheableComponent) + model_count_.InternalSize();
  }

  // TODO make this correct?
  unsigned long sys_overhead_SizeInBytes() const {
    return sizeof(CacheableComponent) + model_count_.InternalSize() + 48;
  }

  // BEGIN Cache Pollution Management

  void set_father(CacheEntryID f) {
    father_ = f;
  }
  const CacheEntryID father() const {
    return father_;
  }

  void set_next_sibling(CacheEntryID sibling) {
    next_sibling_ = sibling;
  }
  CacheEntryID next_sibling() {
    return next_sibling_;
  }

  void set_first_descendant(CacheEntryID descendant) {
    first_descendant_ = descendant;
  }
  CacheEntryID first_descendant() {
    return first_descendant_;
  }

  void set_next_bucket_element(CacheEntryID entry) {
    next_bucket_element_ = entry;
  }

  CacheEntryID next_bucket_element() {
      return next_bucket_element_;
  }

private:
  HashType clhashkey_;
  unsigned num_variables_ = 0;

  T_num model_count_;

  unsigned creation_time_ = 1;
  // this is:  length_solution_period = length_solution_period_and_flags_ >> 1
  // length_solution_period == 0 means unsolved
  // and the first bit is "delete_permitted"
  unsigned length_solution_period_and_flags_ = 0;

  // Cache stuff
  CacheEntryID next_bucket_element_ = 0;
  // theFather and theDescendants:
  // each CCacheEntry is a Node in a tree which represents the relationship
  // of the components stored
  CacheEntryID father_ = 0;
  CacheEntryID first_descendant_ = 0;
  CacheEntryID next_sibling_ = 0;

};


#endif /* CACHEABLE_COMPONENT_H_ */