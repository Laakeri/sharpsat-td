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

#include "difference_packed_component.h"
//#include "simple_unpacked_component.h"


using namespace std;

#define NIL_ENTRY 0

class Component;

template<class T_num>
class ComponentArchetype;


// GenericCacheableComponent Adds Structure to PackedComponent that is
// necessary to store it in the cache
// namely, the descendant tree structure that
// allows for the removal of cache pollutions

template< class T_Component>
class GenericCacheableComponent: public T_Component {
public:
  GenericCacheableComponent() {
  }

  GenericCacheableComponent(Component &comp, const Hasher& hasher) :
      T_Component(comp, hasher) {
  }

  unsigned long SizeInBytes() const {
    return sizeof(GenericCacheableComponent<T_Component>)
        + T_Component::raw_data_byte_size();
  }

  // the 48 = 16*3 in overhead stems from the three parts of the component
  // being dynamically allocated (i.e. the GenericCacheableComponent itself,
  // the data_ and the model_count data
  unsigned long sys_overhead_SizeInBytes() const {
      return sizeof(GenericCacheableComponent<T_Component>)
          + T_Component::sys_overhead_raw_data_byte_size()
         // + 24;
          +32;
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

  CacheEntryID next_bucket_element_ = 0;

  // theFather and theDescendants:
  // each CCacheEntry is a Node in a tree which represents the relationship
  // of the components stored
  CacheEntryID father_ = 0;
  CacheEntryID first_descendant_ = 0;
  CacheEntryID next_sibling_ = 0;

};


template <class T_num>
using CacheableComponent = GenericCacheableComponent<DifferencePackedComponent<T_num>>;
//typedef GenericCacheableComponent<DifferencePackedComponent> CacheableComponent;
//typedef GenericCacheableComponent<SimplePackedComponent> CacheableComponent;





#endif /* CACHEABLE_COMPONENT_H_ */
