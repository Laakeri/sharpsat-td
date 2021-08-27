/*
 * component_cache.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

#ifndef COMPONENT_CACHE_H_
#define COMPONENT_CACHE_H_


#include "component_types/cacheable_component.h"
#include "statistics.h"

#include <gmpxx.h>

#include "component_types/component.h"

#include "stack.h"
#include "hasher.h"

template<class T_num>
class ComponentCache {
public:

  ComponentCache(DataAndStatistics<T_num> &statistics);

  ~ComponentCache() {
   // debug_dump_data();
    for (auto &pentry : entry_base_)
          if (pentry != nullptr)
            delete pentry;
  }

  void init(Component &super_comp, const Hasher& hasher);

  // compute the size in bytes of the component cache from scratch
  // the value is stored in bytes_memory_usage_
  uint64_t compute_byte_size_infrasture();

  CacheableComponent<T_num> &entry(CacheEntryID id) {
    assert(entry_base_.size() > id);
    assert(entry_base_[id] != nullptr);
    return *entry_base_[id];
  }

  CacheableComponent<T_num> &entry(const Component& comp) {
      return entry(comp.id());
  }

  bool hasEntry(CacheEntryID id) {
    assert(entry_base_.size() > id);
    return entry_base_[id];
  }

  // removes the entry id from the hash table
  // but not from the entry base
  inline void removeFromHashTable(CacheEntryID id);

  // we delete the Component with ID id
  // and all its descendants from the cache
  inline void cleanPollutionsInvolving(CacheEntryID id);

  // creates a CCacheEntry in the entry base
  // which contains a packed copy of comp
  // returns the id of the entry created
  // stores in the entry the position of
  // comp which is a part of the component stack
  inline CacheEntryID storeAsEntry(CacheableComponent<T_num> &ccomp,
                            CacheEntryID super_comp_id);

  // check quickly if the model count of the component is cached
  // if so, incorporate it into the model count of top
  // if not, store the packed version of it in the entry_base of the cache
  bool manageNewComponent(StackLevel<T_num> &top, CacheableComponent<T_num> &packed_comp) {
   statistics_.num_cache_look_ups_++;
   unsigned table_ofs =  packed_comp.hashkey() & table_size_mask_;

   CacheEntryID act_id = table_[table_ofs];
   if (!act_id) {
    return false;
   }
   while(act_id){
     if (entry(act_id).equals(packed_comp)) {
       statistics_.incorporate_cache_hit(packed_comp);
       top.includeSolution(entry(act_id).model_count());
       return true;
     }
     act_id = entry(act_id).next_bucket_element();
   }
   return false;
  }


  // unchecked erase of an entry from entry_base_
  void eraseEntry(CacheEntryID id) {
    statistics_.incorporate_cache_erase(*entry_base_[id]);
    delete entry_base_[id];
    entry_base_[id] = nullptr;
    free_entry_base_slots_.push_back(id);
  }


  // store the number in model_count as the model count of CacheEntryID id
  inline void storeValueOf(CacheEntryID id, const T_num &model_count);

  bool deleteEntries();

  // delete entries, keeping the descendants tree consistent
  inline void removeFromDescendantsTree(CacheEntryID id);

  // test function to ensure consistency of the descendant tree
  inline void test_descendantstree_consistency();
  void debug_dump_data();
private:

  void considerCacheResize(){
    if (entry_base_.size() > table_.size()) {
      reHashTable(2*table_.size());
    }
  }
  void reHashTable(unsigned size){

    table_.clear();
    table_.resize(size,0);
    // we assert that table size is a power of 2
    // otherwise the table_size_mask_ doesn't work
    assert((table_.size() & (table_.size() - 1)) == 0);
    table_size_mask_ = table_.size() - 1;
    cout << "c o ts " << table_.size() << " " << table_size_mask_ << endl;
    unsigned collisions = 0;
    for (unsigned id = 2; id < entry_base_.size(); id++)
      if (entry_base_[id] != nullptr ){
        entry_base_[id]->set_next_bucket_element(0);
       if(entry_base_[id]->modelCountFound()) {
        unsigned table_ofs=tableEntry(id);
        collisions += (table_[table_ofs] > 0 ? 1 : 0);
        entry_base_[id]->set_next_bucket_element(table_[table_ofs]);
        table_[table_ofs] = id;
       }
    }
    cout << "c o coll " << collisions << endl;
  }

  unsigned tableEntry(CacheEntryID id){
    return entry(id).hashkey() & table_size_mask_;
  }
  void add_descendant(CacheEntryID compid, CacheEntryID descendantid) {
      assert(descendantid != entry(compid).first_descendant());
      entry(descendantid).set_next_sibling(entry(compid).first_descendant());
      entry(compid).set_first_descendant(descendantid);
    }

  void remove_firstdescendantOf(CacheEntryID compid) {
      CacheEntryID desc = entry(compid).first_descendant();
      if (desc != 0)
        entry(compid).set_first_descendant(entry(desc).next_sibling());
    }

  vector<CacheableComponent<T_num> *> entry_base_;
  vector<CacheEntryID> free_entry_base_slots_;

  // the actual hash table
  // by means of which the cache is accessed
  vector<CacheEntryID> table_;

  unsigned table_size_mask_;

  DataAndStatistics<T_num> &statistics_;

  unsigned long my_time_ = 0;
};


template<class T_num>
CacheEntryID ComponentCache<T_num>::storeAsEntry(CacheableComponent<T_num> &ccomp, CacheEntryID super_comp_id){
    CacheEntryID id;

    if (statistics_.cache_full())
        deleteEntries();

    assert(!statistics_.cache_full());

    ccomp.set_creation_time(my_time_++);

    if (free_entry_base_slots_.empty()) {
        if (entry_base_.capacity() == entry_base_.size()) {
            entry_base_.reserve(2 * entry_base_.size());
        }
        entry_base_.push_back(&ccomp);
        id = entry_base_.size() - 1;
    } else {
        id = free_entry_base_slots_.back();
        assert(id < entry_base_.size());
        assert(entry_base_[id] == nullptr);
        free_entry_base_slots_.pop_back();
        entry_base_[id] = &ccomp;
    }

    entry(id).set_father(super_comp_id);
    add_descendant(super_comp_id, id);

    assert(hasEntry(id));
    assert(hasEntry(super_comp_id));

    statistics_.incorporate_cache_store(ccomp);

  #ifdef DEBUG
      for (unsigned u = 2; u < entry_base_.size(); u++)
            if (entry_base_[u] != nullptr) {
              assert(entry_base_[u]->father() != id);
              assert(entry_base_[u]->first_descendant() != id);
              assert(entry_base_[u]->next_sibling() != id);
            }
  #endif
    return id;
}

template <class T_num>
void ComponentCache<T_num>::cleanPollutionsInvolving(CacheEntryID id) {
  CacheEntryID father = entry(id).father();
  if (entry(father).first_descendant() == id) {
    entry(father).set_first_descendant(entry(id).next_sibling());
  } else {
    CacheEntryID act_sibl = entry(father).first_descendant();
    while (act_sibl) {
      CacheEntryID next_sibl = entry(act_sibl).next_sibling();
      if (next_sibl == id) {
        entry(act_sibl).set_next_sibling(entry(next_sibl).next_sibling());
        break;
      }
      act_sibl = next_sibl;
    }
  }
  CacheEntryID next_child = entry(id).first_descendant();
  entry(id).set_first_descendant(0);
  while (next_child) {
    CacheEntryID act_child = next_child;
    next_child = entry(act_child).next_sibling();
    cleanPollutionsInvolving(act_child);
  }
  removeFromHashTable(id);
  eraseEntry(id);
}

template <class T_num>
void ComponentCache<T_num>::removeFromHashTable(CacheEntryID id) {
  //assert(false);
  unsigned act_id = table_[tableEntry(id)];
  if(act_id == id){
    table_[tableEntry(id)] = entry(act_id).next_bucket_element();
  }
  else {
  while (act_id) {
        CacheEntryID next_id = entry(act_id).next_bucket_element();
        if (next_id == id) {
          entry(act_id).set_next_bucket_element(entry(next_id).next_bucket_element());
          break;
        }
        act_id = next_id;
      }
  }
//  CacheBucket *p_bucket = bucketOf(entry(id));
//  if(p_bucket)
//    for (auto it = p_bucket->begin(); it != p_bucket->end(); it++)
//      if (*it == id) {
//        *it = p_bucket->back();
//        p_bucket->pop_back();
//        break;
//      }
}

template <class T_num>
void ComponentCache<T_num>::removeFromDescendantsTree(CacheEntryID id) {
  assert(hasEntry(id));
  // we need a father for this all to work
  assert(entry(id).father());
  assert(hasEntry(entry(id).father()));
  // two steps
  // 1. remove id from the siblings list
  CacheEntryID father = entry(id).father();
  if (entry(father).first_descendant() == id) {
    entry(father).set_first_descendant(entry(id).next_sibling());
  } else {
    CacheEntryID act_sibl = entry(father).first_descendant();
    while (act_sibl) {
      CacheEntryID next_sibl = entry(act_sibl).next_sibling();
      if (next_sibl == id) {
        entry(act_sibl).set_next_sibling(entry(next_sibl).next_sibling());
        break;
      }
      act_sibl = next_sibl;
    }
  }

  // 2. add the children of this one as
  //    siblings to the current siblings
  CacheEntryID act_child = entry(id).first_descendant();
  while (act_child) {
    CacheEntryID next_child = entry(act_child).next_sibling();
    entry(act_child).set_father(father);
    entry(act_child).set_next_sibling(entry(father).first_descendant());
    entry(father).set_first_descendant(act_child);
    act_child = next_child;
  }
}

template <class T_num>
void ComponentCache<T_num>::storeValueOf(CacheEntryID id, const T_num &model_count) {
  considerCacheResize();
  unsigned table_ofs = tableEntry(id);
  // when storing the new model count the size of the model count
  // and hence that of the component will change
  statistics_.sum_bytes_cached_components_ -= entry(id).SizeInBytes();
  statistics_.overall_bytes_components_stored_ -= entry(id).SizeInBytes();

  statistics_.sys_overhead_sum_bytes_cached_components_ -= entry(id).sys_overhead_SizeInBytes();
  statistics_.sys_overhead_overall_bytes_components_stored_ -= entry(id).sys_overhead_SizeInBytes();

  entry(id).set_model_count(model_count,my_time_);
  entry(id).set_creation_time(my_time_);

  entry(id).set_next_bucket_element(table_[table_ofs]);
  table_[table_ofs] = id;

  statistics_.sum_bytes_cached_components_ += entry(id).SizeInBytes();
  statistics_.overall_bytes_components_stored_ += entry(id).SizeInBytes();

  statistics_.sys_overhead_sum_bytes_cached_components_ += entry(id).sys_overhead_SizeInBytes();
   statistics_.sys_overhead_overall_bytes_components_stored_ += entry(id).sys_overhead_SizeInBytes();

}

#include <algorithm>

#ifdef __linux__

#include <sys/sysinfo.h>
#include <cstdint>

inline uint64_t freeram() {

  struct sysinfo info;
      sysinfo(&info);

  return info.freeram *(uint64_t) info.mem_unit;
}

#elif __APPLE__ && __MACH__

#include <sys/types.h>
#include <sys/sysctl.h>


inline uint64_t freeram() {

  int mib[2];
  int64_t physical_memory;
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  size_t length = sizeof(int64_t);
  sysctl(mib, 2, &physical_memory, &length, NULL, 0);

  return physical_memory;
}

#else

#endif


template <class T_num>
ComponentCache<T_num>::ComponentCache(DataAndStatistics<T_num> &statistics) :
    statistics_(statistics) {
}

template <class T_num>
void ComponentCache<T_num>::init(Component &super_comp, const Hasher& hasher) {

   // cout << sizeof(CacheableComponent<T_num>) << " " << sizeof(T_num) << endl;
    CacheableComponent<T_num> &packed_super_comp = *new CacheableComponent<T_num>(super_comp, hasher);
  my_time_ = 1;

  entry_base_.clear();
  entry_base_.reserve(2000000);
  entry_base_.push_back(new CacheableComponent<T_num>()); // dummy Element
  table_.clear();
  table_.resize(1024*1024, 0);
  table_size_mask_ = table_.size() - 1;

  free_entry_base_slots_.clear();
  free_entry_base_slots_.reserve(10000);

  uint64_t free_ram = freeram();
  uint64_t max_cache_bound = 95 * (free_ram / 100);

  if (statistics_.maximum_cache_size_bytes_ == 0) {
    statistics_.maximum_cache_size_bytes_ = max_cache_bound;
  }

  if (statistics_.maximum_cache_size_bytes_ > free_ram) {
    cout <<"c o WARNING: Maximum cache size larger than free RAM available" << endl;
    cout << "c o Free RAM " << free_ram / 1000000 << "MB" << endl;
  }

  cout << "c o Maximum cache size:\t"
      << statistics_.maximum_cache_size_bytes_ / 1000000 << " MB" << endl;

  assert(!statistics_.cache_full());

  if (entry_base_.capacity() == entry_base_.size())
    entry_base_.reserve(2 * entry_base_.size());

  entry_base_.push_back(&packed_super_comp);

  statistics_.incorporate_cache_store(packed_super_comp);

  super_comp.set_id(1);
}

template <class T_num>
void ComponentCache<T_num>::test_descendantstree_consistency() {
  for (unsigned id = 2; id < entry_base_.size(); id++)
    if (entry_base_[id] != nullptr) {
      CacheEntryID act_child = entry(id).first_descendant();
      while (act_child) {
        CacheEntryID next_child = entry(act_child).next_sibling();
        assert(entry(act_child).father() == id);

        act_child = next_child;
      }
      CacheEntryID father = entry(id).father();
      CacheEntryID act_sib = entry(father).first_descendant();

      bool found = false;

      while (act_sib) {
        CacheEntryID next_sib = entry(act_sib).next_sibling();
        if (act_sib == id)
          found = true;
        act_sib = next_sib;
      }
      assert(found);
    }
}


template <class T_num>
bool ComponentCache<T_num>::deleteEntries() {
  assert(statistics_.cache_full());

  vector<double> scores;
  for (auto it = entry_base_.begin() + 1; it != entry_base_.end(); it++)
    if (*it != nullptr && (*it)->isDeletable()) {
      scores.push_back((double) (*it)->creation_time());
    }
  sort(scores.begin(), scores.end());
  double cutoff = scores[scores.size() / 2];

  //cout << "cutoff" << cutoff  << " entries: "<< entry_base_.size()<< endl;

  // first : go through the EntryBase and mark the entries to be deleted as deleted (i.e. EMPTY
  // note we start at index 2,
  // since index 1 is the whole formula,
  // should always stay here!
  for (unsigned id = 2; id < entry_base_.size(); id++)
    if (entry_base_[id] != nullptr &&
        entry_base_[id]->isDeletable() &&
          (double) entry_base_[id]->creation_time() <= cutoff) {
        removeFromDescendantsTree(id);
        eraseEntry(id);

        }
  // then go through the Hash Table and erase all Links to empty entries


#ifdef DEBUG
  test_descendantstree_consistency();
#endif

  reHashTable(table_.size());
  statistics_.sum_size_cached_components_ = 0;
  statistics_.sum_bytes_cached_components_ = 0;
   statistics_.sys_overhead_sum_bytes_cached_components_ =0;

  statistics_.sum_bytes_pure_cached_component_data_ = 0;

  for (unsigned id = 2; id < entry_base_.size(); id++)
    if (entry_base_[id] != nullptr) {
      statistics_.sum_size_cached_components_ +=
          entry_base_[id]->num_variables();
      statistics_.sum_bytes_cached_components_ +=
          entry_base_[id]->SizeInBytes();
       statistics_.sys_overhead_sum_bytes_cached_components_ +=
           entry_base_[id]->sys_overhead_SizeInBytes();
    }

  statistics_.num_cached_components_ = entry_base_.size();
  compute_byte_size_infrasture();

  //cout << " \t entries: "<< entry_base_.size() - free_entry_base_slots_.size()<< endl;
  return true;
}

template <class T_num>
uint64_t ComponentCache<T_num>::compute_byte_size_infrasture() {
  statistics_.cache_infrastructure_bytes_memory_usage_ =
      sizeof(ComponentCache)
      + sizeof(CacheEntryID)* table_.capacity()
      + sizeof(CacheableComponent<T_num> *)* entry_base_.capacity()
      + sizeof(CacheEntryID) * free_entry_base_slots_.capacity();
  return statistics_.cache_infrastructure_bytes_memory_usage_;
}

template <class T_num>
void ComponentCache<T_num>::debug_dump_data(){
    cout << "sizeof (CacheableComponent *, CacheEntryID) "
         << sizeof(CacheableComponent<T_num> *) << ", "
         << sizeof(CacheEntryID) << endl;
    cout << "table (size/capacity) " << table_.size()
         << "/" << table_.capacity() << endl;
    cout << "entry_base_ (size/capacity) " << entry_base_.size()
             << "/" << entry_base_.capacity() << endl;
    cout << "free_entry_base_slots_ (size/capacity) " << free_entry_base_slots_.size()
             << "/" << free_entry_base_slots_.capacity() << endl;

//    uint64_t size_model_counts = 0;
    uint64_t alloc_model_counts = 0;
    for (auto &pentry : entry_base_)
              if (pentry != nullptr){
//                size_model_counts += pentry->size_of_model_count();
                alloc_model_counts += pentry->alloc_of_model_count();
              }
    cout << "model counts size " << alloc_model_counts << endl;
}

#endif /* COMPONENT_CACHE_H_ */
