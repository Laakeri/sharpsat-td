/*
 * difference_packed_component.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

#ifndef DIFFERENCE_PACKED_COMPONENT_H_
#define DIFFERENCE_PACKED_COMPONENT_H_

#include "base_packed_component.h"
#include "component.h"
#include <type_traits>

#include <math.h>

#include "../hasher.h"


template <class T_num>
class DifferencePackedComponent : public BasePackedComponent<T_num> {
public:

  DifferencePackedComponent() {
  }

  inline DifferencePackedComponent(Component &rComp, const Hasher& hasher);

  unsigned num_variables() const{
    return num_variables_;
  }

  unsigned data_size() const {
         return 0;
    }

  unsigned data_only_byte_size() const {
        return data_size()* sizeof(unsigned);
    }

    unsigned raw_data_byte_size() const;

    // raw data size with the overhead
    // for the supposed 16byte alignment of malloc
    unsigned sys_overhead_raw_data_byte_size() const;

  bool equals(const DifferencePackedComponent &comp) const {
    return clhashkey_ == comp.clhashkey();
  }

  HashType clhashkey() const {
    return clhashkey_;
  }

  unsigned hashkey() {
    return clhashkey_[0];
  }

private:
  HashType clhashkey_;
  unsigned num_variables_ = 0;
};

template <>
inline unsigned DifferencePackedComponent<LogNum>::raw_data_byte_size() const {
  return data_size()* sizeof(unsigned);
}

template <>
inline unsigned DifferencePackedComponent<double>::raw_data_byte_size() const {
  return data_size()* sizeof(unsigned);
}

template <>
inline unsigned DifferencePackedComponent<mpz_class>::raw_data_byte_size() const {
  return data_size()* sizeof(unsigned)
           + BasePackedComponent<mpz_class>::model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
}

template<>
inline unsigned DifferencePackedComponent<LogNum>::sys_overhead_raw_data_byte_size() const {
  unsigned ds = data_size()* sizeof(unsigned);
  unsigned ms = 0;
//      unsigned mask = 0xfffffff8;
//      return (ds & mask) + ((ds & 7)?8:0)
//            +(ms & mask) + ((ms & 7)?8:0);
  unsigned mask = 0xfffffff0;
        return (ds & mask) + ((ds & 15)?16:0)
              +(ms & mask) + ((ms & 15)?16:0);
}

template<>
inline unsigned DifferencePackedComponent<double>::sys_overhead_raw_data_byte_size() const {
  unsigned ds = data_size()* sizeof(unsigned);
  unsigned ms = 0;
//      unsigned mask = 0xfffffff8;
//      return (ds & mask) + ((ds & 7)?8:0)
//            +(ms & mask) + ((ms & 7)?8:0);
  unsigned mask = 0xfffffff0;
        return (ds & mask) + ((ds & 15)?16:0)
              +(ms & mask) + ((ms & 15)?16:0);
}

template<>
inline unsigned DifferencePackedComponent<mpz_class>::sys_overhead_raw_data_byte_size() const {
  unsigned ds = data_size()* sizeof(unsigned);
  unsigned ms = BasePackedComponent<mpz_class>::model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
//      unsigned mask = 0xfffffff8;
//      return (ds & mask) + ((ds & 7)?8:0)
//            +(ms & mask) + ((ms & 7)?8:0);
  unsigned mask = 0xfffffff0;
        return (ds & mask) + ((ds & 15)?16:0)
              +(ms & mask) + ((ms & 15)?16:0);
}

template <class T_num>
DifferencePackedComponent<T_num>::DifferencePackedComponent(Component &rComp, const Hasher& hasher) {
  num_variables_ = rComp.num_variables();
  clhashkey_ = hasher.Hash(rComp.RawData());
}

#endif /* DIFFERENCE_PACKED_COMPONENT_H_ */
