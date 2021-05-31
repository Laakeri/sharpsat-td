#ifndef HASHER_H_
#define HASHER_H_

#include "clhash/clhash.h"

#include <array>
#include <cstdio>
#include <random>
#include <cassert>

using namespace std;

const size_t HashM = 2;
using HashType = std::array<uint64_t, HashM>;

class Hasher {
 public:
 	inline Hasher(std::mt19937_64& gen) {
 		assert(HashM >= 1);
 		for (size_t i = 0; i < HashM; i++) {
 			random_data_[i] = get_random_key_for_clhash(gen(), gen());
 		}
 	}
 	inline HashType Hash(const vector<unsigned>& data) const {
 		HashType ret;
 		for (size_t i = 0; i < HashM; i++) {
 			ret[i] = clhash(random_data_[i], (const char *)data.data(), data.size() * sizeof(unsigned));
 		}
 		return ret;
 	}
 	inline ~Hasher() {
 		for (size_t i = 0; i < HashM; i++) {
 			std::free((void*)random_data_[i]);
 		}
 	}
 private:
  const void* random_data_[HashM];
};

#endif /* HASHER_H_ */