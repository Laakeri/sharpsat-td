/*
 * structures.h
 *
 *  Created on: Jun 25, 2012
 *      Author: Marc Thurley
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <cassert>
#include "primitive_types.h"
#include <gmpxx.h>
#include "mpfr/mpreal.h"
using namespace std;

struct SDouble {
 public:
  SDouble() {
    n = 0;
    has = false;
  }
  void Init(double d) {
    assert(d != 0);
    n = d;
    has = true;
  }
  SDouble(const SDouble& other) {
    n = other.n;
    has = other.has;
  }
  SDouble& operator=(const SDouble& other) {
    n = other.n;
    has = other.has;
    return *this;
  }
  bool IsAlgZero() const {
    return !has;
  }
  SDouble operator*(SDouble other) const {
    SDouble ret = other;
    ret.n *= n;
    ret.has &= has;
    return ret;
  }
  SDouble operator+(SDouble other) const {
    SDouble ret = other;
    ret.n += n;
    ret.has |= has;
    return ret;
  }
  SDouble& operator*=(const SDouble& other) {
    n *= other.n;
    has &= other.has;
    return *this;
  }
  SDouble& operator/=(const SDouble& other) {
    assert(other.n != 0);
    assert(other.has);
    n /= other.n;
    return *this;
  }
  size_t InternalSize() const {
    return 0;
  }
  double Get() const {
    return n;
  }
  static SDouble Zero() {
    SDouble ret;
    return ret;
  }
  static SDouble One() {
    SDouble ret;
    ret.n = 1;
    ret.has = true;
    return ret;
  }
 private:
  double n = 0;
  bool has = false;
};

struct Smpr {
 public:
  Smpr() {
    n = 0;
    has = false;
  }
  void Init(double d) {
    assert(d != 0);
    n = d;
    has = true;
  }
  Smpr(const Smpr& other) {
    n = other.n;
    has = other.has;
  }
  Smpr& operator=(const Smpr& other) {
    n = other.n;
    has = other.has;
    return *this;
  }
  bool IsAlgZero() const {
    return !has;
  }
  Smpr operator*(Smpr other) const {
    Smpr ret = other;
    ret.n *= n;
    ret.has &= has;
    return ret;
  }
  Smpr operator+(Smpr other) const {
    Smpr ret = other;
    ret.n += n;
    ret.has |= has;
    return ret;
  }
  Smpr& operator*=(const Smpr& other) {
    n *= other.n;
    has &= other.has;
    return *this;
  }
  Smpr& operator/=(const Smpr& other) {
    assert(other.n != 0);
    assert(other.has);
    n /= other.n;
    return *this;
  }
  size_t InternalSize() const {
    return 0;
  }
  mpfr::mpreal Get() const {
    return n;
  }
  static Smpr Zero() {
    Smpr ret;
    return ret;
  }
  static Smpr One() {
    Smpr ret;
    ret.n = 1;
    ret.has = true;
    return ret;
  }
 private:
  mpfr::mpreal n = 0;
  bool has = false;
};

struct Smpz {
 public:
  Smpz() {
    n = 0;
    has = false;
  }
  void Init(double d) {
    assert(d == 1 || d == -1);
    n = (int)d;
    has = true;
  }
  Smpz(const Smpz& other) {
    n = other.n;
    has = other.has;
  }
  Smpz& operator=(const Smpz& other) {
    n = other.n;
    has = other.has;
    return *this;
  }
  bool IsAlgZero() const {
    return !has;
  }
  Smpz operator*(Smpz other) const {
    Smpz ret = other;
    ret.n *= n;
    ret.has &= has;
    return ret;
  }
  Smpz operator+(Smpz other) const {
    Smpz ret = other;
    ret.n += n;
    ret.has |= has;
    return ret;
  }
  Smpz& operator*=(const Smpz& other) {
    n *= other.n;
    has &= other.has;
    return *this;
  }
  Smpz& operator/=(const Smpz& other) {
    assert(other.has);
    if (other.n == -1) {
      n = -n;
    } else {
      assert(other.n == 1);
    }
    return *this;
  }
  size_t InternalSize() const {
    return n.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
  }
  mpz_class Get() const {
    return n;
  }
  static Smpz Zero() {
    Smpz ret;
    return ret;
  }
  static Smpz One() {
    Smpz ret;
    ret.n = 1;
    ret.has = true;
    return ret;
  }
 private:
  mpz_class n = 0;
  bool has = false;
};

#define INVALID_DL -1

typedef unsigned char TriValue;
#define   F_TRI  0
#define   T_TRI  1
#define   X_TRI  2

class LiteralID {
public:

  LiteralID() {
    value_ = 0;
  }
  LiteralID(int lit) {
    value_ = (abs(lit) << 1) + (unsigned) (lit > 0);
  }

  LiteralID(VariableIndex var, bool sign) {
    value_ = (var << 1) + (unsigned) sign;
  }

  VariableIndex var() const {
    return (value_ >> 1);
  }

  int toInt() const {
    return ((int) value_ >> 1) * ((sign()) ? 1 : -1);
  }

  void inc(){++value_;}

  void copyRaw(unsigned int v) {
    value_ = v;
  }

  bool sign() const {
    return (bool) (value_ & 0x01);
  }

  bool operator!=(const LiteralID &rL2) const {
    return value_ != rL2.value_;
  }

  bool operator==(const LiteralID &rL2) const {
    return value_ == rL2.value_;
  }

  const LiteralID neg() const {
    return LiteralID(var(), !sign());
  }

  void print() const {
    cout << (sign() ? " " : "-") << var() << " ";
  }

  unsigned raw() const { return value_;}

private:
  unsigned value_;

  template <class _T> friend class LiteralIndexedVector;
};

static const LiteralID NOT_A_LIT(0, false);
#define SENTINEL_LIT NOT_A_LIT

class Literal {
public:
  vector<LiteralID> binary_links_ = vector<LiteralID>(1,SENTINEL_LIT);
  vector<ClauseOfs> watch_list_ = vector<ClauseOfs>(1,SENTINEL_CL);
  float activity_score_ = 0.0f;

  void increaseActivity(unsigned u = 1){
    activity_score_+= u;
  }

  void removeWatchLinkTo(ClauseOfs clause_ofs) {
    for (auto it = watch_list_.begin(); it != watch_list_.end(); it++)
          if (*it == clause_ofs) {
            *it = watch_list_.back();
            watch_list_.pop_back();
            return;
          }
  }

  void replaceWatchLinkTo(ClauseOfs clause_ofs, ClauseOfs replace_ofs) {
        for (auto it = watch_list_.begin(); it != watch_list_.end(); it++)
          if (*it == clause_ofs) {
            *it = replace_ofs;
            return;
          }
  }

  void addWatchLinkTo(ClauseIndex clause_ofs) {
    watch_list_.push_back(clause_ofs);
  }

  void addBinLinkTo(LiteralID lit) {
    binary_links_.back() = lit;
    binary_links_.push_back(SENTINEL_LIT);
  }

  void resetWatchList(){
        watch_list_.clear();
        watch_list_.push_back(SENTINEL_CL);
  }

  bool hasBinaryLinkTo(LiteralID lit) {
    for (auto l : binary_links_) {
      if (l == lit)
        return true;
    }
    return false;
  }

  bool hasBinaryLinks() {
    return !binary_links_.empty();
  }
};

class Antecedent {
  unsigned int val_;

public:
  Antecedent() {
    val_ = 1;
  }

  Antecedent(const ClauseOfs cl_ofs) {
     val_ = (cl_ofs << 1) | 1;
   }
  Antecedent(const LiteralID idLit) {
    val_ = (idLit.raw() << 1);
  }

  bool isAClause() const {
    return val_ & 0x01;
  }

  ClauseOfs asCl() const {
      return val_ >> 1;
    }

  LiteralID asLit() {
    LiteralID idLit;
    idLit.copyRaw(val_ >> 1);
    return idLit;
  }
  // A NON-Antecedent will only be A NOT_A_CLAUSE Clause Id
  bool isAnt() {
    return val_ != 1; //i.e. NOT a NOT_A_CLAUSE;
  }
};


struct Variable {
  Antecedent ante;
  int decision_level = INVALID_DL;
};

// for now Clause Header is just a dummy
// we keep it for possible later changes
class ClauseHeader {
  unsigned creation_time_; // number of conflicts seen at creation time
  unsigned glue_;
  unsigned length_;
public:

  void setGlue(unsigned glue) {
    glue_ = glue;
  }
  void see() {
    glue_ += 3;
  }
  unsigned glue() const {
      return glue_;
  }
  bool isLearned() {
    return glue_ >= 1;
  }
  unsigned creation_time() {
      return creation_time_;
  }
  unsigned length(){ return length_;}
  void set_length(unsigned length){ length_ = length;}

  void set_creation_time(unsigned time) {
    creation_time_ = time;
  }
  static unsigned overheadInLits(){return sizeof(ClauseHeader)/sizeof(LiteralID);}
};

#endif /* STRUCTURES_H_ */
