#include "solver.h"
#include "preprocessor/instance.hpp"
#include "preprocessor/preprocessor.hpp"
#include "preprocessor/treewidth.hpp"

#include "solver_config.h"

#include <iostream>

#include <vector>
#include <limits>

#include <string>

#include <sys/time.h>
#include <sys/resource.h>
#include <gmpxx.h>
#include <limits>
#include "mpfr/mpreal.h"

#include <random>

using namespace std;


void PrintSat(bool sat) {
  if (sat) {
    cout<<"s SATISFIABLE"<<endl;
  } else {
    cout<<"s UNSATISFIABLE"<<endl;
  }
}

void PrintType(const sspp::Instance& ins) {
  if (ins.weighted) {
    cout<<"c s type wmc"<<endl;
  } else {
    cout<<"c s type mc"<<endl;
  }
}

mpfr::mpreal Log10(const mpz_class& num) {
  assert(num >= 0);
  if (num == 0) {
    return -std::numeric_limits<double>::infinity();
  }
  mpfr::mpreal num1(num.get_mpz_t());
  return mpfr::log10(num1);
}

void PrintLog10(const mpz_class& num) {
  cout<<"c s log10-estimate "<<Log10(num)<<endl;
}

void PrintLog10(double num, double logwf) {
  cout<<"c s log10-estimate "<<log10(num)+logwf<<endl;
}

void PrintLog10(const mpfr::mpreal& num) {
  cout<<"c s log10-estimate "<<mpfr::log10(num)<<endl;
}

void PrintExact(const mpz_class& num) {
  cout<<"c s exact arb int "<<num<<endl;
}

void PrintExact(const mpfr::mpreal& num) {
  cout<<"c s exact arb float "<<num<<endl;
}

void PrintDouble(double num) {
  cout<<"c s exact double float "<<num<<endl;
}

int main(int argc, char *argv[]) {
  cout<<std::setprecision(16);
  sspp::Timer glob_timer;
  glob_timer.start();
  string input_file;

  // Randomness used only for component hashing
  std::mt19937_64 gen(1337);

  string tmp_dir;
  double decot = -1;

  SolverConfiguration config_;

  uint64_t max_cache = 0;

  int weighted = 0;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-WD") == 0) {
      assert(weighted == 0);
      weighted = 1;
    } else if (strcmp(argv[i], "-WE") == 0) {
      assert(weighted == 0);
      weighted = 2;
    } else if (strcmp(argv[i], "-tmpdir") == 0) {
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      tmp_dir = string(argv[i+1]);
    } else if (strcmp(argv[i], "-cs") == 0) {
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      max_cache = atol(argv[i + 1]) * (uint64_t) 1000000;
      // theSolver.statistics().maximum_cache_size_bytes_ = atol(argv[i + 1]) * (uint64_t) 1000000;
    } else if (strcmp(argv[i], "-nofreq") == 0) {
      assert(config_.vsads_freq == true);
      config_.vsads_freq = false;
    } else if (strcmp(argv[i], "-noact") == 0) {
      assert(config_.vsads_act == true);
      config_.vsads_act = false;
    } else if (strcmp(argv[i], "-noactinit") == 0) {
      assert(config_.vsads_act_init == true);
      config_.vsads_act_init = false;
    } else if (strcmp(argv[i], "-decot") == 0) {
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      decot = atof(argv[++i]);
    } else if (strcmp(argv[i], "-decow") == 0) {
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      double decow = atof(argv[++i]);
      assert(decow > 0.001);
      config_.decomp_weight = decow;
    } else if (strcmp(argv[i], "-wemod") == 0) {
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      int wemod = atoi(argv[++i]);
      assert(wemod >= 1 && wemod <= 3);
      config_.weight_mode = wemod;
    } else if (strcmp(argv[i], "-prec") == 0) {
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      int prec = atoi(argv[++i]);
      assert(prec >= 1);
      cout<<std::setprecision(prec);
    } else {
      input_file = argv[i];
    }
  }

  assert(!tmp_dir.empty());
  assert(decot > 0.0001 && decot < 10000);

  if (weighted == 0) {
    sspp::Instance ins(input_file, false);
    sspp::Preprocessor ppp;
    ppp.SetMaxGTime(150);
    ppp.SetMaxSparsTime(120);
    ins = ppp.Preprocess(ins, "FPVSEGV");
    ins.UpdClauseInfo();
    cout<<"c o Preprocessed. "<<glob_timer.get()<<"s Vars: "<<ins.vars<<" Clauses: "<<ins.clauses.size()<<" Free vars: "<<ppp.FreeVars()<<endl;
    if (ins.vars == 1) {
      assert(ins.clauses.size() == 2);
      PrintSat(false);
      PrintType(ins);
      PrintLog10((mpz_class)0);
      PrintExact((mpz_class)0);
      return 0;
    }
    mpz_class ans0 = sspp::Power2<mpz_class>(ppp.FreeVars());
    if (ins.vars == 0) {
      PrintSat(true);
      PrintType(ins);
      PrintLog10(ans0);
      PrintExact(ans0);
      return 0;
    }
    sspp::Graph primal(ins.vars, ins.clauses);
    sspp::TreeDecomposition tdecomp = sspp::decomp::Treedecomp(primal, decot, tmp_dir);
    cout<<"c o Now solving. "<<glob_timer.get()<<endl;
    Solver<Smpz> theSolver(gen);
    theSolver.config() = config_;
    if (max_cache > 0) {
      theSolver.statistics().maximum_cache_size_bytes_ = max_cache;
    }
    mpz_class ans = theSolver.solve(ins, tdecomp).Get();
    cout<<"c o Solved. "<<glob_timer.get()<<endl;
    ans *= ans0;
    PrintSat(true);
    PrintType(ins);
    PrintLog10(ans);
    PrintExact(ans);
    return 0;
  } else if (weighted == 1 || weighted == 2) {
    sspp::Instance ins(input_file, true);
    sspp::Preprocessor ppp;
    ins = ppp.Preprocess(ins, "FPVE");
    ins.UpdClauseInfo();
    cout<<"c o Preprocessed. "<<glob_timer.get()<<"s Vars: "<<ins.vars<<" Clauses: "<<ins.clauses.size()<<" Free vars: "<<ppp.FreeVars()<<endl;
    if (ins.vars == 1 && ins.clauses.size() == 2) {
      PrintSat(false);
      PrintType(ins);
      PrintLog10((mpz_class)0);
      PrintExact((mpfr::mpreal)0);
      return 0;
    }
    mpfr::mpreal ans0 = ins.weight_factor;
    cout<<"c o wf "<<ans0<<endl;
    if (ins.vars == 0) {
      PrintSat(true);
      PrintType(ins);
      PrintLog10(ans0);
      PrintExact(ans0);
      return 0;
    }
    sspp::Graph primal(ins.vars, ins.clauses);
    sspp::TreeDecomposition tdecomp = sspp::decomp::Treedecomp(primal, decot, tmp_dir);
    cout<<"c o Now solving. "<<glob_timer.get()<<endl;
    if (weighted == 1) {
      Solver<SDouble> theSolver(gen);
      theSolver.config() = config_;
      if (max_cache > 0) {
        theSolver.statistics().maximum_cache_size_bytes_ = max_cache;
      }
      double ans1 = theSolver.solve(ins, tdecomp).Get();
      cout<<"c o Solved. "<<glob_timer.get()<<endl;
      PrintSat(true);
      PrintType(ins);
      PrintLog10(ans1, (double)mpfr::log10(ans0));
      PrintDouble(ans1*(double)ans0);
    } else {
      Solver<Smpr> theSolver(gen);
      theSolver.config() = config_;
      if (max_cache > 0) {
        theSolver.statistics().maximum_cache_size_bytes_ = max_cache;
      }
      mpfr::mpreal ans1 = theSolver.solve(ins, tdecomp).Get();
      cout<<"c o Solved. "<<glob_timer.get()<<endl;
      PrintSat(true);
      PrintType(ins);
      PrintLog10(ans1*ans0);
      PrintExact(ans1*ans0);
    }
    return 0;
  } else {
    assert(0);
  }
}