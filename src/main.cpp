#include "solver.h"
#include "preprocessor/instance.hpp"
#include "preprocessor/preprocessor.hpp"
#include "preprocessor/treewidth.hpp"

#include "solver_config.h"

#include <iostream>

#include <vector>
#include <limits>

//#include <malloc.h>
#include <string>

#include <sys/time.h>
#include <sys/resource.h>
#include <gmpxx.h>

#include <random>

using namespace std;

#ifdef WEIGHTED
const bool weighted = true;
using NumType = double;
#else
const bool weighted = false;
using NumType = mpz_class;
#endif


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

void PrintLog10(mpz_class num) {
  if (num == 0) {
    cout<<"c s log10-estimate -inf"<<endl;
    return;
  }
  long double out = 0;
  mpz_class billion = 1000000000;
  mpz_class billion4 = billion * billion * billion * billion;
  while (num > billion4) {
    out += 9;
    num /= billion;
  }
  double dnum = num.get_d();
  out += log10(dnum);
  cout<<"c s log10-estimate "<<std::setprecision(16)<<(out)<<endl;
}

void PrintLog10(double num, double logwf) {
  cout<<"c s log10-estimate "<<std::setprecision(16)<<log10(num)+logwf<<endl;
}

void PrintExact(const mpz_class& num) {
  cout<<"c s exact arb int "<<num<<endl;
}

void PrintDouble(double num) {
  cout<<"c s exact double float "<<std::setprecision(16)<<num<<endl;
}

int main(int argc, char *argv[]) {
  sspp::Timer glob_timer;
  glob_timer.start();
  string input_file;

  // Randomness used only for component hashing
  std::mt19937_64 gen(1337);

  string tmp_dir;
  double decot = -1;

  SolverConfiguration config_;

  uint64_t max_cache = 0;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-tmpdir") == 0) {
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
      assert(decow > 0.09);
      config_.decomp_weight = decow;
    } else if (strcmp(argv[i], "-wemod") == 0) {
      if (argc <= i + 1) {
        cout << " wrong parameters" << endl;
        return -1;
      }
      int wemod = atoi(argv[++i]);
      assert(wemod >= 1 && wemod <= 3);
      config_.weight_mode = wemod;
    } else
      input_file = argv[i];
  }

  assert(!tmp_dir.empty());
  assert(decot > 0.0001 && decot < 10000);

  sspp::Instance ins(input_file, weighted);

  #ifndef WEIGHTED
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
    PrintLog10(0);
    PrintExact(0);
    return 0;
  }
  NumType ans0 = sspp::Power2<NumType>(ppp.FreeVars());
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
  Solver<NumType> theSolver(gen);
  theSolver.config() = config_;
  if (max_cache > 0) {
    theSolver.statistics().maximum_cache_size_bytes_ = max_cache;
  }
  NumType ans = theSolver.solve(ins, tdecomp);
  cout<<"c o Solved. "<<glob_timer.get()<<endl;
  ans *= ans0;
  PrintSat(true);
  PrintType(ins);
  PrintLog10(ans);
  PrintExact(ans);
  return 0;
  #else
  sspp::Preprocessor ppp;
  ins = ppp.Preprocess(ins, "FPVE");
  ins.UpdClauseInfo();
  cout<<"c o Preprocessed. "<<glob_timer.get()<<"s Vars: "<<ins.vars<<" Clauses: "<<ins.clauses.size()<<" Free vars: "<<ppp.FreeVars()<<endl;
  if (ins.vars == 1 && ins.clauses.size() == 2) {
    PrintSat(false);
    PrintType(ins);
    PrintLog10((double)0);
    PrintDouble(0);
    return 0;
  }
  double ans0 = ins.weight_factor;
  double ans0log = ins.weight_factor_log;
  cout<<"c o wf "<<ans0<<" "<<ans0log<<endl;
  if (ins.vars == 0) {
    PrintSat(true);
    PrintType(ins);
    PrintLog10(1, ans0log);
    PrintDouble(ans0);
    return 0;
  }
  sspp::Graph primal(ins.vars, ins.clauses);
  sspp::TreeDecomposition tdecomp = sspp::decomp::Treedecomp(primal, decot, tmp_dir);
  cout<<"c o Now solving. "<<glob_timer.get()<<endl;
  NumType ans1;
  {
    Solver<NumType> theSolver(gen);
    theSolver.config() = config_;
    if (max_cache > 0) {
      theSolver.statistics().maximum_cache_size_bytes_ = max_cache;
    }
    ans1 = theSolver.solve(ins, tdecomp);
  }
  cout<<"c o Solved. "<<glob_timer.get()<<endl;
  // ans *= ans0;
  if (ans1 < 1e-200) {
    cout<<"c o Too small"<<endl;
    Solver<LogNum> theSolver2(gen);
    theSolver2.config() = config_;
    if (max_cache > 0) {
      theSolver2.statistics().maximum_cache_size_bytes_ = max_cache;
    }
    LogNum ans2 = theSolver2.solve(ins, tdecomp);
    cout<<"c o Solved again. "<<glob_timer.get()<<endl;
    PrintSat(true);
    PrintType(ins);
    PrintLog10(1, ans0log + ans2.get());
    PrintDouble(ans1*ans0);
  } else {
    PrintSat(true);
    PrintType(ins);
    PrintLog10(ans1, ans0log);
    PrintDouble(ans1*ans0);
  }
  return 0;
  #endif
}