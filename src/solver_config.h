/*
 * basic_types.h
 *
 *  Created on: Jun 24, 2012
 *      Author: Marc Thurley
 */

#ifndef SOLVER_CONFIG_H_
#define SOLVER_CONFIG_H_


struct SolverConfiguration {

  bool perform_non_chron_back_track = true;

  bool verbose = false;

  // quiet = true will override verbose;
  bool quiet = false;

  // The following settings are treewidth mod
  double decomp_weight = -1;
  int weight_mode = 1;
  // Added settings for experimenting with VSADS scores
  bool vsads_freq = true;
  bool vsads_act = true;
  bool vsads_act_init = true;
};

#endif /* SOLVER_CONFIG_H_ */
