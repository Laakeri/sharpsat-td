/*
 * solver.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: marc
 */
#include "solver.h"

#include <algorithm>


StopWatch::StopWatch() {
  interval_length_.tv_sec = 60;
  gettimeofday(&last_interval_start_, NULL);
  start_time_ = stop_time_ = last_interval_start_;
}

timeval StopWatch::getElapsedTime() {
  timeval r;
  timeval other_time = stop_time_;
  if (stop_time_.tv_sec == start_time_.tv_sec
      && stop_time_.tv_usec == start_time_.tv_usec)
    gettimeofday(&other_time, NULL);
  long int ad = 0;
  long int bd = 0;

  if (other_time.tv_usec < start_time_.tv_usec) {
    ad = 1;
    bd = 1000000;
  }
  r.tv_sec = other_time.tv_sec - ad - start_time_.tv_sec;
  r.tv_usec = other_time.tv_usec + bd - start_time_.tv_usec;
  return r;
}