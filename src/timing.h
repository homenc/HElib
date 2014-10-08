/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
/**
 * @file timing.h
 * @brief Utility functions for measuering time
 * 
 * This module contains some utility functions for measuring the time that
 * various methods take to execute. To use it, we insert the macro
 * FHE_TIMER_START at the beginning of the method(s) that we want to time and
 * FHE_TIMER_STOP at the end, then the main program needs to call the function
 * setTimersOn() to activate the timers and setTimersOff() to pause them. 
 * To obtain the value of a given timer (in seconds), the application can
 * use the function getTime4func(const char *fncName), and the function
 * printAllTimers() prints the values of all timers to an output stream.
 *
 * Using this method we can have at most one timer per method/function, and
 * the timer is called by the same name as the function itself (using the
 * built-in macro \_\_func\_\_). We can also use the "lower level" methods
 * startFHEtimer(name), stopFHEtimer(name), and resetFHEtimer(name) to add
 * timers with arbitrary names (not necessarily associated with functions).
 **/
#ifndef _TIMING_H_
#define _TIMING_H_
#include <ctime>


//! A simple class to toggle timing information on and off
class FHEtimer {
public:
  const char *name;
  const char *loc;
  bool isOn;  // a broken semaphore
  std::clock_t counter;
  long numCalls;

  FHEtimer(const char *_name, const char *_loc) :
    name(_name), loc(_loc), isOn(false), counter(0), numCalls(0) 
  { }

  void reset();
  void start();
  void stop();
  double getTime() const;
  long getNumCalls() const;
};


// Activate/deactivate/check-status of timing
extern bool FHEtimersOn;

inline void setTimersOn()  { FHEtimersOn=true; }
inline void setTimersOff() { FHEtimersOn=false; }
inline bool areTimersOn()  { return FHEtimersOn; }

void registerTimer(FHEtimer *timer);

inline void buildTimer(FHEtimer*& timer, const char *name, const char *loc)
{
  if (areTimersOn() && !timer) {
    FHEtimer *tmp = new FHEtimer(name, loc);
    registerTimer(tmp);
    timer = tmp;
  }
}



void resetAllTimers();
//! Print the value of all timers to stream
void printAllTimers(std::ostream& str=std::cerr);

// return true if timer was found, false otherwise
bool getTimerByName(FHEtimer& timer, const char* name);
bool printNamedTimer(ostream& str, const char* name);


class auto_timer {
public:
  FHEtimer *timer;
  auto_timer(FHEtimer *_timer) : timer(_timer) { if (timer && areTimersOn()) timer->start(); }
  void stop() { if (timer && areTimersOn()) timer->stop(); }
  ~auto_timer() { stop(); }
};


// NOTE: the STOP functions below are not really needed,
// but are provided for backward compatibility

#define FHE_STRINGIFY(x) #x
#define FHE_TOSTRING(x) FHE_STRINGIFY(x)
#define FHE_AT __FILE__ ":" FHE_TOSTRING(__LINE__)

#define FHE_stringify_aux(s) #s
#define FHE_stringify(s) FHE_stringify_aux(s)

#define FHE_TIMER_START \
  static FHEtimer *_local_timer = 0; \
  buildTimer(_local_timer, __func__, FHE_AT ); \
  auto_timer _local_auto_timer(_local_timer)
  
#define FHE_TIMER_STOP  _local_auto_timer.stop()


#define FHE_NTIMER_START(n) \
  static FHEtimer *_named_local_timer ## n = 0; \
  buildTimer(_named_local_timer ## n, # n, FHE_AT ); \
  auto_timer _named_local_auto_timer ## n(_named_local_timer ## n)

#define FHE_NTIMER_STOP(n)    _named_local_auto_timer ## n.stop();

#endif // _TIMING_H_
