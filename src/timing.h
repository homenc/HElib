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
#include <iostream>

// Activate/deactivate/check-status of timing
extern bool FHEtimersOn;
inline void setTimersOn()  { FHEtimersOn=true; }
inline void setTimersOff() { FHEtimersOn=false; }
inline bool areTimersOn()  { return FHEtimersOn; }

//! Start a timer
void startFHEtimer(const char *fncName); // start counting
//! Stop a timer
void stopFHEtimer(const char *fncName);  // pause counting
//! Reset a timer for some label to zero
void resetFHEtimer(const char *fncName); // reset value to zero
//! Read the value of a timer (in seconds)
double getTime4func(const char *fncName); // returns time in seconds
//! Returns number of calls for that timer
long getNumCalls4func(const char *fncName); // returns # of calls

void resetAllTimers();
//! Print the value of all timers to stream
void printAllTimers(std::ostream& str=std::cerr);

#define FHE_TIMER_START {if (areTimersOn()) startFHEtimer(__func__);}
#define FHE_TIMER_STOP  {if (areTimersOn()) stopFHEtimer(__func__);}

#define FHE_NTIMER_START(n) {if (areTimersOn()) startFHEtimer(n);}
#define FHE_NTIMER_STOP(n)  {if (areTimersOn()) stopFHEtimer(n);}

#endif // _TIMING_H_
