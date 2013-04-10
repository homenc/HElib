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
/* timing.h - utility functions for measuering time
 */
#ifndef _TIMING_H_
#define _TIMING_H_
#include <iostream>

// Activate/deactivate/check-status of timing
extern bool FHEtimersOn;
inline void setTimersOn()  { FHEtimersOn=true; }
inline void setTimersOff() { FHEtimersOn=false; }
inline bool areTimersOn()  { return FHEtimersOn; }

void startFHEtimer(const char *fncName); // start counting
void stopFHEtimer(const char *fncName);  // pause counting
void resetFHEtimer(const char *fncName); // reset value to zero
double getTime4func(const char *fncName); // returns time in seconds
long getNumCalls4func(const char *fncName); // returns # of calls

void resetAllTimers();
void printAllTimers(std::ostream& str=std::cerr);

#define FHE_TIMER_START {if (areTimersOn()) startFHEtimer(__func__);}
#define FHE_TIMER_STOP  {if (areTimersOn()) stopFHEtimer(__func__);}

#define FHE_NTIMER_START(n) {if (areTimersOn()) startFHEtimer(n);}
#define FHE_NTIMER_STOP(n)  {if (areTimersOn()) stopFHEtimer(n);}

#endif // _TIMING_H_
