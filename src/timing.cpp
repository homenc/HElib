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
#include "timing.h"

#include <algorithm>
#include <utility>
#include <cstring>


#if 1

static
clock_t GetClock()
{
  timespec ts;

  clock_gettime(CLOCK_MONOTONIC, &ts);

  return ts.tv_sec*1000000L + ts.tv_nsec/1000L; 
}

const long CLOCK_SCALE = 1000000L;

#else

static
clock_t GetClock()
{
  return clock();
}

const long CLOCK_SCALE = CLOCKS_PER_SEC;

#endif


bool FHEtimersOn=false;

bool timer_compare(const FHEtimer *a, const FHEtimer *b)
{
  return strcmp(a->name, b->name) < 0;
}



static vector<FHEtimer *> timerMap;

void registerTimer(FHEtimer *timer)
{
  timerMap.push_back(timer);
}

// Reset a timer for some label to zero
void FHEtimer::reset()
{
  numCalls = 0;
  counter = 0;
  if (isOn) counter -= GetClock();
}

// Start a timer
void FHEtimer::start()
{
  if (!isOn) {
    isOn = true;
    numCalls++;
    counter -= GetClock();
  }
}

// Stop a timer
void FHEtimer::stop()
{
  if (isOn) {
    isOn = false;
    counter += GetClock();
  }
}

// Read the value of a timer (in seconds)
double FHEtimer::getTime() const // returns time in seconds
{
  // If the counter is currently counting, add the clock() value
  clock_t c = isOn? (counter + GetClock()) : counter;
  return ((double)c)/CLOCK_SCALE;
}

// Returns number of calls for that timer
long FHEtimer::getNumCalls() const
{
    return numCalls;
}

void resetAllTimers()
{
  for (long i = 0; i < long(timerMap.size()); i++) 
    timerMap[i]->reset();
}

// Print the value of all timers to stream
void printAllTimers(ostream& str)
{

  sort(timerMap.begin(), timerMap.end(), timer_compare);

  for (long i = 0; i < long(timerMap.size()); i++) {
    const char *name = timerMap[i]->name;
    const char *loc = timerMap[i]->loc;
    double t =  timerMap[i]->getTime();
    long n = timerMap[i]->getNumCalls();
    double ave;
    if (n > 0) { 
      ave = t/n;
    }
    else {
      continue;
    }

    str << "  " << name << ": " << t << " / " << n << " = " << ave << "   [" << loc << "]\n";
  }
}

bool getTimerByName(FHEtimer& timer, const char* name)
{
  for (long i = 0; i < long(timerMap.size()); i++) {
    if (strcmp(name, timerMap[i]->name) == 0) {
      timer = *timerMap[i];
      return true;
    }
  }
  timer.numCalls = timer.counter = 0;
  return false;
}

bool printNamedTimer(ostream& str, const char* name)
{
  FHEtimer timer(NULL,NULL);
  getTimerByName(timer, name);
  long n = timer.getNumCalls();
  if (n>0) {
    double t = timer.getTime();
    double ave = t/n;

    str << "  " << name << ": " << t << " / " << n << " = " 
	<< ave << "   [" << timer.loc << "]\n";
    return true;
  }
  return false;
}
