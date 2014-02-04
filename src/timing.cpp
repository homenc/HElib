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
#include <ctime>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>
#include <cstring>

using namespace std;
#if (__cplusplus>199711L)
#include <unordered_map>
#else
#include <tr1/unordered_map>
using namespace tr1;
#endif

#include "timing.h"

//! A simple class to toggle timing information on and off
class FHEtimer {
public:
  bool isOn;  // a broken semaphore
  clock_t counter;

  long numCalls;

  FHEtimer() { isOn=false; counter=0; numCalls=0; }
};

bool string_compare(const char *a, const char *b)
{
  return strcmp(a, b) < 0;
}

bool FHEtimersOn=false;

typedef unordered_map<const char*,FHEtimer>timerMap;
static timerMap timers;

// Reset a timer for some label to zero
void resetFHEtimer(const char *fncName)
{
  FHEtimer& t = timers[fncName];   // insert to map if not already there
  t.numCalls = 0;
  t.counter = 0;
  if (t.isOn) t.counter -= std::clock();
}

// Start a timer
void startFHEtimer(const char *fncName)
{
  FHEtimer& t = timers[fncName];   // insert to map if not already there
  if (!t.isOn) {
    t.isOn = true;
    t.numCalls++;
    t.counter -= std::clock();
  }
}

// Stop a timer
void stopFHEtimer(const char *fncName)
{
  FHEtimer& t = timers[fncName];   // insert to map if not already there
  if (t.isOn) {
    t.isOn = false;
    t.counter += std::clock();
  }
}

// Read the value of a timer (in seconds)
double getTime4func(const char *fncName) // returns time in seconds
{
  FHEtimer& t = timers[fncName];   // insert to map if not already there

  // If the counter is currently counting, add the clock() value
  clock_t c = t.isOn? (t.counter + std::clock()) : t.counter;
  return ((double)c)/CLOCKS_PER_SEC;
}

// Returns number of calls for that timer
long getNumCalls4func(const char *fncName) 
{
    FHEtimer& t = timers[fncName];   // insert to map if not already there
    return t.numCalls;
}

void resetAllTimers()
{
  for (timerMap::iterator it = timers.begin(); it != timers.end(); ++it)
    resetFHEtimer(it->first);
}

// Print the value of all timers to stream
void printAllTimers(ostream& str)
{
  vector<const char *> vec;
  for (timerMap::iterator it = timers.begin(); it != timers.end(); ++it) {
    vec.push_back(it->first);
  }

  sort(vec.begin(), vec.end(), string_compare);

  for (vector<const char *>::iterator it = vec.begin(); it != vec.end(); ++it) {
    double t =  getTime4func(*it);
    long n = getNumCalls4func(*it);
    double ave;
    if (n > 0) { 
      ave = t/n;
    }
    else {
      continue;
    }

    str << "  " << (*it) << ": " << t << " / " << n << " = " << ave << "\n";
  }
}
