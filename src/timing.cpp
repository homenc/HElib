/* Copyright (C) 2012-2017 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */
#include <algorithm>
#include <utility>
#include <cstring>
#include <ctime>
#include "timing.h"

#ifdef CLOCK_MONOTONIC
unsigned long GetTimerClock()
{
  timespec ts;

  clock_gettime(CLOCK_MONOTONIC_RAW, &ts);

  // Here are some other clocks, but they are not very useful
  // clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  // clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);

  return ((unsigned long)ts.tv_sec)*1000000UL + ((unsigned long)ts.tv_nsec/1000); 
}

const unsigned long CLOCK_SCALE = 1000000UL;

#else
#warning "using low-resolution clock"

// NOTE: the clock used by clock_gettime seems to be much higher
// resolution than that used by clock.

unsigned long GetTimerClock()
{
  return clock();
}

const unsigned long CLOCK_SCALE = (unsigned long) CLOCKS_PER_SEC;
#endif



bool timer_compare(const FHEtimer *a, const FHEtimer *b)
{
  return strcmp(a->name, b->name) < 0;
}



static vector<FHEtimer *> timerMap;
static FHE_MUTEX_TYPE timerMapMx;

void registerTimer(FHEtimer *timer)
{
  FHE_MUTEX_GUARD(timerMapMx);
  timerMap.push_back(timer);
}

// Reset a timer for some label to zero
void FHEtimer::reset()
{
  numCalls = 0;
  counter = 0;
}


// Read the value of a timer (in seconds)
double FHEtimer::getTime() const // returns time in seconds
{
  // If the counter is currently counting, add the clock() value
  return ((double)counter)/CLOCK_SCALE;
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

const FHEtimer *getTimerByName(const char *name)
{
  for (long i = 0; i < long(timerMap.size()); i++) {
    if (strcmp(name, timerMap[i]->name) == 0)
      return timerMap[i];
  }

  return 0;
}

bool printNamedTimer(ostream& str, const char* name)
{
  for (long i = 0; i < long(timerMap.size()); i++) {
    if (strcmp(name, timerMap[i]->name) == 0) {
      
      long n = timerMap[i]->getNumCalls();
      if (n>0) {
        double t = timerMap[i]->getTime();
        double ave = t/n;
    
        str << "  " << name << ": " << t << " / " << n << " = " 
    	<< ave << "   [" << timerMap[i]->loc << "]\n";
      }
      else {
        str << "  " << name << " -- [" << timerMap[i]->loc << "]\n";
      }
      return true;
    }
  }
  return false;
}
