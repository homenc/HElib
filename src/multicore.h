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
 * @file multicore.h
 * @brief Support for multi-threaded implementations
 **/

#ifndef FHE_multicore_H
#define FHE_multicore_H

#ifdef FHE_THREADS

#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <NTL/SmartPtr.h>

using namespace std;
using namespace NTL;


/*************************************************************

Some simple thread pooling.

You create a thread pool by constructing a MultiTask object.
For example:

   long nthreads = 4;
   MultiTask multiTask(nthreads);

creates a thread pool of 4 threads.  These threads will exist
until the destructor for multiTask is called.  

The simplest way to use a thread pools is as follows.
Suppose you have a task that consists of N subtasks,
indexed 0..N-1.  Then you can write:


   multiTask.exec1(N, 
      [&](long first, long last) {
         for (long i = first; i < last; i++) {
            ... code to process subtask i ...
         }
      }
   );

The second argument to exec1 is a C++11 "lambda".
The "[&]" indicates that all local variables in the calling
context are captured by reference, so the lambda body can 
reference all visible local variables directly.

A lower-level interface is also provided.
One can write:

   multiTask.exec(n,
      [&](long index) {
         ... code to process index i ...
      }
   );

This will activate n threads with indices 0..n-1, and execute
the given code on each index.  The parameter n must be
in the range 1..nthreads, otherwise an error is raised.

This lower-level interface is useful in some cases,
especially when memory is managed in some special way.
For convenience, a method is provided to break
subtasks up into smaller, almost-equal-sized groups
of subtasks:

   Vec<long> pvec;
   long n = multiTask.SplitProblems(N, pvec);

can be used for this.  N is the number of subtasks, indexed 0..N-1.
This method will compute n as needed by exec, and 
the range of subtasks to be processed by a given index in the range
0..n-1 is pvec[index]..pvec[index+1]-1
Thus, the logic of the above exec1 example can be written
using the lower-level exec interface as follows:

   
   Vec<long> pvec;
   long n = multiTask.SplitProblems(N, pvec);
   multiTask.exec(n,
      [&](long index) {
         long first = pvec[index];
         long last = pvec[index+1];
         for (long i = first; i < last; i++) {
            ... code to process subtask i ...
         }
      }
   );

However, with this approach, memory or other resources can be
assigned to each index = 0..n-1, and managed externally. 


TODO: think about how to manage exceptions thrown in one 
of the threads in the pool.  Presumably, before launch,
we set a "dirty bit" to false...if any thread throws an exception,
it is caught, and the "dirty bit" is set to true.
After we wait for all threads to complete, we check the
"dirty bit", and throw an exception.
Unfortunately, it is not clear how to transmit the exception...
and there could be several...it looks like in C++11, std::exception_ptr
is the tool to use, in conjunction with std::current_exception
and std::rethrow_exception.  So I could remember the first exception
thrown, and that is all.  I would need to mutex-protect storing this
exception info, but performance is not critical.



*************************************************************/


template<class T>
class SimpleSignal {
private:
  T val; 
  mutex m;
  condition_variable cv;

  SimpleSignal(const SimpleSignal&); // disabled
  void operator=(const SimpleSignal&); // disabled

public:
  SimpleSignal() : val(0) { }

  T wait() 
  {
    unique_lock<mutex> lock(m);
    cv.wait(lock, [&]() { return val; } );
    T old_val = val;
    val = 0;
    return old_val;
  }

  void send(T new_val)
  {
    lock_guard<mutex> lock(m);
    val = new_val;
    cv.notify_one();
  }
};



class ConcurrentTask {
public:
  virtual void run(long index) = 0;
};


// dummy class, used for signalling termination
class ConcurrentTaskTerminate : public ConcurrentTask {
public:
  void run(long index) { }
};

template<class Fct>
class ConcurrentTaskFct : public ConcurrentTask {
public:
  Fct fct;
  ConcurrentTaskFct(Fct&& _fct) : fct(std::move(_fct)) { }

  void run(long index) { fct(index); }
};

template<class Fct>
class ConcurrentTaskFct1 : public ConcurrentTask {
public:
  Fct fct;
  const Vec<long>& pvec;
  ConcurrentTaskFct1(Fct&& _fct, const Vec<long>& _pvec) : 
    fct(std::move(_fct)), pvec(_pvec)  { }

  void run(long index) { fct(pvec[index], pvec[index+1]); }
};


struct AutomaticThread;
struct MultiTask;
inline void doMultiTask(AutomaticThread *autoThread);

struct AutomaticThread {
   MultiTask *multiTask;
   long index;
   SimpleSignal< ConcurrentTask * > localSignal;
   thread t;

   ConcurrentTaskTerminate termTask;
   // dummy task used to signal termination

   AutomaticThread(MultiTask *_multiTask, long _index) 
     : multiTask(_multiTask), index(_index),
       t(doMultiTask, this) 
   { 
      // cerr << "starting thread " << t.get_id() << "\n";
   }

   inline ~AutomaticThread();
};

struct MultiTask {
  long nthreads;


  atomic<long> counter;
  SimpleSignal<bool> globalSignal;

  Vec< UniquePtr<AutomaticThread> > threadVec;

  MultiTask(const MultiTask&); // disabled
  void operator=(const MultiTask&); // disabled

  MultiTask(long _nthreads) : nthreads(_nthreads), counter(0)
  {
    assert(nthreads > 0);

    threadVec.SetLength(nthreads);

    for (long i = 0; i < nthreads; i++) {
      threadVec[i].make(this, i);
    }
  }

  void begin(long cnt)
  {
    assert(cnt > 0 && cnt <= nthreads);

    counter = cnt;
  }

  void end()
  {
    globalSignal.wait();
  }

  void launch(ConcurrentTask *task, long index)
  {
    assert(task != 0 && index >= 0 && index < nthreads);

    threadVec[index]->localSignal.send(task);
  }

  

  // High level interfaces, intended to be used with lambdas

  // In this version, fct takes one argument, which is
  // an index in [0..cnt)

  template<class Fct>
  void exec(long cnt, Fct fct) 
  {
    ConcurrentTaskFct<Fct> task(std::move(fct));

    begin(cnt);
    for (long t = 0; t < cnt; t++) launch(&task, t);
    end();
  }

  // even higher level version: sz is the number of subproblems,
  // and fct takes two args, first and last, so that subproblems
  // [first..last) are processed.

  template<class Fct>
  void exec1(long sz, Fct fct) 
  {
    Vec<long> pvec;
    long cnt = SplitProblems(sz, pvec);
    ConcurrentTaskFct1<Fct> task(std::move(fct), pvec);

    begin(cnt);
    for (long t = 0; t < cnt; t++) launch(&task, t);
    end();
  }


  // splits nproblems problems among (at most) nthreads threads.
  // returns the actual number of threads nt to be used, and 
  // initializes pvec to have length nt+1, so that for t = 0..nt-1,
  // thread t processes subproblems pvec[t]..pvec[t+1]-1
  long SplitProblems(long nproblems, Vec<long>& pvec) const
  {
    long blocksz = (nproblems + nthreads - 1)/nthreads;
    long nt = (nproblems + blocksz - 1)/blocksz;
  
    pvec.SetLength(nt+1);
  
    for (long t = 0; t < nt; t++) pvec[t] = blocksz*t;
    pvec[nt] = nproblems;
  
    return nt;
  }


  long getNumThreads() const { return nthreads; }


};


inline AutomaticThread::~AutomaticThread()
{
  // cerr << "stopping thread " << t.get_id() << "...";
  localSignal.send(&termTask);
  t.join();
  // cerr << "\n";
}


inline void doMultiTask(AutomaticThread *autoThread)
{
  MultiTask *multiTask = autoThread->multiTask;
  long index = autoThread->index;

  for (;;) {
    ConcurrentTask *task = autoThread->localSignal.wait();

    if (task == &autoThread->termTask) return;
    // special test for termination

    task->run(index);
    if (--(multiTask->counter) == 0) multiTask->globalSignal.send(true);
  }
}

#define FHE_atomic_long atomic_long
#define FHE_atomic_ulong atomic_ulong

#define FHE_MUTEX_TYPE mutex
#define FHE_MUTEX_GUARD(mx) lock_guard<mutex> _lock ## __LINE__ (mx)

#else

#define FHE_atomic_long long
#define FHE_atomic_ulong unsigned long

#define FHE_MUTEX_TYPE int
#define FHE_MUTEX_GUARD(mx) ((void) mx)

#endif



#endif
