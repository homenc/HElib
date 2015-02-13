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

template<class Fct>
class ConcurrentTaskFct : public ConcurrentTask {
public:
  Fct fct;
  ConcurrentTaskFct(Fct&& _fct) : fct(std::move(_fct)) { }

  void run(long index) { fct(index); }
};


struct MultiTask;
inline void doMultiTask(MultiTask *multiTask, long index);

struct MultiTask {
  long nthreads;

  atomic<long> counter;
  SimpleSignal<bool> globalSignal;
  Vec< SimpleSignal< ConcurrentTask * > > localSignal;

  MultiTask(const MultiTask&); // disabled
  void operator=(const MultiTask&); // disabled

  MultiTask(long _nthreads) : nthreads(_nthreads), counter(0)
  {
    assert(nthreads > 0);

    localSignal.SetLength(nthreads);

    for (long i = 0; i < nthreads; i++) {
      thread t(doMultiTask, this, i);
      t.detach();

      // I may want to think about leaving these threads
      // attached and joining with them in the destructor
      // somehow
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

    localSignal[index].send(task);
  }

  

  // High level interface, intended to be used with lambdas
  template<class Fct>
  void exec(long cnt, Fct fct) 
  {
    ConcurrentTaskFct<Fct> task(std::move(fct));

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

inline void doMultiTask(MultiTask *multiTask, long index)
{
  for (;;) {
    ConcurrentTask *task = multiTask->localSignal[index].wait();
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
