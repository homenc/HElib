#ifndef FHE_multicore_H
#define FHE_multicore_H

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

};

inline void doMultiTask(MultiTask *multiTask, long index)
{
  for (;;) {
    ConcurrentTask *task = multiTask->localSignal[index].wait();
    task->run(index);
    if (--(multiTask->counter) == 0) multiTask->globalSignal.send(true);
  }
}



#endif
