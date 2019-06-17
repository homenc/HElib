
#include "fhe_stats.h"
#include "multicore.h"
#include <algorithm>
#include <utility>
#include <cstring>

bool fhe_stats = false;

static std::vector<fhe_stats_record*> stats_map;
static FHE_MUTEX_TYPE stats_mutex;

fhe_stats_record::fhe_stats_record(const char *_name)
  : name(_name), count(0), sum(0), max(0) 
// FIXME: maybe max should be -infty?
{
  FHE_MUTEX_GUARD(stats_mutex);
  stats_map.push_back(this);
}

void
fhe_stats_record::update(double val)
{
  FHE_MUTEX_GUARD(stats_mutex);
  count++;
  sum += val;
  if (val > max) max = val;
}


static bool 
stats_compare(const fhe_stats_record *a, const fhe_stats_record *b)
{
  return strcmp(a->name, b->name) < 0;
}

void 
print_stats(std::ostream& s)
{
  s << "||||| stats |||||\n";
  sort(stats_map.begin(), stats_map.end(), stats_compare);
  for (long i = 0; i < long(stats_map.size()); i++) {
    const char *name = stats_map[i]->name;
    long count = stats_map[i]->count;
    double sum = stats_map[i]->sum;
    double max = stats_map[i]->max;

    if (count > 0) {
      s << name << " ave=" << (sum/count) << " max=" << max << "\n";
    }
  }
}
