/* Copyright (C) 2019-2020 IBM Corp.
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

#include <helib/fhe_stats.h>
#include <helib/multicore.h>
#include <algorithm>
#include <utility>
#include <cstring>

namespace helib {

bool fhe_stats = false;

static std::vector<fhe_stats_record*> stats_map;
static HELIB_MUTEX_TYPE stats_mutex;

fhe_stats_record::fhe_stats_record(const char* _name) :
    name(_name), count(0), sum(0), max(0)
// FIXME: maybe max should be -infty?
{
  HELIB_MUTEX_GUARD(stats_mutex);
  stats_map.push_back(this);
}

void fhe_stats_record::update(double val)
{
  HELIB_MUTEX_GUARD(stats_mutex);
  count++;
  sum += val;
  if (val > max)
    max = val;
}

void fhe_stats_record::save(double val)
{
  HELIB_MUTEX_GUARD(stats_mutex);
  saved_values.push_back(val);
}

static bool stats_compare(const fhe_stats_record* a, const fhe_stats_record* b)
{
  return strcmp(a->name, b->name) < 0;
}

void print_stats(std::ostream& s)
{
  s << "||||| stats |||||\n";
  sort(stats_map.begin(), stats_map.end(), stats_compare);
  for (long i = 0; i < long(stats_map.size()); i++) {
    const char* name = stats_map[i]->name;
    long count = stats_map[i]->count;
    double sum = stats_map[i]->sum;
    double max = stats_map[i]->max;

    if (count > 0) {
      s << name << " ave=" << (sum / count) << " max=" << max << "\n";
    }
  }
}

const std::vector<double>* fetch_saved_values(const char* name)
{
  for (long i = 0; i < long(stats_map.size()); i++) {
    if (strcmp(name, stats_map[i]->name) == 0)
      return &stats_map[i]->saved_values;
  }

  return 0;
}

} // namespace helib
