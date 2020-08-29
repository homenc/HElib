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

#ifndef HELIB_STATS_H
#define HELIB_STATS_H

#include <vector>
#include <iostream>

namespace helib {

struct fhe_stats_record
{
  const char* name;
  long count;
  double sum;
  double max;

  std::vector<double> saved_values;
  // save all values --- only used if explicitly requested

  static std::vector<fhe_stats_record*> map;

  fhe_stats_record(const char* _name);
  void update(double val);
  void save(double val);
};

#define HELIB_STATS_UPDATE(name, val)                                          \
  do {                                                                         \
    if (fhe_stats) {                                                           \
      static fhe_stats_record _local_stats_record(name);                       \
      _local_stats_record.update(val);                                         \
    }                                                                          \
  } while (0)

#define HELIB_STATS_SAVE(name, val)                                            \
  do {                                                                         \
    if (fhe_stats) {                                                           \
      static fhe_stats_record _local_stats_record(name);                       \
      _local_stats_record.save(val);                                           \
    }                                                                          \
  } while (0)

void print_stats(std::ostream& s);

const std::vector<double>* fetch_saved_values(const char*);

extern bool fhe_stats;

} // namespace helib

#endif
