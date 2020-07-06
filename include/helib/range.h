/* Copyright (C) 2012-2020 IBM Corp.
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

#ifndef HELIB_RANGE_H
#define HELIB_RANGE_H

namespace helib {

template <typename T>
class general_range
{
public:
  class iterator
  {
    friend class general_range;

  public:
    T operator*() const { return i_; }
    iterator& operator++()
    {
      ++i_;
      return *this;
    }

    bool operator==(const iterator& other) const { return i_ == other.i_; }
    bool operator!=(const iterator& other) const { return i_ != other.i_; }

  protected:
    iterator(T start) : i_(start) {}

  private:
    T i_;
  };

  iterator begin() const { return begin_; }
  iterator end() const { return end_; }
  general_range(T begin, T end) : begin_(begin), end_(end)
  {
    if (end < begin)
      end = begin;
  }

private:
  iterator begin_;
  iterator end_;
};

inline general_range<long> range(long n) { return general_range<long>(0, n); }

inline general_range<long> range(long m, long n)
{
  return general_range<long>(m, n);
}

} // namespace helib

#endif // ifndef HELIB_RANGE_H
