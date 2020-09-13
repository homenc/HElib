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
/**
 * @file multicore.h
 * @brief Support for multi-threaded implementations
 **/

#ifndef HELIB_MULTICORE_H
#define HELIB_MULTICORE_H

#ifdef HELIB_THREADS

#include <atomic>
#include <mutex>

namespace helib {

#define HELIB_atomic_long std::atomic_long
#define HELIB_atomic_ulong std::atomic_ulong

#define HELIB_MUTEX_TYPE std::mutex
#define HELIB_MUTEX_GUARD(mx) std::lock_guard<std::mutex> _lock##__LINE__(mx)

} // namespace helib

#else

namespace helib {

#define HELIB_atomic_long long
#define HELIB_atomic_ulong unsigned long

#define HELIB_MUTEX_TYPE int
#define HELIB_MUTEX_GUARD(mx) ((void)mx)

} // namespace helib

#endif // ifdef HELIB_THREADS

#endif // ifndef HELIB_MULTICORE_H
