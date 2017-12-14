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
/**
 * @file multicore.h
 * @brief Support for multi-threaded implementations
 **/

#ifndef FHE_multicore_H
#define FHE_multicore_H

#ifdef FHE_THREADS

#include <atomic>
#include <mutex>





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
