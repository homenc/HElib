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

#include <atomic>
#include <mutex>

using namespace std;




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
