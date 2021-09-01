/* Copyright (C) 2021 Intel Corporation
 * SPDX-License-Identifier: Apache-2.0
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

// This program performs basic operations n times and is used for profiling.

#include <iostream>

#include <helib/helib.h>
#include <helib/EncryptedArray.h>
#include <helib/ArgMap.h>
#include <NTL/BasicThreadPool.h>

void add(helib::Ctxt& ctxt, int iter)
{
  for (int i = 0; i < iter; ++i) {
    // Create tmp to avoid noise growth.
    helib::Ctxt tmp(ctxt);
    tmp += ctxt;
  }
}

void mult(helib::Ctxt& ctxt, int iter)
{
  for (int i = 0; i < iter; ++i) {
    // Create tmp to avoid noise growth.
    helib::Ctxt tmp(ctxt);
    tmp *= ctxt;
  }
}

void rotate(helib::Ctxt& ctxt, const helib::EncryptedArray& ea, int iter)
{
  for (int i = 0; i < iter; ++i) {
    // Create tmp to avoid noise growth.
    helib::Ctxt tmp(ctxt);
    ea.rotate(tmp, 1);
  }
}

void scalarAdd(helib::Ctxt ctxt, int iter)
{
  for (int i = 0; i < iter; ++i) {
    ctxt += 1l;
  }
}

void scalarMult(helib::Ctxt ctxt, int iter)
{
  for (int i = 0; i < iter; ++i) {
    ctxt *= 2l;
  }
}

int main(int argc, char* argv[])
{
  // Plaintext prime modulus
  unsigned long p = 131;
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m = 130; // this will give 48 slots
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of bits of the modulus chain
  unsigned long bits = 1000;
  // Number of columns of Key-Switching matrix (default = 2 or 3)
  unsigned long c = 2;
  // Size of NTL thread pool (default =1)
  unsigned long nthreads = 1;

  helib::ArgMap()
      .optional()
      .named()
      .arg("m", m, "Cyclotomic polynomial ring")
      .arg("p", p, "Plaintext prime modulus")
      .arg("r", r, "Hensel lifting")
      .arg("bits", bits, "# of bits in the modulus chain")
      .arg("c", c, "# fo columns of Key-Switching matrix")
      .arg("nthreads", nthreads, "Size of NTL thread pool")
      .parse(argc, argv);

  // set NTL Thread pool size
  if (nthreads > 1)
    NTL::SetNumThreads(nthreads);

  helib::Context context = helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .bits(bits)
                               .c(c)
                               .build();

  helib::SecKey secret_key = helib::SecKey(context);
  secret_key.GenSecKey();
  helib::addSome1DMatrices(secret_key);
  const helib::PubKey& public_key = secret_key;
  const helib::EncryptedArray& ea = context.getEA();

  std::cout << std::endl;
  context.printout();

  long nslots = ea.size();

  helib::PtxtArray pa(context);
  pa.random();
  helib::Ctxt ctxt(public_key);
  pa.encrypt(ctxt);

  int opt = 0;
  std::cout << "What operation(s) do you want to run? [0:all, 1:add, 2:mult, "
               "3:rotate, 4:scalarAdd, 5:scalarMult]: ";
  std::cin >> opt;

  int iter = 0;
  std::cout << "How many iterations?: ";
  std::cin >> iter;
  iter = (iter < 1) ? 10 : iter;

  std::cout << "Option: " << opt << ", " << iter << " times.\n";

  switch (opt) {
  case 0:
    std::cout << "Running option " << opt << ": all operations...\n";
    add(ctxt, iter);
    mult(ctxt, iter);
    rotate(ctxt, ea, iter);
    scalarAdd(ctxt, iter);
    scalarMult(ctxt, iter);
    std::cout << "Done!\n";
    break;
  case 1:
    std::cout << "Running option " << opt << ": additions...\n";
    add(ctxt, iter);
    std::cout << "Done!\n";
    break;
  case 2:
    std::cout << "Running option " << opt << ": multiplications...\n";
    mult(ctxt, iter);
    std::cout << "Done!\n";
    break;
  case 3:
    std::cout << "Running option " << opt << ": rotations...\n";
    rotate(ctxt, ea, iter);
    std::cout << "Done!\n";
    break;
  case 4:
    std::cout << "Running option " << opt << ": scalar additions...\n";
    scalarAdd(ctxt, iter);
    std::cout << "Done!\n";
    break;
  case 5:
    std::cout << "Running option " << opt << ": scalar multiplications...\n";
    scalarMult(ctxt, iter);
    std::cout << "Done!\n";
    break;
  default:
    std::cout << "Received invalid option: " << opt << std::endl;
  }

  return 0;
}
