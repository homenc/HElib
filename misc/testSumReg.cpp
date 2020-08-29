/* Copyright (C) 2020 IBM Corp.
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
#include <iostream>
#include "../include/helib/SumRegister.h"

using namespace helib;

int main(int argc, char** argv)
{

  if (argc < 2) {
    std::cout << "Give integer\n";
    exit(1);
  }

  int howMany = atoi(argv[1]);
  SumRegister<int> sr(howMany);
  std::cout << "Depth: " << sr.getDepth() << std::endl;

  std::unique_ptr<int> num = std::unique_ptr<int>(new int);

  for (int i = 0; i < howMany; i++) {
    num.reset(new int(i + 1));
    std::cout << "Adding: " << *num << '\n';
    sr.add(num);
    sr.print();
  }

  num.reset(new int(0));
  sr.add(num);
  sr.print();

  if (sr.hasResult())
    std::cout << "Result:" << *sr.getResult() << std::endl;

  sr.print();

  return 0;
}
