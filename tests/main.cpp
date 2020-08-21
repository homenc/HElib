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

#include "gtest/gtest.h"
#include "test_common.h"

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  // Now argc and argv will have had their gtest_* entries
  // processed and removed.  Parse for our own purposes.
  helib_test::parse_common_args(argc, argv);
  return RUN_ALL_TESTS();
}
