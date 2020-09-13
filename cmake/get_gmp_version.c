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

#include <stdio.h>
#include <gmp.h>

int main(void)
{
  const unsigned gmp_version = __GNU_MP_VERSION;
  const unsigned gmp_version_minor = __GNU_MP_VERSION_MINOR;
  const unsigned gmp_version_patchlevel = __GNU_MP_VERSION_PATCHLEVEL;
  printf("%d.%d.%d", gmp_version, gmp_version_minor, gmp_version_patchlevel);
  return 0;
}
