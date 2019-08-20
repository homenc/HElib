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
