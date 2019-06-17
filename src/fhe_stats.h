#ifndef FHE_STATS_H
#define FHE_STATS_H

#include <vector>
#include <iostream>

struct fhe_stats_record {
   const char *name;
   long count;
   double sum;
   double max;

   static std::vector<fhe_stats_record*> map;

   fhe_stats_record(const char *_name);
   void update(double val);
};


#define FHE_STATS_UPDATE(name, val) \
  do { \
     if (fhe_stats) { \
       static fhe_stats_record _local_stats_record(name); \
       _local_stats_record.update(val);  \
     } \
  } while (0) 


void
print_stats(std::ostream& s);

extern bool fhe_stats;


#endif
