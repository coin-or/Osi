// Benchmark for OsiFeatures extraction.
// Usage: osi_features_bench <instance.mps[.gz]>
//
// Output (tab-separated, one line, easy to parse):
//   instance_name\tread_time\tfeatures_time\tfeature1_name\tfeature1_val\t...

#include "OsiClpSolverInterface.hpp"
#include "OsiFeatures.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>

int main(int argc, char *argv[])
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <instance.mps[.gz]>\n", argv[0]);
    return 1;
  }

  const char *path = argv[1];

  // Extract instance name from path (strip directory and extensions)
  const char *base = strrchr(path, '/');
  base = base ? base + 1 : path;
  char name[512];
  snprintf(name, sizeof(name), "%s", base);
  // Strip .gz then .mps/.lp
  char *dot = strrchr(name, '.');
  if (dot && strcmp(dot, ".gz") == 0) {
    *dot = '\0';
    dot = strrchr(name, '.');
  }
  if (dot)
    *dot = '\0';

  using Clock = std::chrono::high_resolution_clock;

  // Read
  OsiClpSolverInterface si;
  si.getModelPtr()->setLogLevel(0);

  auto t0 = Clock::now();
  int rc = si.readMps(path);
  auto t1 = Clock::now();

  if (rc) {
    fprintf(stderr, "ERROR: failed to read %s (rc=%d)\n", path, rc);
    return 1;
  }

  double readTime = std::chrono::duration<double>(t1 - t0).count();

  // Extract features
  double f[OFCount];

  auto t2 = Clock::now();
  OsiFeatures::compute(f, &si);
  auto t3 = Clock::now();

  double featTime = std::chrono::duration<double>(t3 - t2).count();

  // Print header + values (tab-separated)
  printf("instance\tread_sec\tfeatures_sec");
  for (int i = 0; i < OFCount; ++i)
    printf("\t%s", OsiFeatures::name(i));
  printf("\n");

  printf("%s\t%.6f\t%.6f", name, readTime, featTime);
  for (int i = 0; i < OFCount; ++i)
    printf("\t%g", f[i]);
  printf("\n");

  return 0;
}
