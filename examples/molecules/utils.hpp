////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2016
//
// Example using the basic_mesh class to genererate solvent excluded
// surfaces for molecules.
//


#undef min

inline std::ostream & operator<<(std::ostream &os, const glm::vec3 &v) {
  return os << "vec3(" << v.x << ", " << v.y << ", " << v.z << ")";
}

template <class... P>
const char *fmt(const char *str, P... params) {
  static char buf[256];
  snprintf(buf, sizeof(buf), str, params...);
  return buf;
}

template <class F>
void par_for(int begin, int end, F fn) {
  std::atomic<int> idx;
  idx = begin;
  int num_cpus = std::thread::hardware_concurrency();
  std::vector<std::future<void>> futures(num_cpus);
  for (int cpu = 0; cpu != num_cpus; ++cpu) {
    futures[cpu] = std::async(std::launch::async, [cpu, &idx, end, &fn]() {
      for (;;) {
        int i = idx++;
        printf("[%d %d]", i, cpu);
        fflush(stdout);
        if (i >= end) break;
        fn(i);
      }
    });
  }
  for (int cpu = 0; cpu != num_cpus; ++cpu) {
    futures[cpu].get();
  }
  printf("\n");
}

