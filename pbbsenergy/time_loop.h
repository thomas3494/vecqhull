#include "../parlay/internal/get_time.h"

#include <energy_bench.h>

template<class F, class G, class H>
void time_loop(int rounds, double delay, F initf, G runf, H endf) {
  parlay::internal::timer t;
  // run for delay seconds to "warm things up"
  // will skip if delay is zero
  while (t.total_time() < delay) {
    initf(); runf(); endf();
  }
  for (int i=0; i < rounds; i++) {
    initf();

    EnergyInfo *start = start_energy_measure();

    t.start();
    runf();
    t.next("");

    EnergyResult *stop = stop_energy_measure(start);
    print_energy_results(&stop);
    free_energy_results(stop);

    endf();
  }
}
