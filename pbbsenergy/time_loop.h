#include "../parlay/internal/get_time.h"

#include <rapl_energy.h>

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

    struct RaplEnergy *rapl;
    rapl = rapl_start();

    t.start();
    runf();
    t.next("");

    struct RaplElapsed *elapsed;
    elapsed = rapl_elapsed(rapl);

    for (uintptr_t i = 0; i < elapsed->len; i++) {
        printf(" %lf", elapsed->energy[i]);
    }

    elapsed_free(elapsed);
    rapl_free(rapl);

    endf();
  }
}
