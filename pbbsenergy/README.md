
PBBS does not measure energy consumption.

We can hack it in by replacing `time_loop.h` in PBBS with the file
provided in this folder, which includes RAPL energy measuring.

```bash
cp pbbsenergy/time_loop.h pbbsbench/benchmarks/convexHull/bench/common/time_loop.h
cp pbbsenergy/parallelDefs pbbsbench/common/parallelDefs
```
