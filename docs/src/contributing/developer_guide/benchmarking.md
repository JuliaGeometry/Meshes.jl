# Benchmarking

## Overview

We use [PkgBenchmark](https://github.com/JuliaCI/PkgBenchmark.jl) to conduct benchmarks, coupled with [BenchmarkCI](https://github.com/tkf/BenchmarkCI.jl) to automate benchmark runs with [GitHub Actions](https://github.com/features/actions) on certain PRs. Currently, PRs will _not_ run benchmarks by default. Only PRs with the label `run benchmark` will trigger benchmarks on push, adding a comment with the benchmark results.

## Run benchmarks locally

#### With PkgBenchmark

To locally run benchmarks, you can use PkgBenchmark with the following code:

```julia
using Meshes, PkgBenchmark # make sure you have PkgBenchmark installed, e.g. globally
benchmarkpkg(Meshes)
```

It will include `benchmark/benchmarks.jl` and look for a [`BenchmarkGroup`](https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/doc/manual.md#the-benchmarkgroup-type) variable named `SUITE`, that it will run for you with a nice printing.

The full list of options is detailed in the [PkgBenchmark documentation](https://juliaci.github.io/PkgBenchmark.jl/stable/run_benchmarks/#PkgBenchmark.benchmarkpkg).

#### Manually

Sometimes you may prefer to run the suite manually, especially during interactive development. All you have to do is include the `benchmark/benchmarks.jl` file, and run the `BenchmarkGroup` suite:

```julia
include("benchmark/benchmarks.jl")
run(SUITE)
```

## Limitations

Note that because GitHub Actions may use different runners between benchmarks, you are likely to see fluctuations and performance changes that may not always be relevant. This makes it somewhat unreliable for tracking regressions.

Benchmark results can be very useful for validation, for example when optimizing or modifying existing features, but you are encouraged to benchmark important changes locally.
