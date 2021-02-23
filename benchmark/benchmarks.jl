using BenchmarkTools
using BenchmarkTools: tune!, allocs,
      prettytime, time, prettymemory, memory, generate_benchmark_definition,
      Parameters, benchmarkable_parts, collectvars, warmup, Trial

include("utils.jl")

SUITE = BenchmarkGroup()

include("points.jl")

t = @elapsed using Meshes
println("Julia version is $VERSION")
println(string("Meshes.jl loading time: \e[33;1;1m$t\e[m seconds"))
println()
println("Benchmarking Meshes.jl...")
println()

tune!(SUITE)
result = run(SUITE)
display_group(result)
