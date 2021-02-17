using BenchmarkTools: hasevals, tune!, run_result, minimum, allocs, prettytime, time, prettymemory, memory, generate_benchmark_definition, Parameters, benchmarkable_parts, collectvars, warmup

"""
Non-macro adaptation from `BenchmarkTools.@benchmarkable`
"""
function benchmarkable(ex, params)
  core, setup, teardown, _ = benchmarkable_parts([:(Ref($ex)[])])

  # extract any variable bindings shared between the core and setup expressions
  setup_vars = isa(setup, Expr) ? collectvars(setup) : []
  core_vars = isa(core, Expr) ? collectvars(core) : []
  out_vars = filter(var -> var in setup_vars, core_vars)

  # generate the benchmark definition
  bench = generate_benchmark_definition(@__MODULE__,
                                  out_vars,
                                  setup_vars,
                                  core,
                                  setup,
                                  teardown,
                                  Parameters(; params...))
end

"""
Non-macro Adaptation from `BenchmarkTools.@btime`.
"""
function btime(ex; params...)
  bench = benchmarkable(ex, params)
  warmup(bench)
  !hasevals(params) && tune!(bench)
  trial, result = run_result(bench)
  trialmin = minimum(trial)
  trialallocs = allocs(trialmin)
  string("  \e[33;1;1m",
          prettytime(time(trialmin)),
          "\e[m (", trialallocs , " allocation",
          trialallocs == 1 ? "" : "s", ": ",
          prettymemory(memory(trialmin)), ")")
end

function run_benchmarks(exs, message = nothing)
  !isnothing(message) && println("\e[34;1;1m", message, "\e[m")
  nchars_left = maximum(length, string.(exs))
  map(exs) do ex
    println(' '^2, rpad(ex, nchars_left), " :", btime(ex))
  end
end

println("Julia version is $VERSION")
t = @elapsed using Meshes
println(string("Meshes.jl loading time: \e[33;1;1m$t\e[m seconds"))

const point_exs = [
  :(Point(0., 1.)),
  :(Point([0., 1.])),
  :(Point([0., 1])),
  :(Point(0., convert(Float64, 1))),
]

println()
println("Benchmarking Meshes.jl...")
println()
run_benchmarks(point_exs, "Points")
