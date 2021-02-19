using BenchmarkTools
using BenchmarkTools: tune!, allocs,
      prettytime, time, prettymemory, memory, generate_benchmark_definition,
      Parameters, benchmarkable_parts, collectvars, warmup, Trial

"""
Non-macro adaptation from `BenchmarkTools.@benchmarkable`.
"""
function benchmarkable(ex, params = [])
  core, setup, teardown, _ = benchmarkable_parts([:(Ref($ex)[])])

  # extract any variable bindings shared between the core and setup expressions
  setup_vars = isa(setup, Expr) ? collectvars(setup) : []
  core_vars = isa(core, Expr) ? collectvars(core) : []
  out_vars = filter(in(setup_vars), core_vars)

  # generate the benchmark definition
  bench = generate_benchmark_definition(
    @__MODULE__,
    out_vars,
    setup_vars,
    core,
    setup,
    teardown,
    Parameters(; params...),
  )
end

function add_benchmark!(group::BenchmarkGroup, ex::Expr)
  group[string(ex)] = benchmarkable(ex)
end

function add_benchmark!(group::BenchmarkGroup, exs)
  add_benchmark!.(Ref(group), exs)
  group
end

function display_trial(trial::Trial)
  trialallocs = allocs(trial)
  trialmin = minimum(trial)

  println(
    "\e[33;1;1m",
    prettytime(time(trial)),
    "\e[m (", trialallocs , " allocation",
    trialallocs == 1 ? "" : "s", ": ",
    prettymemory(memory(trial)), ")"
  )
end

function display_group(group::BenchmarkGroup)
    nchars_left = maximum(length, keys(group.data))
    foreach(group.data) do (str, val)
      if val isa Trial
        print(' '^2, rpad(str, nchars_left), " :  ")
        display_trial(val)
      else
        println("\e[34;1;1m", str, "\e[m")
        display_group(val)
      end
    end
end

suite = BenchmarkGroup()

suite["points"] = BenchmarkGroup(["microbenchmark", "points"]) # add tags for eventual filtering

add_benchmark!(suite["points"], [
  :(Point(0., 1.)),
  :(Point([0., 1.])),
  :(Point([0., 1])),
  :(Point(0., convert(Float64, 1))),
])

println("Julia version is $VERSION")
t = @elapsed using Meshes
println(string("Meshes.jl loading time: \e[33;1;1m$t\e[m seconds"))
println()
println("Benchmarking Meshes.jl...")
println()

tune!(suite)
result = run(suite)

display_group(result)
