using BenchmarkTools
t = @elapsed using Meshes
@info "Loading time: " * string(t) * "s"

# using SnoopCompile
# using MethodAnalysis

# meths = []
# visit(Meshes) do item
#     item isa Method && push!(meths, item)
#     true
# end

p = Polygon(p_outer, [])

@btime precompile(Polygon, (Vector{Tuple{Float64, Float64}}, Vector{Any}))