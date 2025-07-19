# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using DataStructures: BinaryHeap
using LinearAlgebra: norm
import IterTools

struct AdaptiveSampling{TV} <: ContinuousSamplingMethod
    tol::TV
    errfun::Function
    min_points::Int
    max_points::Int
    function AdaptiveSampling(tol::TV, errfun, min_points, max_points)  where {TV}
        # assertion(min_points >= 2, "Need at least two initial points")
        new{TV}(tol, errfun, min_points, max_points)
    end
end

AdaptiveSampling(tol=5e-4, min_points=20.0, max_points=4_000.0) = AdaptiveSampling(tol, err_relative_range, min_points, max_points)

# SOMEDAY should we make a module to isolate these helpers from the rest of the package namespace?

# --------------
# Helper structs
# --------------

# using Test

"A segment of three `(t, v)` values where `t` is the parameter value and `v` is the point. Require `t1 < t_mid < t2` but not checked."
# SOMEDAY is there a better struct for this somewhere already?
struct ParametrizedSegment{T,V}
    t1::T
    v1::V
    t_mid::T
    v_mid::V
    t2::T
    v2::V
end

function ParametrizedSegment(f, t1, v1, t2, v2)
    t_mid = 0.5 * (t1 + t2)
    v_mid = f(t_mid)
    ParametrizedSegment(t1, v1, t_mid, v_mid, t2, v2)
end

function ParametrizedSegment(f, t1, t2)
    ParametrizedSegment(f, t1, f(t1), t2, f(t2))
end

"A range of vectors, from the 'lower-left' corner to the 'upper-right' one."
# TODO replace by Box
mutable struct Range{V}
    # Implicit: All lengths are equal
    mins::Vector{V}
    maxs::Vector{V}
end

"Iterator for the widths of each dimension. An iterator of Float64."
widths(ra::Range) = (ma-mi for (ma, mi) in zip(ra.maxs, ra.mins))

function Range(vs)
    vs_vectors = map(to, vs)
    ixs = eachindex(first(vs_vectors)) |> collect
    mins = [minimum(v[i] for v in vs_vectors) for i in ixs]
    maxs = [maximum(v[i] for v in vs_vectors) for i in ixs]
    Range(mins, maxs)
end

function empty_range(n::Integer)
    Range(fill(0.0, n), fill(0.0, n))
end

"All lengths must be the same and v must use linear indexing"
function push_range!(range::Range, v)
    v_vector = to(v)
    for i in eachindex(v_vector)
        range.mins[i] = min(range.mins[i], v_vector[i])
        range.maxs[i] = max(range.maxs[i], v_vector[i])
    end
end

"Check error value against tolerance and push if above."
function maybe_push_errfun!(queue, errfun, tol, x, args...)
    err = errfun(x, args...)
    if err > tol
        push!(queue, (-err, x))
    end
end

"An `errfun` for the sampling functions. l-infinity error relative to the range, component-wise."
function err_relative_range(s::ParametrizedSegment, ra::Range)
    # unitful handling; probably not the most efficient way.
    unit_ = unit(ra.mins[1])
    # eps() is for handling of zero widths.
    # SOMEDAY is this good enough or actually still kinda unstable?
    widths1 = max.(widths(ra), eps() * unit_)

    errs = @. (to(s.v_mid) - 0.5 * (to(s.v1) + to(s.v2))) / widths1
    norm(errs, Inf)
end

# @test err_relative_range(
#     Segment(t -> (t, t^2), 0.0, 10.0),
#     Range([0.0, 0.0], [1000.0, 1000.0]),
# ) == 0.025

# --------
# sampling
# --------

function sample(::AbstractRNG, geom::Geometry, method::AdaptiveSampling)
    assertion(paramdim(geom) == 1, "Not implemented: Adaptive sampling for > 1 parameter dimension")

    T = numtype(lentype(geom))
    t_min = zero(T)
    t_max = one(T)

    ts_init = range(t_min, t_max, length = method.min_points)
    vs_init = geom.(ts_init)
    
    # We store initial points here and append to it as we go (out of order, sorted later)
    # SOMEDAY smarter accounting to avoid the final sort?
    # We could use the Mesh structure here to insert points as we go, also towards a multi-dim generalization.
    ps = zip(ts_init, vs_init) |> collect

    # Initialize range. We'll update as we go.
    ra = Range(vs_init)

    # Queue of pending splits, in error order, so we can stick to `max_points`. We sort again later.
    # Initialize with pairs of successive points.
    queue = BinaryHeap{Tuple{Float64,ParametrizedSegment}}(Base.By(first), [])
    for ((t1, v1), (t2, v2)) in IterTools.partition(ps, 2, 1)
        s = ParametrizedSegment(geom, t1, v1, t2, v2)
        maybe_push_errfun!(queue, method.errfun, method.tol, s, ra)
    end

    # Refinement.
    while !isempty(queue) && length(ps) < method.max_points
        (_, s) = pop!(queue)
        push!(ps, (s.t_mid, s.v_mid))
        s1 = ParametrizedSegment(geom, s.t1, s.v1, s.t_mid, s.v_mid)
        s2 = ParametrizedSegment(geom, s.t_mid, s.v_mid, s.t2, s.v2)
        push_range!(ra, s1.v_mid)
        push_range!(ra, s2.v_mid)

        maybe_push_errfun!(queue, method.errfun, method.tol, s1, ra)
        maybe_push_errfun!(queue, method.errfun, method.tol, s2, ra)
    end

    if !isempty(queue)
        max_err = first(queue)[1]
        @warn "max_points=$(method.max_points) reached with errors above tolerance. Max error = $(-max_err) > $(method.tol) = tolerance"
    end

    sort!(ps; by = first)
    return [p[2] for p in ps]
end

