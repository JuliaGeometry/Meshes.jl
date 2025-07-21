# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using DataStructures: BinaryHeap
using LinearAlgebra: norm
import IterTools

"""
    AdaptiveSampling(; [options])

Incrementally sample points such that the error from linear interpolation between these points is small. This is useful for plotting parametric functions or Bezier curves where curvature is very non-uniform.

This is *only* implemented for geometries with a one-dimensional parameter domain (i.e., realistically, [`ParametrizedCurve`](@ref) or [`BezierCurve`](@ref))!

## Options

* `tol` - Tolerance for the linear interpolation error. The meaning depends on `errfun`.
* `errfun` - Function to compute the interpolation error. The default is [`errRelativeBBox`](@ref) = the l-infinity distance component-wise relative to (an estimate of) the range of the bounding box; this is a good default for plotting. If provided, it must satisfy `errfun(::ParametrizedSegment, ::MutableBox)::E`
* `minPoints` - Number of points sampled uniformly initially.
* `maxPoints` - Maximum number of points to sample. We do our best to, if this is hit, get rid of the worst approximation errors first.

### Operation

This function works as follows: first, it performs an initial presampling step. Then, it considers the segment where the error of linear interpolation across the segment vs the midpoint value (midpoint in t space) is worst and splits it into two halves. Then repeat until all errors are below the tolerance or we run out of points. The estimate of the value range is updated as we add new points.

### Caveats

If `minPoints` is low, `range` in `errfun` may be catastrophically small, at least for initial refinement. You probably don't want that. More generally, if the initial samples do not represent the value range reasonably well, in some cases, refinement may add excessive points, and may even reach `maxPoints` without any meaningful progress.

Recursion is by splitting into halves in t space, so if your function has details that are missed by this recursive grid, they won't show up. In this case, increasing `minPoints` may help. 

Polar points are not handled well (see below).

### Features not yet implemented (SOMEDAY)

- Only geometries with one-dimensional parameter domains are supported. It would be cool to implement higher dimensions. We then probably want a discretization, not just sampling.
- Infinite or open t ranges where we would adaptively sample larger/smaller (towards infinity or the edges) t values. This may be an instance of a more general hook into refinement.
- Some way of detecting and avoiding polar points. Basically a way to say "it's not useful to sample this part, better to cut it out". Not clear how we'd do this in a robust way without affecting some use cases negatively.
- Optional penalty for recursion depth (or something like this) to avoid excessive concentration at polar points if `maxPoints` gets exhausted. This would have to be configured by the user with knowledge of the function, can be detrimental otherwise. Unclear if we want this.
- Option to split into more than two parts per recursion step, or to split not into halves but by some other proportions. This could help when details are missed by the 2-split recursion (see Caveats).
- Point density target in value space. Could help in situations where the midway point is interpolated well linearly, but there is additional variation at some other point.
- Option to drop points that are ultimately not needed if the range expands during refinement (may reduce number of points, might be useful for some later computation steps).
- A way to understand that the range may actually be the wrong thing, e.g., when we use
  `aspect_ratio=1` in the plot and the plot therefore creates additional space, or some xlim or
  ylim. Hard to do generically. Maybe make an option to pass the plot range explicitly.
- Is there some smartness we can do if the (second) derivative of f is available? Choosing the splitting point better than at the interval midpoint seems attractive.
"""
struct AdaptiveSampling{E} <: ContinuousSamplingMethod
    tol::E
    errfun::Function
    minPoints::Int
    maxPoints::Int
    function AdaptiveSampling(tol::E, errfun, minPoints, maxPoints) where {E}
        # assertion(minPoints >= 2, "Need at least two initial points")
        new{E}(tol, errfun, minPoints, maxPoints)
    end
end

AdaptiveSampling(;
    tol = 5e-4,
    errfun = errRelativeBBox,
    minPoints = 20,
    maxPoints = 4_000,
) = AdaptiveSampling(tol, errfun, minPoints, maxPoints)

# SOMEDAY should we make a module to isolate these helpers from the rest of the package namespace?

# --------------
# Helper structs
# --------------

# SOMEDAY do we want these inline tests?
# using Test

"A segment of three `(t, v)` values where `t` is the parameter value and `v` is the point. Require `t1 < tMid < t2` but not checked."
# SOMEDAY is there a better struct for this somewhere already?
struct ParametrizedSegment{T,V}
    t1::T
    v1::V
    tMid::T
    vMid::V
    t2::T
    v2::V
end

function ParametrizedSegment(f, t1, v1, t2, v2)
    tMid = 0.5 * (t1 + t2)
    vMid = f(tMid)
    ParametrizedSegment(t1, v1, tMid, vMid, t2, v2)
end

function ParametrizedSegment(f, t1, t2)
    ParametrizedSegment(f, t1, f(t1), t2, f(t2))
end

"Like `Box` but mutable."
mutable struct MutableBox{V}
    # Implicit: All lengths are equal
    mins::Vector{V}
    maxs::Vector{V}
end

"Iterator for the widths of each dimension. An iterator of unitful numbers."
_widths(bbox::MutableBox) = (ma - mi for (ma, mi) in zip(bbox.maxs, bbox.mins))

function MutableBox(vs)
    vsVectors = map(to, vs)
    ixs = eachindex(first(vsVectors)) |> collect
    mins = [minimum(v[i] for v in vsVectors) for i in ixs]
    maxs = [maximum(v[i] for v in vsVectors) for i in ixs]
    MutableBox(mins, maxs)
end

"All lengths must be the same and v must use linear indexing"
function _pushBBox!(bbox::MutableBox, v)
    vVector = to(v)
    for i in eachindex(vVector)
        bbox.mins[i] = min(bbox.mins[i], vVector[i])
        bbox.maxs[i] = max(bbox.maxs[i], vVector[i])
    end
end

"Check error value against tolerance and push if above."
function _maybePushErrfun!(queue, errfun, tol, x, args...)
    err = errfun(x, args...)
    if err > tol
        push!(queue, (-err, x))
    end
end

"An `errfun` for the sampling functions. l-infinity error relative to the range of the bounding box, component-wise."
function errRelativeBBox(s::ParametrizedSegment, bbox::MutableBox)
    # unitful handling; probably not the most efficient way.
    theUnit = unit(bbox.mins[1])
    # eps() is for handling of zero widths.
    # SOMEDAY is this good enough or actually still kinda unstable?
    ws = max.(_widths(bbox), eps() * theUnit)

    errs = @. (to(s.vMid) - 0.5 * (to(s.v1) + to(s.v2))) / ws
    norm(errs, Inf)
end

# @test errRelativeBBox(
#     Segment(t -> (t, t^2), 0.0, 10.0),
#     Range([0.0, 0.0], [1000.0, 1000.0]),
# ) == 0.025

# --------
# sampling
# --------

function sample(::AbstractRNG, geom::Geometry, method::AdaptiveSampling{E}) where {E}
    assertion(
        paramdim(geom) == 1,
        "Not implemented: Adaptive sampling for > 1 parameter dimension",
    )

    T = numtype(lentype(geom))
    tMin = zero(T)
    tMax = one(T)

    tsInit = range(tMin, tMax, length = method.minPoints)
    vsInit = geom.(tsInit)

    # We store initial points here and append to it as we go (out of order, sorted later)
    # SOMEDAY smarter accounting to avoid the final sort?
    # We could use the Mesh structure here to insert points as we go, also towards a multi-dim generalization.
    ps = zip(tsInit, vsInit) |> collect

    # Initialize bounding box. We'll update as we go.
    bbox = MutableBox(vsInit)

    # Queue of pending splits, in error order, so we can stick to `maxPoints`. We sort again later.
    # Initialize with pairs of successive points.
    # NB we can almost use RegularSampling for this but we also need the t values.
    queue = BinaryHeap{Tuple{E,ParametrizedSegment}}(Base.By(first), [])
    for ((t1, v1), (t2, v2)) in IterTools.partition(ps, 2, 1)
        s = ParametrizedSegment(geom, t1, v1, t2, v2)
        _maybePushErrfun!(queue, method.errfun, method.tol, s, bbox)
    end

    # Refinement.
    while !isempty(queue) && length(ps) < method.maxPoints
        (_, s) = pop!(queue)
        push!(ps, (s.tMid, s.vMid))
        s1 = ParametrizedSegment(geom, s.t1, s.v1, s.tMid, s.vMid)
        s2 = ParametrizedSegment(geom, s.tMid, s.vMid, s.t2, s.v2)
        _pushBBox!(bbox, s1.vMid)
        _pushBBox!(bbox, s2.vMid)

        _maybePushErrfun!(queue, method.errfun, method.tol, s1, bbox)
        _maybePushErrfun!(queue, method.errfun, method.tol, s2, bbox)
    end

    if !isempty(queue)
        maxErr = first(queue)[1]
        @warn "maxPoints=$(method.maxPoints) reached with errors above tolerance. Max error = $(-maxErr) > $(method.tol) = tolerance"
    end

    sort!(ps; by = first)
    return [p[2] for p in ps]
end

