# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    AdaptiveDiscretization(; [options])  

Discretize one-dimensional parameteric geometries using [`AdaptiveSampling`](@ref). See there for details.
"""
struct AdaptiveDiscretization{TV} <: DiscretizationMethod
    adaptive_sampling::AdaptiveSampling{TV}
end

AdaptiveDiscretization(; tol=5e-4, errfun=err_relative_range, min_points=20.0, max_points=4_000.0) = AdaptiveDiscretization(AdaptiveSampling(; tol, errfun, min_points, max_points))

function discretize(geom::Geometry, method::AdaptiveDiscretization)
    # NB this is trivial right now b/c `AdaptiveSampling` we only support 1-dimensional lines right now.
    # TODO I should probably connect beginning and end if isperidic(geom), like in wrapgrid().
    points = sample(geom, method.adaptive_sampling)
    connec = [connect((i, i+1)) for i in 1:length(points)-1]
    SimpleMesh(points, connec)
end

