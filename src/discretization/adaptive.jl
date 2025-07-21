# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    AdaptiveDiscretization(; [options])  

Discretize one-dimensional parameteric geometries using [`AdaptiveSampling`](@ref). See there for details.
"""
struct AdaptiveDiscretization{TV} <: DiscretizationMethod
    adaptiveSampling::AdaptiveSampling{TV}
end

AdaptiveDiscretization(;
    tol = 5e-4,
    errfun = errRelativeBBox,
    minPoints = 20.0,
    maxPoints = 4_000.0,
) = AdaptiveDiscretization(AdaptiveSampling(; tol, errfun, minPoints, maxPoints))

function discretize(geom::Geometry, method::AdaptiveDiscretization)
    # NB this is trivial right now b/c `AdaptiveSampling` we only support 1-dimensional lines right now.
    # TODO I should probably connect beginning and end if isperidic(geom), like in wrapgrid().
    ps = sample(geom, method.adaptiveSampling)
    connec = [connect((i, i + 1)) for i = 1:length(ps)-1]
    SimpleMesh(ps, connec)
end

