"""
    ParametrizedCurve(a, b, func)

TODO
"""
struct ParametrizedCurve{M<:Meshes.Manifold,C<:Meshes.CRS,T<:Real,F<:Function} <: Primitive{M,C}
    a::T
    b::T
    func::F

    function ParametrizedCurve(a::T, b::T, func) where {T<:Real}
        p = func(a)
        M = manifold(p)
        C = crs(p)
        new{M, C, T, typeof(func)}(a, b, func)
    end
end

paramdim(::ParametrizedCurve) = 1
start_parameter(curve::ParametrizedCurve) = curve.a
start_point(curve::ParametrizedCurve) = curve(curve.a)
end_parameter(curve::ParametrizedCurve) = curve.b
end_point(curve::ParametrizedCurve) = curve(curve.b)

function (curve::ParametrizedCurve)(t)
    if t < start_parameter(curve) || t > end_parameter(curve)
        throw(DomainError(t, "c(t) is not defined for t outside [$(start_parameter(curve)), $(end_parameter(curve))]."))
    end
    curve.func(t)
end
