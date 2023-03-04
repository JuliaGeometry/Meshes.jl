# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# -------------
# 2D ROTATIONS
# -------------

# temporary fix for https://github.com/JuliaSpace/ReferenceFrameRotations.jl/issues/22 

"""
    ClockwiseAngle(θ)

Clockwise rotation in 2D space by angle `θ`.
"""
struct ClockwiseAngle{T}
  θ::T
end

Base.inv(cw::ClockwiseAngle) = CounterClockwiseAngle(cw.θ)

function Base.convert(::Type{<:DCM}, cw::ClockwiseAngle)
  s, c = sincos(cw.θ)
  T = typeof(s)
  return SMatrix{2,2,T}([c s; -s c])
end

"""
    CounterClockwiseAngle(θ)

Counter-clockwise rotation in 2D space by angle `θ`.
"""
struct CounterClockwiseAngle{T}
  θ::T
end

Base.inv(ccw::CounterClockwiseAngle) = ClockwiseAngle(ccw.θ)

function Base.convert(::Type{<:DCM}, ccw::CounterClockwiseAngle)
  return convert(DCM, ClockwiseAngle(-ccw.θ))
end

# -------------
# 3D ROTATIONS
# -------------

TaitBryanAngles(θs...) = inv(EulerAngles((-1 .* θs)..., :ZXY))
