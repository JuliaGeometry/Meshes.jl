# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Leapfrog

Rotation convention of the Leapfrog software.
"""
struct Leapfrog <: RotationConvention end
axesseq(::Type{Leapfrog}) = :ZXZ
orientation(::Type{Leapfrog}) = (:CW,:CW,:CW)
angleunits(::Type{Leapfrog}) = :DEG
mainaxis(::Type{Leapfrog}) = :X
isextrinsic(::Type{Leapfrog}) = false
