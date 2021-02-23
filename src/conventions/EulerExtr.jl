# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EulerExtr

Extrinsic right-handed rotation by the ZXZ axes.
"""
struct EulerExtr <: RotationConvention end
axesseq(::Type{EulerExtr}) = :ZXZ
orientation(::Type{EulerExtr}) = (:CCW,:CCW,:CCW)
angleunits(::Type{EulerExtr}) = :RAD
mainaxis(::Type{EulerExtr}) = :X
isextrinsic(::Type{EulerExtr}) = true
