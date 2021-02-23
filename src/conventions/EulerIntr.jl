# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EulerIntr

Intrinsic right-handed rotation by the ZXZ axes.
"""
struct EulerIntr <: RotationConvention end
axesseq(::Type{EulerIntr}) = :ZXZ
orientation(::Type{EulerIntr}) = (:CCW,:CCW,:CCW)
angleunits(::Type{EulerIntr}) = :RAD
mainaxis(::Type{EulerIntr}) = :X
isextrinsic(::Type{EulerIntr}) = false
