# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TaitBryanIntr

Intrinsic right-handed rotation by the ZXY axes.
"""
struct TaitBryanIntr <: RotationConvention end
axesseq(::Type{TaitBryanIntr}) = :ZXY
orientation(::Type{TaitBryanIntr}) = (:CCW,:CCW,:CCW)
angleunits(::Type{TaitBryanIntr}) = :RAD
mainaxis(::Type{TaitBryanIntr}) = :X
isextrinsic(::Type{TaitBryanIntr}) = false
