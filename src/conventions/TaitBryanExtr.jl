# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TaitBryanExtr

Extrinsic right-handed rotation by the ZXY axes.
"""
struct TaitBryanExtr <: RotationConvention end
axesseq(::Type{TaitBryanExtr}) = :ZXY
orientation(::Type{TaitBryanExtr}) = (:CCW,:CCW,:CCW)
angleunits(::Type{TaitBryanExtr}) = :RAD
mainaxis(::Type{TaitBryanExtr}) = :X
isextrinsic(::Type{TaitBryanExtr}) = true
