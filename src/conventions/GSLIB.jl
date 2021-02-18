# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GSLIB

Rotation convention of the GSLIB software.
"""
struct GSLIB <: RotationConvention end
axesseq(::Type{GSLIB}) = :ZXY
orientation(::Type{GSLIB}) = (:CW,:CCW,:CCW)
angleunits(::Type{GSLIB}) = :DEG
mainaxis(::Type{GSLIB}) = :Y
isextrinsic(::Type{GSLIB}) = false
