# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Datamine

Rotation convention of the Datamine software (fixed to ZXZ axes).
"""
struct Datamine <: RotationConvention end
axesseq(::Type{Datamine}) = :ZXZ
orientation(::Type{Datamine}) = (:CW,:CW,:CW)
angleunits(::Type{Datamine}) = :DEG
mainaxis(::Type{Datamine}) = :X
isextrinsic(::Type{Datamine}) = false
