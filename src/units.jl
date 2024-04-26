# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# helper type alias
const Len{T} = Quantity{T,u"𝐋"}
const Area{T} = Quantity{T,u"𝐋^2"}
const Vol{T} = Quantity{T,u"𝐋^3"}
const Met{T} = Quantity{T,u"𝐋",typeof(u"m")}

"""
    addunit(x, u)

Adds the unit only if the argument is not a quantity, otherwise an error is thrown.
"""
addunit(x::Number, u) = x * u
addunit(x::AbstractArray{<:Number}, u) = x * u
addunit(::Quantity, _) = throw(ArgumentError("invalid units, please check the documentation"))
addunit(::AbstractArray{<:Quantity}, _) = throw(ArgumentError("invalid units, please check the documentation"))
