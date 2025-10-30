# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# helper type alias
const Len{T} = Quantity{T,u"𝐋"}
const Area{T} = Quantity{T,u"𝐋^2"}
const Vol{T} = Quantity{T,u"𝐋^3"}
const Met{T} = Quantity{T,u"𝐋",typeof(u"m")}
const Deg{T} = Quantity{T,NoDims,typeof(u"°")}

"""
    aslen(x)

Adds meter unit if the argument is not a quantity
with length unit, otherwise, returns `x` as is.
"""
aslen(x::Len) = x
aslen(x::Number) = x * u"m"
aslen(::Quantity) = throw(ArgumentError("invalid length unit"))

"""
    aslentype(T)

If `T` has length unit, return it. If `T` is number, return it with meter unit.
Otherwise, throw an error.
"""
aslentype(::Type{T}) where {T<:Len} = T
aslentype(::Type{T}) where {T<:Number} = Met{T}
aslentype(::Type{<:Quantity}) = throw(ArgumentError("invalid length unit"))

"""
    numconvert(T, x)

Converts the number type of quantity `x` to `T`.
"""
numconvert(::Type{T}, x::Quantity{S,D,U}) where {T,S,D,U} = convert(Quantity{T,D,U}, x)

"""
    withunit(x, u)

Adds the unit if the argument is not a quantity, 
otherwise, converts the unit of `x` to `u`.
"""
withunit(x::Number, u) = x * u
withunit(x::Quantity, u) = uconvert(u, x)
