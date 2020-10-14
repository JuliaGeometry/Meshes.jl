
"""
    OffsetInteger{O, T}

OffsetInteger type mainly for indexing.
* `O` - The offset relative to Julia arrays. This helps reduce copying when
communicating with 0-indexed systems such as OpenGL.
"""
struct OffsetInteger{O,T<:Integer} <: Integer
    i::T
    OffsetInteger{O,T}(x::Integer) where {O,T<:Integer} = new{O,T}(T(x + O))
end

const ZeroIndex{T<:Integer} = OffsetInteger{-1,T}
const OneIndex{T<:Integer} = OffsetInteger{0,T}
const GLIndex = ZeroIndex{Cuint}

raw(x::OffsetInteger) = x.i
raw(x::Integer) = x
value(x::OffsetInteger{O,T}) where {O,T} = raw(x) - O
value(x::Integer) = x

function show(io::IO, oi::OffsetInteger)
    return print(io, "|$(raw(oi)) (indexes as $(value(oi))|")
end

Base.eltype(::Type{OffsetInteger{O,T}}) where {O,T} = T
Base.eltype(oi::OffsetInteger) = eltype(typeof(oi))

# constructors and conversion
function OffsetInteger{O1,T1}(x::OffsetInteger{O2,T2}) where {O1,O2,T1<:Integer,T2<:Integer}
    return OffsetInteger{O1,T1}(convert(T2, x))
end

OffsetInteger{O}(x::Integer) where {O} = OffsetInteger{O,eltype(x)}(x)
OffsetInteger{O}(x::OffsetInteger) where {O} = OffsetInteger{O,eltype(x)}(x)
# This constructor has a massive method invalidation as a consequence,
# and doesn't seem to be needed, so let's remove it!

Base.convert(::Type{IT}, x::OffsetInteger) where {IT<:Integer} = IT(value(x))

Base.promote_rule(::Type{IT}, ::Type{<:OffsetInteger}) where {IT<:Integer} = IT

function Base.promote_rule(::Type{OffsetInteger{O1,T1}},
                           ::Type{OffsetInteger{O2,T2}}) where {O1,O2,T1<:Integer,
                                                                T2<:Integer}
    return OffsetInteger{pure_max(O1, O2),promote_type(T1, T2)}
end

Base.@pure pure_max(x1, x2) = x1 > x2 ? x1 : x2

# Need to convert to Int here because of: https://github.com/JuliaLang/julia/issues/35038
Base.to_index(I::OffsetInteger) = convert(Int, raw(OneIndex(I)))
Base.to_index(I::OffsetInteger{0}) = convert(Int, raw(I))

# basic operators
for op in (:(-), :abs)
    @eval Base.$op(x::T) where {T<:OffsetInteger} = T($(op)(value(x)))
end

for op in (:(+), :(-), :(*), :(/), :div)
    @eval begin
        @inline function Base.$op(x::OffsetInteger{O}, y::OffsetInteger{O}) where {O}
            return OffsetInteger{O}($op(value(x), value(y)))
        end
    end
end

for op in (:(==), :(>=), :(<=), :(<), :(>), :sub_with_overflow)
    @eval begin
        @inline function Base.$op(x::OffsetInteger{O}, y::OffsetInteger{O}) where {O}
            return $op(x.i, y.i)
        end
    end
end
