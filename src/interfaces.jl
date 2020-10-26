# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    decompose(T, meshable)

Decompose a `meshable` object into elements of type `T`.
"""
function decompose(::Type{T}, primitive) where {T}
    return collect_with_eltype(T, primitive)
end

function decompose(::Type{P}, primitive) where {P<:Point}
    return convert.(P, coordinates(primitive))
end

function decompose(::Type{F}, primitive) where {F<:AbstractFace}
    f = faces(primitive)
    f === nothing && return nothing
    return collect_with_eltype(F, f)
end

function collect_with_eltype(::Type{T}, iter) where {T}
    result = T[]
    for element in iter
        for telement in convert_simplex(T, element)
            push!(result, telement)
        end
    end
    return result
end
