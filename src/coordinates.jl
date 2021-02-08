# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    coordinates!(buff, object, inds)

Non-allocating version of [`coordinates`](@ref)
"""
function coordinates!(buff, obj, inds::AbstractVector{Int})
  for j in eachindex(inds)
    coordinates!(view(buff,:,j), obj, inds[j])
  end
  buff
end

"""
    coordinates(object, ind)

Return the coordinates of the `ind`-th element in the `object`.
"""
function coordinates(obj, ind::Int)
  buff = MVector{embeddim(obj),coordtype(obj)}(undef)
  coordinates!(buff, obj, ind)
end

"""
    coordinates(object, inds)

Return the coordinates of `inds` in the `object`.
"""
function coordinates(obj, inds::AbstractVector{Int})
  buff = Matrix{coordtype(obj)}(undef, embeddim(obj), length(inds))
  coordinates!(buff, obj, inds)
end
