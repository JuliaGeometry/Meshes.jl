# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# LinearAlgebra function wrappers that handle units
# The result units of some operations, such as dot and cross, 
# are treated in a special way to handle Meshes.jl use cases

function usvd(A)
  u = unit(eltype(A))
  F = svd(ustrip.(A))
  SVD(F.U * u, F.S * u, F.Vt * u)
end

uinv(A) = inv(ustrip.(A)) * unit(eltype(A))^-1

unormalize(a::Vec{Dim,ℒ}) where {Dim,ℒ} = Vec(normalize(a) * unit(ℒ))

udot(a::Vec{Dim,ℒ}, b::Vec{Dim,ℒ}) where {Dim,ℒ} = ustrip(a ⋅ b) * unit(ℒ)

ucross(a::Vec{Dim,ℒ}, b::Vec{Dim,ℒ}) where {Dim,ℒ} = Vec(ustrip.(a × b) * unit(ℒ))
ucross(a::Vec{Dim,ℒ}, b::Vec{Dim,ℒ}, c::Vec{Dim,ℒ}) where {Dim,ℒ} = Vec(ustrip.(a × b × c) * unit(ℒ))
