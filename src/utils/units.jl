# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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
udot(a::Vec{Dim,ℒ₁}, b::Vec{Dim,ℒ₂}) where {Dim,ℒ₁,ℒ₂} = udot(promote(a, b)...)

ucross(a::Vec{Dim,ℒ}, b::Vec{Dim,ℒ}) where {Dim,ℒ} = Vec(ustrip.(a × b) * unit(ℒ))
ucross(a::Vec{Dim,ℒ₁}, b::Vec{Dim,ℒ₂}) where {Dim,ℒ₁,ℒ₂} = ucross(promote(a, b)...)

ucross(a::Vec{Dim,ℒ}, b::Vec{Dim,ℒ}, c::Vec{Dim,ℒ}) where {Dim,ℒ} = Vec(ustrip.(a × b × c) * unit(ℒ))

urotbetween(u::Vec, v::Vec) = rotation_between(ustrip.(u), ustrip.(v))

urotapply(R::Rotation, v::Vec{Dim,ℒ}) where {Dim,ℒ} = Vec(R * ustrip.(v) * unit(ℒ))
