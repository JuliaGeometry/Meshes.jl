# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ellipsoid(semiaxes, angles; convention=TaitBryanExtr)

An ellipsoid with `semiaxes` and `angles` according to the rotation `convention`.

- For 2D ellipses, there are two semiaxes and one rotation angle.
- For 3D ellipsoids, there are three semiaxes and three rotation angles.

The list of available conventions can be found with `subtypes(RotationConvention)`.
"""
struct Ellipsoid{S,A,C,M} <: MetricBall
  # input fields
  semiaxes::S
  angles::A
  convention::C

  # state fields
  metric::M
end

function Ellipsoid(semiaxes::S, angles::A; convention::C=TaitBryanExtr) where {S,A,C}
  Dim, nangles = length(semiaxes), length(angles)
  valid = (Dim == 3 && nangles == 3) || (Dim == 2 && nangles == 1)
  @assert valid "invalid number of semiaxes/angles"

  # invert semiaxes if necessary
  invert = mainaxis(convention) == :Y
  ranges = invert ? [semiaxes[i] for i in reverse(1:Dim,1,2)] : semiaxes

  # scaling matrix
  Λ = Diagonal(SVector{Dim}(one(eltype(ranges))./ranges.^2))

  # rotation matrix
  R = rotmat(angles, convention)

  # ellipsoid metric
  metric = Mahalanobis(Symmetric(R'*Λ*R))

  Ellipsoid{S,A,C,typeof(metric)}(semiaxes, angles, convention, metric)
end

metric(ellips::Ellipsoid) = ellips.metric

function Base.show(io::IO, ellips::Ellipsoid)
  s = ellips.semiaxes
  a = ellips.angles
  c = ellips.convention
  T = eltype(s)
  print(io, "Ellipsoid{$T}($s, $a, $c)")
end
