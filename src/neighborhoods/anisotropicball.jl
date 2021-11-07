# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    AnisotropicBall(radii, angles; convention=TaitBryanExtr)

An anisotropic ball with `radii` and `angles` according to the rotation `convention`.

- For 2D balls (i.e. ellipses), there are two radii and one rotation angle.
- For 3D balls (i.e. ellipsoids), there are three radii and three rotation angles.

The list of available conventions can be found with `subtypes(RotationConvention)`.
"""
struct AnisotropicBall{V,A,C,M} <: MetricBall
  # input fields
  radii::V
  angles::A
  convention::C

  # state fields
  metric::M
end

function AnisotropicBall(radii::V, angles::A; convention::C=TaitBryanExtr) where {V,A,C}
  Dim, nangles = length(radii), length(angles)
  valid = (Dim == 3 && nangles == 3) || (Dim == 2 && nangles == 1)
  @assert valid "invalid number of radii/angles"

  # invert radii if necessary
  invert = mainaxis(convention) == :Y
  ranges = invert ? [radii[i] for i in reverse(1:Dim,1,2)] : radii

  # scaling matrix
  Λ = Diagonal(SVector{Dim}(one(eltype(ranges))./ranges.^2))

  # rotation matrix
  R = rotmat(angles, convention)

  # ellipsoid metric
  metric = Mahalanobis(Symmetric(R'*Λ*R))

  AnisotropicBall{V,A,C,typeof(metric)}(radii, angles, convention, metric)
end

metric(ball::AnisotropicBall) = ball.metric

radii(ball::AnisotropicBall) = ball.radii

function Base.show(io::IO, ball::AnisotropicBall)
  r = ball.radii
  a = ball.angles
  c = ball.convention
  T = eltype(r)
  print(io, "AnisotropicBall{$T}($r, $a, $c)")
end
