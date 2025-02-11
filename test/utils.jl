@testitem "Utilities" setup = [Setup] begin
  a, b, c = cart(0, 0), cart(1, 0), cart(0, 1)
  @test signarea(a, b, c) == T(0.5) * u"m^2"
  a, b, c = cart(0, 0), cart(0, 1), cart(1, 0)
  @test signarea(a, b, c) == T(-0.5) * u"m^2"

  normals = [
    vector(1, 0, 0),
    vector(0, 1, 0),
    vector(0, 0, 1),
    vector(-1, 0, 0),
    vector(0, -1, 0),
    vector(0, 0, -1),
    vector(ntuple(i -> rand() - 0.5, 3))
  ]
  for n in normals
    u, v = householderbasis(n)
    @test u isa Vec{3}
    @test v isa Vec{3}
    @test ustrip.(u × v) ≈ n ./ norm(n)
  end
  n = Vec(T(1) * u"cm", T(1) * u"cm", T(1) * u"cm")
  u, v = householderbasis(n)
  @test unit(eltype(u)) == u"cm"
  @test unit(eltype(v)) == u"cm"
  n = Vec(T(1) * u"km", T(1) * u"km", T(1) * u"km")
  u, v = householderbasis(n)
  @test unit(eltype(u)) == u"km"
  @test unit(eltype(v)) == u"km"

  @test Meshes.mayberound(1.1, 1.0, 0.2) ≈ 1.0
  @test Meshes.mayberound(1.1, 1.0, 0.10000000000000001) ≈ 1.1
  @test Meshes.mayberound(1.1, 1.0, 0.05) ≈ 1.1

  # intersect parameters
  p1, p2 = cart(0, 0), cart(1, 1)
  p3, p4 = cart(1, 0), cart(0, 1)
  @inferred Meshes.intersectparameters(p1, p2, p3, p4)
  @inferred Meshes.intersectparameters(p1, p3, p2, p4)
  @inferred Meshes.intersectparameters(p1, p2, p1, p2)

  p1, p2 = cart(0, 0, 0), cart(1, 1, 1)
  p3, p4 = cart(1, 0, 0), cart(0, 1, 1)
  @inferred Meshes.intersectparameters(p1, p2, p3, p4)
  @inferred Meshes.intersectparameters(p1, p3, p2, p4)
  @inferred Meshes.intersectparameters(p1, p2, p1, p2)

  # withcrs
  c = (T(1), T(1))
  p = merc(c)
  v = to(p)
  @inferred Meshes.withcrs(p, v)
  @inferred Meshes.withcrs(p, c)
  c = (T(30), T(60))
  p = latlon(c) |> Proj(Cartesian)
  @inferred Meshes.withcrs(p, c, LatLon)
end
