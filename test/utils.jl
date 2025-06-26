@testitem "signarea" setup = [Setup] begin
  a, b, c = cart(0, 0), cart(1, 0), cart(0, 1)
  @test signarea(a, b, c) == T(0.5) * u"m^2"
  a, b, c = cart(0, 0), cart(0, 1), cart(1, 0)
  @test signarea(a, b, c) == T(-0.5) * u"m^2"
end

@testitem "householderbasis" setup = [Setup] begin
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
end

@testitem "mayberound" setup = [Setup] begin
  @test Meshes.mayberound(1.1, 1.0, 0.2) ≈ 1.0
  @test Meshes.mayberound(1.1, 1.0, 0.10000000000000001) ≈ 1.1
  @test Meshes.mayberound(1.1, 1.0, 0.05) ≈ 1.1
end

@testitem "intersectparameters" setup = [Setup] begin
  # https://github.com/JuliaGeometry/Meshes.jl/issues/1218
  if T === Float64
    p1 = cart(387843.1300172474, 7.648008470021072e6)
    p2 = cart(387526.44396928686, 7.647621555327687e6)
    p3 = cart(387732.29, 7.64787305e6)
    p4 = cart(387676.87, 7.64780534e6)
    @test Meshes.intersectparameters(p1, p2, p3, p4) == (T(0), T(0), 1, 1)
    @test Meshes.intersectparameters(p3, p4, p1, p2) == (T(0), T(0), 1, 1)
  end

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
end

@testitem "withcrs" setup = [Setup] begin
  c = (T(1), T(1))
  p = merc(c)
  v = to(p)
  @inferred Meshes.withcrs(p, v)
  @inferred Meshes.withcrs(p, c)
  c = (T(30), T(60))
  p = latlon(c) |> Proj(Cartesian)
  @inferred Meshes.withcrs(p, c, LatLon)
end

@testitem "coordround" setup = [Setup] begin
  p₁ = cart(1, 1)
  p₂ = cart(1.0000000000004, 0.9999999999996)
  @test Meshes.coordround(p₁, sigdigits=5) == p₁
  @test Meshes.coordround(p₂, digits=10) == p₁
  @inferred Meshes.coordround(p₁, digits=10)
end
