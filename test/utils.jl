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

@testitem "bentleyottmann" setup = [Setup] begin
  # simple endpoint case
  segs = Segment.([(cart(0, 0), cart(2, 2)), (cart(0, 2), cart(2, 0)), (cart(0, 1), cart(0.5, 1))])
  points, seginds = Meshes.bentleyottmann(segs)
  @test length(points) == 7
  @test length(seginds) == 7

  # small number of segments, handling endpoints and precision
  segs = Segment.([(cart(0, 0), cart(2, 2)), (cart(1.5, 1), cart(2, 1)), (cart(1.51, 1.3), cart(2, 0.9))])
  points, seginds = Meshes.bentleyottmann(segs)

  @test length(points) == 7

  # box case with one segment outside
  segs =
    Segment.([
      (cart(0, 0), cart(1.1, 1.1)),
      (cart(1, 0), cart(0, 1)),
      (cart(0, 0), cart(0, 1)),
      (cart(0, 0), cart(1, 0)),
      (cart(0, 1), cart(1, 1)),
      (cart(1, 0), cart(1, 1))
    ])
  points, seginds = Meshes.bentleyottmann(segs)
  @test length(
    setdiff(Set([cart(0, 0), cart(1, 1), cart(1.1, 1.1), cart(0, 1), cart(1, 0), cart(0.5, 0.5)]), Set(points))
  ) == 0
  @test length(points) == 6
  @test length(seginds) == 6
  inds = Dict(p => i for (i, p) in enumerate(points))
  @test Set(seginds[inds[cart(0.5, 0.5)]]) == Set([1, 2])
  @test Set(seginds[inds[cart(1, 1)]]) == Set([1, 6, 5])

  # multiple intersections, endpoints as intersections
  # in FP32, the outputs are correct, but spread over multiple points
  # off by 1e-7. All segments are still found
  segs =
    Segment.([
      (cart(9, 13), cart(6, 9)),
      (cart(2, 12), cart(9, 4.8)),
      (cart(12, 11), cart(4, 7)),
      (cart(2.5, 10), cart(12.5, 2)),
      (cart(13, 6), cart(10, 4)),
      (cart(10.5, 5.5), cart(9, 1)),
      (cart(10, 4), cart(11, -1)),
      (cart(10, 3), cart(10, 5))
    ])
  points, seginds = Meshes.bentleyottmann(segs)
  if T != Float32
    @test length(points) == 17
    @test length(seginds) == 17
    @test Set(seginds[findfirst(p -> p ≈ cart(10, 4), points)]) == Set([4, 5, 6, 7, 8])
    @test Set(seginds[findfirst(p -> p ≈ cart(9, 4.8), points)]) == Set([4, 2])
  end
  @test Set(reduce(vcat, seginds)) == Set(1:8)
  @test points[findfirst(p -> p ≈ cart(10, 4), points)] ≈ cart(10, 4)

  # finds all intersections in a grid
  n = 10
  horizontal = [Segment(cart(1, i), cart(n, i)) for i in 1:n]
  vertical = [Segment(cart(i, 1), cart(i, n)) for i in 1:n]
  segs = [horizontal; vertical]
  points, seginds = Meshes.bentleyottmann(segs)
  @test length(points) == 100
  @test length(seginds) == 100
  @test Set(length.(seginds)) == Set([2])

  # result is invariant under rotations
  for θ in T(π / 6):T(π / 6):T(2π - π / 6)
    θsegs = segs |> Rotate(θ)
    # Rotation by π in FP32 is not robust, skips test
    # technically the bentley-ottmann algorithm as implemented cant handle
    # infinitessimally off vertical segments with middle intersections.
    # so we snap the segments to the nearest integer coordinates
    if θ isa Float32 && θ == T(π)
      θsegs = map(segs) do s
        a, b = vertices(s)
        â = Meshes.coordround(a; digits=0)
        b̂ = Meshes.coordround(b; digits=0)
        Segment(â, b̂)
      end
    end
    θpoints, θseginds = Meshes.bentleyottmann(θsegs)
    @test length(θpoints) == 100
    @test length(θseginds) == 100
    @test Set(length.(θseginds)) == Set([2])
  end

  # inference test
  segs = facets(cartgrid(10, 10))
  @inferred Meshes.bentleyottmann(segs)
end

@testitem "isthreaded" setup = [Setup] begin
  if Threads.nthreads() > 1
    @test Meshes.isthreaded()
  end
end
