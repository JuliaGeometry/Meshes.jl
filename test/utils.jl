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

@testitem "pairwiseintersect" setup = [Setup] begin
  # helper to sort points and seginds by point coordinates
  function sortedintersection(segs)
    points, inds = Meshes.pairwiseintersect(segs)
    perm = sortperm(points)
    points[perm], inds[perm]
  end

  # simple endpoint case
  segs = Segment.([(cart(0, 0), cart(2, 2)), (cart(0, 2), cart(2, 0)), (cart(0, 1), cart(0.5, 1))])
  points, seginds = sortedintersection(segs)
  @test length(points) == 1
  @test length(seginds) == 1

  # small number of segments, handling endpoints and precision
  segs = Segment.([(cart(0, 0), cart(2, 2)), (cart(1.5, 1), cart(2, 1)), (cart(1.51, 1.3), cart(2, 0.9))])
  points, seginds = sortedintersection(segs)
  @test length(points) == 1

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
  points, seginds = sortedintersection(segs)
  @test length(points) == 2
  @test length(seginds) == 2
  @test Set(seginds[1]) == Set([1, 2])
  @test Set(seginds[2]) == Set([1, 6, 5])

  # multiple intersections, endpoints as intersections
  if T === Float64
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
    points, seginds = sortedintersection(segs)
    @test length(points) == 4
    @test length(seginds) == 4
    @test points[3] ≈ cart(9, 4.8)
    @test points[4] ≈ cart(10, 4)
    @test Set(seginds[1]) == Set([4, 3])
    @test Set(seginds[2]) == Set([2, 3])
    @test Set(seginds[3]) == Set([4, 2])
    @test Set(seginds[4]) == Set([4, 5, 6, 7, 8])
  end

  # finds all intersections in a grid
  n = 10
  horizontal = [Segment(cart(1, i), cart(n, i)) for i in 1:n]
  vertical = [Segment(cart(i, 1), cart(i, n)) for i in 1:n]
  segs = [horizontal; vertical]
  points, seginds = sortedintersection(segs)
  @test length(points) == n * n - 4
  @test length(seginds) == n * n - 4
  @test Set(length.(seginds)) == Set([2])

  # number of intersections is invariant under rotations
  for θ in T(π / 6):T(π / 6):T(2π - π / 6)
    # rotation by π in Float32 is not robust, skips test
    T === Float32 && θ == T(π) && continue
    θpoints, θseginds = sortedintersection(segs |> Rotate(θ))
    @test length(θpoints) == n * n - 4
    @test length(θseginds) == n * n - 4
    @test Set(length.(θseginds)) == Set([2])
  end

  # tests coverage for when intervals don't overlap
  segs = [
    Segment(cart(0, 2), cart(2, 0)),
    Segment(cart(0, 0), cart(2, 2)),
    Segment(cart(3, 1), cart(3, 3)),
    Segment(cart(3, 3), cart(3, 3))
  ]
  points, seginds = sortedintersection(segs)
  @test length(points) == 1

  # inference test
  segs = facets(cartgrid(10, 10))
  @inferred Nothing (Meshes.pairwiseintersect(segs))
end

@testitem "isthreaded" setup = [Setup] begin
  if Threads.nthreads() > 1
    @test Meshes.isthreaded()
  end
  @test !Meshes.isthreaded(false)
end
