@testset "Boolean Operations" setup = [Setup] begin
  p1 = cart.([(0, 0), (4, 0), (4, 4), (0, 4)])
  p2 = cart.([(2, 2), (6, 2), (6, 6), (2, 6)])
  poly1 = PolyArea([Ring(p1)])
  poly2 = PolyArea([Ring(p2)])

  unionpoly = poly1 ∪ poly2
  @test issimple(unionpoly)
  @test length(rings(unionpoly)) == 1
  @test area(unionpoly) ≈ 28.0u"m^2"

  intersectionpoly = geobooleanop(poly1, poly2, GeoIntersection())
  @test issimple(intersectionpoly)
  @test length(rings(intersectionpoly)) == 1
  @test area(intersectionpoly) ≈ 4.0u"m^2"

  differencepoly = setdiff(poly1, poly2)
  @test issimple(differencepoly)
  @test length(rings(differencepoly)) == 1
  @test area(differencepoly) ≈ 12.0u"m^2"

  symdiffpoly = poly1 ⊻ poly2
  @test !issimple(symdiffpoly)
  @test length(rings(symdiffpoly)) == 2
  @test area(symdiffpoly) ≈ 24.0u"m^2"

  # test with self intersecting polygons
  p3 = cart.([(1, 1), (5, 1), (1, 5), (5, 5)])
  poly3 = PolyArea([Ring(p3)])
  unionpoly2 = poly1 ∪ poly3
  @test issimple(unionpoly2)
  @test area(unionpoly2) ≈ 19.5u"m^2"

  intersectionpoly2 = poly1 ∩ poly3
  @test issimple(intersectionpoly2)
  @test Set(vertices(intersectionpoly2)) == Set(cart.([(1, 1), (4, 1), (4, 2), (3, 3), (2, 4), (4, 4)]))

  differencepoly2 = setdiff(poly1, poly3)
  @test issimple(differencepoly2)
  @test Set(vertices(differencepoly2)) ==
        Set(cart.([(0, 0), (4, 0), (4, 4), (0, 4), (1, 1), (4, 1), (4, 2), (3, 3), (2, 4)]))

  symdiffpoly2 = poly1 ⊻ poly3
  @test !issimple(symdiffpoly2)
  @test Set(vertices(symdiffpoly2)) ==
        Set(cart.([(0, 0), (4, 0), (0, 4), (4, 4), (1, 1), (4, 1), (4, 2), (3, 3), (2, 4), (5, 1), (5, 5), (1, 5)]))

  poly = Triangle(cart(6, 2), cart(3, 5), cart(0, 2))
  other = Quadrangle(cart(5, 0), cart(5, 4), cart(0, 4), cart(0, 0))
  clipped = geobooleanop(poly, other, GeoIntersection())
  @test issimple(clipped)
  @test all(vertices(clipped) .≈ [cart(5, 3), cart(4, 4), cart(2, 4), cart(0, 2), cart(5, 2)])

  # octagon
  poly = Octagon(cart(8, -2), cart(8, 5), cart(2, 5), cart(4, 3), cart(6, 3), cart(4, 1), cart(2, 1), cart(2, -2))
  other = Quadrangle(cart(5, 0), cart(5, 4), cart(0, 4), cart(0, 0))
  clipped = geobooleanop(poly, other, GeoIntersection())
  @test issimple.(rings(clipped)) |> all
  @test Set(vertices(clipped)) == Set([
    cart(3, 4),
    cart(4, 3),
    cart(5, 3),
    cart(5, 4),
    cart(5, 2),
    cart(4, 1),
    cart(2, 1),
    cart(2, 0),
    cart(5, 0)
  ])

  # inside
  poly = Quadrangle(cart(1, 0), cart(1, 1), cart(0, 1), cart(0, 0))
  other = Quadrangle(cart(5, 0), cart(5, 4), cart(0, 4), cart(0, 0))
  clipped = geobooleanop(poly, other, GeoIntersection())
  @test issimple(clipped)
  @test all(vertices(clipped) .≈ vertices(poly))

  # outside
  poly = Quadrangle(cart(7, 6), cart(7, 7), cart(6, 7), cart(6, 6))
  other = Quadrangle(cart(5, 0), cart(5, 4), cart(0, 4), cart(0, 0))
  clipped = geobooleanop(poly, other, GeoIntersection())
  @test isnothing(clipped)

  # surrounded
  poly = Hexagon(cart(0, 2), cart(-2, 2), cart(-2, 0), cart(0, -2), cart(2, -2), cart(2, 0))
  other = Hexagon(cart(1, 0), cart(0, 1), cart(-1, 1), cart(-1, 0), cart(0, -1), cart(1, -1))
  clipped = geobooleanop(poly, other, GeoIntersection())
  @test issimple(clipped)
  @test all(vertices(clipped) .≈ vertices(other))

  # PolyArea with box
  outer = Ring(cart(8, 0), cart(4, 8), cart(2, 8), cart(-2, 0), cart(0, 0), cart(1, 2), cart(5, 2), cart(6, 0))
  inner = Ring(cart(4, 4), cart(2, 4), cart(3, 6))
  poly = PolyArea([outer, inner])
  other = Box(cart(0, 1), cart(3, 7))
  clipped = geobooleanop(poly, other, GeoIntersection())
  crings = rings(clipped)
  @test !issimple(clipped)
  @test all(
    vertices(crings[1]) .≈
    [cart(1.5, 7.0), cart(0.0, 4.0), cart(0.0, 1.0), cart(0.5, 1.0), cart(1.0, 2.0), cart(3.0, 2.0), cart(3.0, 7.0)]
  )
  @test all(vertices(crings[2]) .≈ [cart(3.0, 4.0), cart(2.0, 4.0), cart(3.0, 6.0)])

  # PolyArea with outer ring outside and inner ring inside
  outer = Ring(cart(8, 0), cart(2, 6), cart(-4, 0))
  inner = Ring(cart(1, 3), cart(3, 3), cart(3, 1), cart(1, 1))
  poly = PolyArea([outer, inner])
  other = Quadrangle(cart(4, 4), cart(0, 4), cart(0, 0), cart(4, 0))
  clipped = geobooleanop(poly, other, GeoIntersection())
  @test !issimple(clipped)
  crings = rings(clipped)
  @test Set(vertices(crings[1])) == Set(vertices(other))
  @test all(vertices(crings[2]) .≈ vertices(inner))

  # PolyArea with one inner ring inside `other` and another inner ring outside `other`
  outer = Ring(cart(6, 4), cart(6, 7), cart(1, 6), cart(1, 1), cart(5, 2))
  inner₁ = Ring(cart(3, 3), cart(3, 4), cart(4, 3))
  inner₂ = Ring(cart(2, 5), cart(2, 6), cart(3, 5))
  poly = PolyArea([outer, inner₁, inner₂])
  other = PolyArea(Ring(cart(6, 1), cart(7, 2), cart(6, 5), cart(0, 2), cart(1, 1)))
  clipped = geobooleanop(poly, other, GeoIntersection())
  crings = rings(clipped)
  @test !issimple(clipped)
  @test length(crings) == 2
  @test all(vertices(crings[1]) .≈ [cart(6, 4), cart(6, 5), cart(1, 2.5), cart(1, 1), cart(5, 2)])
  @test all(vertices(crings[2]) .≈ [cart(3.0, 3.0), cart(3.0, 3.5), cart(10 / 3, 11 / 3), cart(4.0, 3.0)])

  # https://github.com/JuliaGeometry/Meshes.jl/issues/1218
  data1 = readdlm(joinpath(datadir, "issue1218-1.dat"), ',')
  data2 = readdlm(joinpath(datadir, "issue1218-2.dat"), ',')
  ring1 = Ring(cart.(data1[:, 1], data1[:, 2]))
  ring2 = Ring(cart.(data2[:, 1], data2[:, 2]))
  cring = geobooleanop(ring1, ring2, GeoIntersection())
  perim = length(cring[1])
  if T === Float32
    @test perim ≈ T(15880.919)u"m"
  elseif T === Float64
    @test perim ≈ T(15887.308996863363)u"m"
  end
end
