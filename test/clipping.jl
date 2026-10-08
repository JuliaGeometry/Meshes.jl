@testitem "SutherlandHodgman" setup = [Setup] begin
  # triangle
  poly = Triangle(cart(6, 2), cart(3, 5), cart(0, 2))
  other = Quadrangle(cart(5, 0), cart(5, 4), cart(0, 4), cart(0, 0))
  clipped = clip(poly, other, SutherlandHodgmanClipping())
  @test issimple(clipped)
  @test all(vertices(clipped) .≈ [cart(5, 3), cart(4, 4), cart(2, 4), cart(0, 2), cart(5, 2)])

  # octagon
  poly = Octagon(cart(8, -2), cart(8, 5), cart(2, 5), cart(4, 3), cart(6, 3), cart(4, 1), cart(2, 1), cart(2, -2))
  other = Quadrangle(cart(5, 0), cart(5, 4), cart(0, 4), cart(0, 0))
  clipped = clip(poly, other, SutherlandHodgmanClipping())
  @test !issimple(clipped)
  @test all(
    vertices(clipped) .≈
    [cart(3, 4), cart(4, 3), cart(5, 3), cart(5, 2), cart(4, 1), cart(2, 1), cart(2, 0), cart(5, 0), cart(5, 4)]
  )

  # inside
  poly = Quadrangle(cart(1, 0), cart(1, 1), cart(0, 1), cart(0, 0))
  other = Quadrangle(cart(5, 0), cart(5, 4), cart(0, 4), cart(0, 0))
  clipped = clip(poly, other, SutherlandHodgmanClipping())
  @test issimple(clipped)
  @test all(vertices(clipped) .≈ vertices(poly))

  # outside
  poly = Quadrangle(cart(7, 6), cart(7, 7), cart(6, 7), cart(6, 6))
  other = Quadrangle(cart(5, 0), cart(5, 4), cart(0, 4), cart(0, 0))
  clipped = clip(poly, other, SutherlandHodgmanClipping())
  @test isnothing(clipped)

  # surrounded
  poly = Hexagon(cart(0, 2), cart(-2, 2), cart(-2, 0), cart(0, -2), cart(2, -2), cart(2, 0))
  other = Hexagon(cart(1, 0), cart(0, 1), cart(-1, 1), cart(-1, 0), cart(0, -1), cart(1, -1))
  clipped = clip(poly, other, SutherlandHodgmanClipping())
  @test issimple(clipped)
  @test all(vertices(clipped) .≈ vertices(other))

  # PolyArea with Quadrangle
  outer = Ring(cart(8, 0), cart(4, 8), cart(2, 8), cart(-2, 0), cart(0, 0), cart(1, 2), cart(5, 2), cart(6, 0))
  inner = Ring(cart(4, 4), cart(2, 4), cart(3, 6))
  poly = PolyArea([outer, inner])
  other = Quadrangle(cart(0, 1), cart(3, 1), cart(3, 7), cart(0, 7))
  clipped = clip(poly, other, SutherlandHodgmanClipping())
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
  clipped = clip(poly, other, SutherlandHodgmanClipping())
  @test !issimple(clipped)
  crings = rings(clipped)
  @test all(vertices(crings[1]) .≈ vertices(other))
  @test all(vertices(crings[2]) .≈ vertices(inner))

  # PolyArea with one inner ring inside `other` and another inner ring outside `other`
  outer = Ring(cart(6, 4), cart(6, 7), cart(1, 6), cart(1, 1), cart(5, 2))
  inner₁ = Ring(cart(3, 3), cart(3, 4), cart(4, 3))
  inner₂ = Ring(cart(2, 5), cart(2, 6), cart(3, 5))
  poly = PolyArea([outer, inner₁, inner₂])
  other = PolyArea(Ring(cart(6, 1), cart(7, 2), cart(6, 5), cart(0, 2), cart(1, 1)))
  clipped = clip(poly, other, SutherlandHodgmanClipping())
  crings = rings(clipped)
  @test !issimple(clipped)
  @test length(crings) == 2
  @test all(vertices(crings[1]) .≈ [cart(6, 4), cart(6, 5), cart(1, 2.5), cart(1, 1), cart(5, 2)])
  @test all(vertices(crings[2]) .≈ [cart(3.0, 3.0), cart(3.0, 3.5), cart(10 / 3, 11 / 3), cart(4.0, 3.0)])

  # https://github.com/JuliaGeometry/Meshes.jl/issues/1218
  data1 = readdlm(joinpath(datadir, "issue1218-1.dat"), ',')
  data2 = readdlm(joinpath(datadir, "issue1218-2.dat"), ',')
  poly1 = PolyArea(cart.(data1[:, 1], data1[:, 2]))
  poly2 = PolyArea(cart.(data2[:, 1], data2[:, 2]))
  cpoly = clip(poly1, poly2, SutherlandHodgmanClipping())
  perim = perimeter(cpoly)
  if T === Float32
    @test perim ≈ T(15880.919)u"m"
  elseif T === Float64
    @test perim ≈ T(15887.308996863363)u"m"
  end
end

@testitem "GreinerHormann" setup = [Setup] begin
  # same result as Sutherland-Hodgman when the clipping geometry is convex
  poly = Triangle(cart(6, 2), cart(3, 5), cart(0, 2))
  other = Quadrangle(cart(5, 0), cart(5, 4), cart(0, 4), cart(0, 0))
  clipped = clip(poly, other, GreinerHormannClipping())
  @test issimple(clipped)
  @test all(vertices(clipped) .≈ [cart(5, 3), cart(4, 4), cart(2, 4), cart(0, 2), cart(5, 2)])

  # non-convex clipping geometry
  poly = PolyArea(cart.([(0, 0), (4, 0), (4, 1), (1, 1), (1, 4), (0, 4)]))
  other = PolyArea(cart.([(0, 0), (4, 0), (4, 4), (3, 4), (3, 1), (0, 1)]))
  clipped = clip(poly, other, GreinerHormannClipping())
  @test all(vertices(clipped) .≈ [cart(4, 1), cart(0, 1), cart(0, 0), cart(4, 0)])

  # clipped polygon with two components
  poly = Quadrangle(cart(0, 0), cart(6, 0), cart(6, 1), cart(0, 1))
  other = PolyArea(cart.([(0, 0.25), (2, 0.25), (2, 2), (4, 2), (4, 0.25), (6, 0.25), (6, 3), (0, 3)]))
  clipped = clip(poly, other, GreinerHormannClipping())
  @test clipped isa Multi
  @test length(parent(clipped)) == 2
  @test measure(clipped) ≈ T(3) * u"m^2"

  # polygon with hole
  outer = Ring(cart.([(0, 0), (10, 0), (10, 10), (0, 10)]))
  inner = Ring(cart.([(3, 3), (3, 7), (7, 7), (7, 3)]))
  poly = PolyArea([outer, inner])
  other = Quadrangle(cart(5, -2), cart(14, -2), cart(14, 12), cart(5, 12))
  clipped = clip(poly, other, GreinerHormannClipping())
  @test measure(clipped) ≈ T(42) * u"m^2"

  # vertices on edges and shared edges
  tri = Triangle(cart(2, 4), cart(6, 2), cart(6, 6))
  quad1 = Quadrangle(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  quad2 = Quadrangle(cart(4, 0), cart(8, 0), cart(8, 4), cart(4, 4))
  quad3 = Quadrangle(cart(4, 4), cart(8, 4), cart(8, 8), cart(4, 8))
  @test measure(clip(quad1, quad1, GreinerHormannClipping())) ≈ T(16) * u"m^2"
  @test isnothing(clip(quad2, quad1, GreinerHormannClipping()))
  @test isnothing(clip(quad3, quad1, GreinerHormannClipping()))
  clipped = clip(tri, quad1, GreinerHormannClipping())
  @test all(vertices(clipped) .≈ [cart(2, 4), cart(4, 3), cart(4, 4)])

  # inside and outside
  quad1 = Quadrangle(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  quad2 = Quadrangle(cart(1, 1), cart(3, 1), cart(3, 3), cart(1, 3))
  quad3 = Quadrangle(cart(10, 10), cart(11, 10), cart(11, 11), cart(10, 11))
  @test all(vertices(clip(quad2, quad1, GreinerHormannClipping())) .≈ vertices(quad2))
  @test isnothing(clip(quad3, quad1, GreinerHormannClipping()))

  # CRS propagation
  poly = Triangle(merc(6, 2), merc(3, 5), merc(0, 2))
  other = Quadrangle(merc(5, 0), merc(5, 4), merc(0, 4), merc(0, 0))
  @test crs(clip(poly, other, GreinerHormannClipping())) === crs(poly)

  # polygons with different machine precision
  tri = Triangle(Point(0.0, 0.0), Point(4.0, 0.0), Point(0.0, 4.0))
  quad = Quadrangle(Point(1.0f0, 1.0f0), Point(3.0f0, 1.0f0), Point(3.0f0, 3.0f0), Point(1.0f0, 3.0f0))
  clipped = clip(tri, quad, GreinerHormannClipping())
  @test Unitful.numtype(Meshes.lentype(clipped)) == Float64
end
