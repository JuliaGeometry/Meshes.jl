@testset "Clipping" begin
  @testset "SutherlandHodgman" begin
    # triangle
    poly = Triangle(point(6, 2), point(3, 5), point(0, 2))
    other = Quadrangle(point(5, 0), point(5, 4), point(0, 4), point(0, 0))
    clipped = clip(poly, other, SutherlandHodgman())
    @test issimple(clipped)
    @test all(vertices(clipped) .≈ [point(5, 3), point(4, 4), point(2, 4), point(0, 2), point(5, 2)])

    # octagon
    poly =
      Octagon(point(8, -2), point(8, 5), point(2, 5), point(4, 3), point(6, 3), point(4, 1), point(2, 1), point(2, -2))
    other = Quadrangle(point(5, 0), point(5, 4), point(0, 4), point(0, 0))
    clipped = clip(poly, other, SutherlandHodgman())
    @test !issimple(clipped)
    @test all(
      vertices(clipped) .≈ [
        point(3, 4),
        point(4, 3),
        point(5, 3),
        point(5, 2),
        point(4, 1),
        point(2, 1),
        point(2, 0),
        point(5, 0),
        point(5, 4)
      ]
    )

    # inside
    poly = Quadrangle(point(1, 0), point(1, 1), point(0, 1), point(0, 0))
    other = Quadrangle(point(5, 0), point(5, 4), point(0, 4), point(0, 0))
    clipped = clip(poly, other, SutherlandHodgman())
    @test issimple(clipped)
    @test all(vertices(clipped) .≈ vertices(poly))

    # outside
    poly = Quadrangle(point(7, 6), point(7, 7), point(6, 7), point(6, 6))
    other = Quadrangle(point(5, 0), point(5, 4), point(0, 4), point(0, 0))
    clipped = clip(poly, other, SutherlandHodgman())
    @test isnothing(clipped)

    # surrounded
    poly = Hexagon(point(0, 2), point(-2, 2), point(-2, 0), point(0, -2), point(2, -2), point(2, 0))
    other = Hexagon(point(1, 0), point(0, 1), point(-1, 1), point(-1, 0), point(0, -1), point(1, -1))
    clipped = clip(poly, other, SutherlandHodgman())
    @test issimple(clipped)
    @test all(vertices(clipped) .≈ vertices(other))

    # PolyArea with box
    outer =
      Ring(point(8, 0), point(4, 8), point(2, 8), point(-2, 0), point(0, 0), point(1, 2), point(5, 2), point(6, 0))
    inner = Ring(point(4, 4), point(2, 4), point(3, 6))
    poly = PolyArea([outer, inner])
    other = Box(point(0, 1), point(3, 7))
    clipped = clip(poly, other, SutherlandHodgman())
    crings = rings(clipped)
    @test !issimple(clipped)
    @test all(
      vertices(crings[1]) .≈ [
        point(1.5, 7.0),
        point(0.0, 4.0),
        point(0.0, 1.0),
        point(0.5, 1.0),
        point(1.0, 2.0),
        point(3.0, 2.0),
        point(3.0, 7.0)
      ]
    )
    @test all(vertices(crings[2]) .≈ [point(3.0, 4.0), point(2.0, 4.0), point(3.0, 6.0)])

    # PolyArea with outer ring outside and inner ring inside
    outer = Ring(point(8, 0), point(2, 6), point(-4, 0))
    inner = Ring(point(1, 3), point(3, 3), point(3, 1), point(1, 1))
    poly = PolyArea([outer, inner])
    other = Quadrangle(point(4, 4), point(0, 4), point(0, 0), point(4, 0))
    clipped = clip(poly, other, SutherlandHodgman())
    @test !issimple(clipped)
    crings = rings(clipped)
    @test all(vertices(crings[1]) .≈ vertices(other))
    @test all(vertices(crings[2]) .≈ vertices(inner))

    # PolyArea with one inner ring inside `other` and another inner ring outside `other`
    outer = Ring(point(6, 4), point(6, 7), point(1, 6), point(1, 1), point(5, 2))
    inner₁ = Ring(point(3, 3), point(3, 4), point(4, 3))
    inner₂ = Ring(point(2, 5), point(2, 6), point(3, 5))
    poly = PolyArea([outer, inner₁, inner₂])
    other = PolyArea(Ring(point(6, 1), point(7, 2), point(6, 5), point(0, 2), point(1, 1)))
    clipped = clip(poly, other, SutherlandHodgman())
    crings = rings(clipped)
    @test !issimple(clipped)
    @test length(crings) == 2
    @test all(vertices(crings[1]) .≈ [point(6, 4), point(6, 5), point(1, 2.5), point(1, 1), point(5, 2)])
    @test all(vertices(crings[2]) .≈ [point(3.0, 3.0), point(3.0, 3.5), point(10 / 3, 11 / 3), point(4.0, 3.0)])
  end
end
