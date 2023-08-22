@testset "Clipping" begin
  @testset "SutherlandHodgman" begin
    # triangle
    poly = Triangle(P2(6, 2), P2(3, 5), P2(0, 2))
    other = Quadrangle(P2(5, 0), P2(5, 4), P2(0, 4), P2(0, 0))
    clipped = clip(poly, other, SutherlandHodgman())
    @test issimple(clipped)
    @test all(vertices(clipped) .≈ [P2(5, 3), P2(4, 4), P2(2, 4), P2(0, 2), P2(5, 2)])

    # octagon
    poly = Octagon(P2(8, -2), P2(8, 5), P2(2, 5), P2(4, 3), P2(6, 3), P2(4, 1), P2(2, 1), P2(2, -2))
    other = Quadrangle(P2(5, 0), P2(5, 4), P2(0, 4), P2(0, 0))
    clipped = clip(poly, other, SutherlandHodgman())
    @test !issimple(clipped)
    @test all(
      vertices(clipped) .≈ [P2(3, 4), P2(4, 3), P2(5, 3), P2(5, 2), P2(4, 1), P2(2, 1), P2(2, 0), P2(5, 0), P2(5, 4)]
    )

    # inside
    poly = Quadrangle(P2(1, 0), P2(1, 1), P2(0, 1), P2(0, 0))
    other = Quadrangle(P2(5, 0), P2(5, 4), P2(0, 4), P2(0, 0))
    clipped = clip(poly, other, SutherlandHodgman())
    @test issimple(clipped)
    @test all(vertices(clipped) .≈ vertices(poly))

    # outside
    poly = Quadrangle(P2(7, 6), P2(7, 7), P2(6, 7), P2(6, 6))
    other = Quadrangle(P2(5, 0), P2(5, 4), P2(0, 4), P2(0, 0))
    clipped = clip(poly, other, SutherlandHodgman())
    @test isnothing(clipped)

    # surrounded
    poly = Hexagon(P2(0, 2), P2(-2, 2), P2(-2, 0), P2(0, -2), P2(2, -2), P2(2, 0))
    other = Hexagon(P2(1, 0), P2(0, 1), P2(-1, 1), P2(-1, 0), P2(0, -1), P2(1, -1))
    clipped = clip(poly, other, SutherlandHodgman())
    @test issimple(clipped)
    @test all(vertices(clipped) .≈ vertices(other))

    # PolyArea with box
    outer = Ring(P2(8, 0), P2(4, 8), P2(2, 8), P2(-2, 0), P2(0, 0), P2(1, 2), P2(5, 2), P2(6, 0))
    inner = Ring(P2(4, 4), P2(2, 4), P2(3, 6))
    poly = PolyArea(outer, [inner])
    other = Box(P2(0, 1), P2(3, 7))
    clipped = clip(poly, other, SutherlandHodgman())
    crings = rings(clipped)
    @test !issimple(clipped)
    @test all(
      vertices(crings[1]) .≈
      [P2(1.5, 7.0), P2(0.0, 4.0), P2(0.0, 1.0), P2(0.5, 1.0), P2(1.0, 2.0), P2(3.0, 2.0), P2(3.0, 7.0)]
    )
    @test all(vertices(crings[2]) .≈ [P2(3.0, 4.0), P2(2.0, 4.0), P2(3.0, 6.0)])

    # PolyArea with outer ring outside and inner ring inside
    outer = Ring(P2(8, 0), P2(2, 6), P2(-4, 0))
    inner = Ring(P2(1, 3), P2(3, 3), P2(3, 1), P2(1, 1))
    poly = PolyArea(outer, [inner])
    other = Quadrangle(P2(4, 4), P2(0, 4), P2(0, 0), P2(4, 0))
    clipped = clip(poly, other, SutherlandHodgman())
    @test !issimple(clipped)
    crings = rings(clipped)
    @test all(vertices(crings[1]) .≈ vertices(other))
    @test all(vertices(crings[2]) .≈ vertices(inner))

    # PolyArea with one inner ring inside `other` and another inner ring outside `other`
    outer = Ring(P2(6, 4), P2(6, 7), P2(1, 6), P2(1, 1), P2(5, 2))
    inner₁ = Ring(P2(3, 3), P2(3, 4), P2(4, 3))
    inner₂ = Ring(P2(2, 5), P2(2, 6), P2(3, 5))
    poly = PolyArea(outer, [inner₁, inner₂])
    other = PolyArea(Ring(P2(6, 1), P2(7, 2), P2(6, 5), P2(0, 2), P2(1, 1)))
    clipped = clip(poly, other, SutherlandHodgman())
    crings = rings(clipped)
    @test !issimple(clipped)
    @test length(crings) == 2
    @test all(vertices(crings[1]) .≈ [P2(6, 4), P2(6, 5), P2(1, 2.5), P2(1, 1), P2(5, 2)])
    @test all(vertices(crings[2]) .≈ [P2(3.0, 3.0), P2(3.0, 3.5), P2(10 / 3, 11 / 3), P2(4.0, 3.0)])
  end

  @testset "WeilerAtherton" begin
    # triangle
    poly = Triangle(P2(6, 2), P2(3, 5), P2(3, 2))
    other = Quadrangle(P2(5, 0), P2(5, 4), P2(0, 4), P2(0, 0))
    clipped = clip(poly, other, WeilerAtherton())
    
    @test issimple(clipped)
    @test all(vertices(clipped) .≈ [P2(5, 3), P2(4, 4), P2(3, 4), P2(3, 2), P2(5,2)])
  end
end
