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
    poly = PolyArea([outer, inner])
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
    poly = PolyArea([outer, inner])
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
    poly = PolyArea([outer, inner₁, inner₂])
    other = PolyArea(Ring(P2(6, 1), P2(7, 2), P2(6, 5), P2(0, 2), P2(1, 1)))
    clipped = clip(poly, other, SutherlandHodgman())
    crings = rings(clipped)
    @test !issimple(clipped)
    @test length(crings) == 2
    @test all(vertices(crings[1]) .≈ [P2(6, 4), P2(6, 5), P2(1, 2.5), P2(1, 1), P2(5, 2)])
    @test all(vertices(crings[2]) .≈ [P2(3.0, 3.0), P2(3.0, 3.5), P2(10 / 3, 11 / 3), P2(4.0, 3.0)])
  end

  @testset "WeilerAtherton" begin
    # Triangle
    poly = Triangle(P2(4,2), P2(8,2), P2(6,4))
    other = Triangle(P2(0,0), P2(8,0), P2(4,4))
    clipped = clip(poly, other, WeilerAtherton())

    @test all(vertices(clipped) .≈ [P2(5,3), P2(4, 2), P2(6, 2)])

    # Hexagon
    poly = Hexagon(P2(2,10), P2(0,10), P2(5,5), P2(0,0), P2(2,0), P2(7,5))
    other = Hexagon(P2(0,5), P2(5,0), P2(7,0), P2(2,5), P2(7,10), P2(5,10))
    clipped = clip(poly, other, WeilerAtherton())
    clippedgeoms = clipped |> collect

    @test all(vertices(clippedgeoms[1]) .≈ [P2(2.5, 7.5), P2(3.5, 6.5), P2(4.5, 7.5), P2(3.5, 8.5)])
    @test all(vertices(clippedgeoms[2]) .≈ [P2(3.5, 3.5), P2(2.5, 2.5), P2(3.5, 1.5), P2(4.5, 2.5)])

    # PolyArea
    poly = PolyArea(P2(3,13), P2(1,13), P2(1,1), P2(3,1))
    other = PolyArea(P2(8, 14), P2(0, 14), P2(0, 0), P2(8, 0), P2(8, 2), P2(2, 2), P2(2, 6), P2(8, 6), P2(8, 8), P2(2, 8), P2(2, 12), P2(8, 12),)
    clipped = clip(poly, other, WeilerAtherton())

    @test all(vertices(clipped) .≈ [P2(3, 6), P2(3, 8), P2(2, 8), P2(2, 12), P2(3, 12), P2(3, 13), P2(1, 13), P2(1, 1), P2(3, 1), P2(3, 2), P2(2, 2), P2(2, 6)])

    # Petagon with PolyArea
    poly = Pentagon(P2(7,13), P2(1,10), P2(1,1), P2(7,1), P2(2,7))
    other = PolyArea(P2(8, 14), P2(0, 14), P2(0, 0), P2(8, 0), P2(8, 2), P2(2, 2), P2(2, 6), P2(8, 6), P2(8, 8), P2(2, 8), P2(2, 12), P2(8, 12))
    clipped = clip(poly, other, WeilerAtherton())
    clippedgeoms = clipped |> collect

    @test all(vertices(clippedgeoms[1]) .≈ [P2(2, 10.5), P2(1, 10), P2(1, 1), P2(7, 1), P2(37/6, 2), P2(2, 2), P2(2, 6), P2(17/6, 6), P2(2, 7), P2(17/6, 8), P2(2, 8)])
    @test all(vertices(clippedgeoms[2]) .≈ [P2(37/6, 12), P2(7, 13), P2(5, 12)])

    # PolyArea
    poly = PolyArea(P2(8,3), P2(5,6), P2(3,4), P2(3,1), P2(7,-1), P2(5,1), P2(7,1), P2(5,3))
    other = PolyArea(P2(10, 3), P2(7, 6), P2(7, 2), P2(2, 2), P2(0, 0), P2(7, 0))
    clipped = clip(poly, other, WeilerAtherton())
    clippedgeoms = clipped |> collect

    @test all(vertices(clippedgeoms[1]) .≈ [P2(3, 2), P2(3, 1), P2(5, 0), P2(6, 0), P2(5, 1), P2(7, 1), P2(6, 2)])
    @test all(vertices(clippedgeoms[2]) .≈ [P2(7, 3), P2(8, 3), P2(7, 4)])

    # outside
    poly = Quadrangle(P2(7, 6), P2(7, 7), P2(6, 7), P2(6, 6))
    other = Quadrangle(P2(5, 0), P2(5, 4), P2(0, 4), P2(0, 0))
    clipped = clip(poly, other, SutherlandHodgman())
    @test isnothing(clipped)

    # inside
    poly = Pentagon(P2(3,1), P2(3,2), P2(2,3), P2(1,2), P2(1,1))
    other = Quadrangle(P2(4,0), P2(4,4), P2(0,4), P2(0,0))
    clipped = clip(poly, other, WeilerAtherton())
    @test all(vertices(clipped) .≈ [P2(3, 1), P2(3, 2), P2(2, 3), P2(1, 2), P2(1, 1)])

    # intersection 
    poly = Triangle(P2(3,1), P2(1,1), P2(1,3))
    other = Triangle(P2(4,0), P2(0,4), P2(0,0))
    clipped = clip(poly, other, WeilerAtherton())
    @test all(vertices(clipped) .≈ [P2(3, 1), P2(1, 3), P2(1, 1)])

    # intersection 
    poly = Hexagon(P2(4,1), P2(4,3), P2(3,4), P2(1,4), P2(1,2), P2(2,1))
    other = Triangle(P2(5,0), P2(0,5), P2(0,0))
    clipped = clip(poly, other, WeilerAtherton())
    @test all(vertices(clipped) .≈ [P2(1, 4),P2(1, 2),P2(2, 1),P2(4, 1)])
    testviz(poly, other, clipped; alpha=0.8)

    # point intersections
    poly = PolyArea(P2(4,9), P2(2,6), P2(4,3), P2(6,6))
    other = PolyArea(P2(7,7), P2(6,8), P2(2,8), P2(0,7), P2(0,2), P2(1,1), P2(2,2), P2(6,2), P2(8,3), P2(6,4), P2(2,4), P2(1,5), P2(2,6), P2(6,6))
    clipped = clip2(poly, other, WeilerAtherton())
    clippedgeoms = clipped |> collect

    r = clippedgeoms[1] |> rings
    @test all(vertices(r[1]) .≈ [P2(10/3, 8.0), P2(2.0, 6.0), P2(6.0, 6.0), P2(14/3, 8.0)])
    r = clippedgeoms[2] |> rings
    @test all(vertices(r[1]) .≈ [P2(10/3, 4.0), P2(4.0, 3.0), P2(14/3, 4.0)])

    # degenerated
    poly = PolyArea(P2(16,16), P2(-4,16), P2(-4,2), P2(16,2), P2(16,4), P2(-2,4), P2(-2,6), P2(16,6), P2(16,8), P2(-2,8), P2(-2,10), P2(16,10), P2(16,12), P2(-2,12), P2(-2,14), P2(16,14))
    other = PolyArea(P2(12,16), P2(8,16), P2(8,10), P2(4,10), P2(4,16), P2(0,16), P2(0,0), P2(4,0), P2(4,6), P2(8,6), P2(8,0), P2(12,0))
    clipped = clip(poly, other, WeilerAtherton())
    clippedgeoms = clipped |> collect

    @test all(vertices(clippedgeoms[1]) .≈ [P2(12, 16), P2(8, 16), P2(8, 14), P2(12, 14)])
    @test all(vertices(clippedgeoms[2]) .≈ [P2(4, 16), P2(0, 16), P2(0, 14), P2(4, 14)])
    @test all(vertices(clippedgeoms[3]) .≈ [P2(0, 2), P2(4, 2), P2(4, 4), P2(0, 4)])
    @test all(vertices(clippedgeoms[4]) .≈ [P2(8, 2), P2(12, 2), P2(12, 4), P2(8, 4)])
    @test all(vertices(clippedgeoms[5]) .≈ [P2(0, 6), P2(4, 6), P2(8, 6), P2(12, 6), P2(12, 8), P2(0, 8)])
    @test all(vertices(clippedgeoms[6]) .≈ [P2(0, 10),P2(4, 10),P2(4, 12),P2(0, 12)])
    @test all(vertices(clippedgeoms[7]) .≈ [P2(8, 10),P2(12, 10),P2(12, 12),P2(8, 12)])

    # PolyArea with hole
    outer = Ring(P2(4,6), P2(0,6), P2(0,0), P2(4,0))
    inner = Ring(P2(3,1), P2(1,1), P2(1,5), P2(3,5))
    poly = PolyArea(outer, inner)
    other = PolyArea(P2(7,6), P2(-1,6), P2(-1,4), P2(5,4), P2(5,2), P2(-1,2), P2(-1,0), P2(7,0))
    clipped = clip(poly, other, WeilerAtherton())
    clippedgeoms = clipped |> collect

    r = clippedgeoms[1] |> rings
    @test all(vertices(r[1]) .≈ [P2(0, 2), P2(0, 0), P2(4, 0), P2(4, 2)])
    @test all(vertices(r[2]) .≈ [P2(3, 2), P2(3, 1), P2(1, 1), P2(1, 2)])

    r = clippedgeoms[2] |> rings
    @test all(vertices(r[1]) .≈ [P2(4, 4), P2(4, 6), P2(0, 6), P2(0, 4)])
    @test all(vertices(r[2]) .≈ [P2(1, 4), P2(1, 5), P2(3, 5), P2(3, 4)])

    # PolyArea with hole inside and outside
    outer = Ring(P2(7,1), P2(7,4), P2(4,5), P2(1,4), P2(1,1))
    inner1 = Ring(P2(6,2), P2(6,3), P2(5,3), P2(5,2))
    inner2 = Ring(P2(3,2), P2(3,3), P2(2,3), P2(2,2))
    poly = PolyArea(outer, inner1, inner2)
    other = PolyArea(P2(4,0), P2(4,4), P2(0,4), P2(0,0))
    clipped = clip(poly, other, WeilerAtherton())
    clippedgeoms = clipped |> collect

    r = clippedgeoms[1] |> rings
    @test all(vertices(r[1]) .≈ [P2(4, 4),P2(4, 6),P2(0, 6),P2(0, 4)])
    @test all(vertices(r[2]) .≈ [P2(1, 4), P2(1, 5), P2(3, 5), P2(3, 4)])

    # PolyArea outside with hole
    outer = Ring(P2(11,1), P2(8,5), P2(5,1))
    inner = Ring(P2(9,2), P2(8,4), P2(7,2))
    poly = PolyArea(outer, inner)
    other = PolyArea(P2(4,0), P2(4,4), P2(0,4), P2(0,0))
    clipped = clip(poly, other, WeilerAtherton())

    @test isnothing(clipped)
  end
end
