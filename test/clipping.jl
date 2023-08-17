@testset "Clipping" begin
  @testset "Triangle" begin
    poly = Triangle(P2(0, 2), P2(3, 5), P2(6, 2))
    other = Quadrangle(P2(0, 0), P2(0, 4), P2(5, 4), P2(5, 0))

    clipped = clip(poly, other, SutherlandHodgman())

    @test all(vertices(clipped) .≈ [
      P2(0, 2),
      P2(2, 4),
      P2(4, 4),
      P2(5, 3),
      P2(5, 2)
    ])
  end

  @testset "Octagon" begin
    poly = Octagon(
      P2(2,-2),
      P2(2, 1),
      P2(4, 1),
      P2(6, 3),
      P2(4, 3),
      P2(2, 5),
      P2(8, 5),
      P2(8,-2)
    )

    other = Quadrangle(P2(0, 0), P2(0, 4), P2(5, 4), P2(5, 0))
    
    clipped = clip(poly, other, SutherlandHodgman())

    @test all(vertices(clipped) .≈ [
      P2(2, 0),
      P2(2, 1),
      P2(4, 1),
      P2(5, 2),
      P2(5, 3),
      P2(4, 3),
      P2(3, 4),
      P2(5, 4),
      P2(5, 0)
    ])
  end

  @testset "Random" begin
    poly = Ngon(rand(P2, 10)...)
    other = Quadrangle(P2(0, 0), P2(0, 1/2), P2(1/2, 1/2), P2(1/2, 0))

    inverts = filter(p -> p ∈ other, vertices(poly))
    outverts = filter(p -> p ∉ other, vertices(poly))

    clipped = clip(poly, other, SutherlandHodgman())
    clipverts = vertices(clipped)

    @test clipverts ⊆ other
    @test inverts ⊆ clipverts
    @test isempty(outverts ∩ clipverts)
  end

  @testset "Inside" begin
    poly = Quadrangle(P2(0, 0), P2(0, 1), P2(1, 1), P2(1, 0))
    other = Quadrangle(P2(0, 0), P2(0, 4), P2(5, 4), P2(5, 0))

    clipped = clip(poly, other, SutherlandHodgman())

    @test length(rings(clipped)) == 1
    @test first(rings(clipped)) ≈ Ring(P2(0, 0), P2(0, 1), P2(1, 1), P2(1, 0))
  end

  @testset "Outside" begin
    poly = Quadrangle(P2(6, 6), P2(6, 7), P2(7, 7), P2(7, 6))
    other = Quadrangle(P2(0, 0), P2(0, 4), P2(5, 4), P2(5, 0))

    clipped = clip(poly, other, SutherlandHodgman())

    @test isnothing(clipped)
  end
end
