@testset "Clipping" begin
  @test "Triangle" begin
    poly = Triangle(P2(0,2), P2(3,5), P2(6,2))
    other = Quadrangle(P2(0,0), P2(0,4), P2(5,4), P2(5,0))

    newpoly = clip(poly, window, SutherlandHodgman())

    @test length(rings(newpoly)) == 1
    @test first(rings(clippedpoly)) ≈ Ring(P2[(0,2), (2,4), (4,4), (5,3), (5,2)])
  end

  @test "Random" begin
    poly = Ngon(rand(P2, 20)...)
    other = Quadrangle(P2(0,0), P2(0,1/2), P2(1/2,1/2), P2(1/2,0))

    inpoints = [p for p in vertices(poly) if p ∈ other]
    outpoints = [p for p in vertices(poly) if p ∉ other]

    newpoly = clip(poly, other, SutherlandHodgman())
    v = vertices(newpoly)

    @test newpoly isa PolyArea
    @test length(rings(newpoly)) == 1
    @test v ⊆ other
    @test inpoints ⊆ v
    @test isempty(outpoints ∩ v)
  end
end
