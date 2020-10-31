for T in (Float32, Float64)
P2 = Point{2,T}
P3 = Point{3,T}
@testset "Sampling" begin
  @testset "RegularSampler" begin
    P = Point{2,T}
    b = Box(P(0, 0), P(2, 2))
    ps = sample(b, RegularSampler(3))
    @test collect(ps) == P[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
    ps = sample(b, RegularSampler(2, 3))
    @test collect(ps) == P[(0,0),(2,0),(0,1),(2,1),(0,2),(2,2)]
  end
end
end
