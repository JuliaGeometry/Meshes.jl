@testset "Sampling" begin
  @testset "RegularSampler" begin
    b = Box(P2(0, 0), P2(2, 2))
    ps = sample(b, RegularSampler(3))
    @test collect(ps) == P2[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
    ps = sample(b, RegularSampler(2, 3))
    @test collect(ps) == P2[(0,0),(2,0),(0,1),(2,1),(0,2),(2,2)]

    s = Sphere(P2(0, 0), T(2))
    ps = sample(s, RegularSampler(4))
    # @test collect(ps) == P2[(2,0),(0,2),(-2,0),(0,-2)]

    s = Sphere(P3(0, 0, 0), T(2))
    ps = sample(s, RegularSampler(2,2))
    # collect(ps) == P3[(1.7320509, 0.0, 0.99999994),
    #                   (1.7320508, 0.0, -1.0000001),
    #                   (-1.7320509, -1.5142069e-7, 0.99999994),
    #                   (-1.7320508, -1.5142069e-7, -1.0000001)]
  end
end
