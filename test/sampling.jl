@testset "Sampling" begin
  @testset "RegularSampler" begin
    b = Box(P2(0, 0), P2(2, 2))
    ps = sample(b, RegularSampler(3))
    @test collect(ps) == P2[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
    ps = sample(b, RegularSampler(2, 3))
    @test collect(ps) == P2[(0,0),(2,0),(0,1),(2,1),(0,2),(2,2)]

    s = Sphere(P2(0, 0), T(2))
    ps = sample(s, RegularSampler(4))
    ts = P2[(2,0),(0,2),(-2,0),(0,-2)]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    s = Sphere(P3(0, 0, 0), T(2))
    ps = sample(s, RegularSampler(2,2))
    ts = P3[(1.7320508075688772, 0.0, 1.0),
            (1.7320508075688772, 0.0, -1.0),
            (-1.7320508075688772, 0.0, 1.0),
            (-1.7320508075688772, 0.0, -1.0)]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end
  end
end
