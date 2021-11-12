@testset "Neighborhoods" begin
  @testset "MetricBall" begin
    # Euclidean metric
    b = MetricBall(T(1/2))
    m = metric(b)
    r = range(b)
    @test evaluate(m, T[0], T[0]) ≤ r
    @test evaluate(m, T[0], T[1]) > r

    b = MetricBall(T(1))
    m = metric(b)
    r = range(b)
    @test evaluate(m, T[0,0], T[0,0]) ≤ r
    @test evaluate(m, T[0,0], T[1,0]) ≤ r
    @test evaluate(m, T[0,0], T[0,1]) ≤ r

    # Chebyshev metric
    b = MetricBall(T(1/2), Chebyshev())
    m = metric(b)
    r = range(b)
    @test evaluate(m, T[0], T[0]) ≤ r
    @test evaluate(m, T[0], T[1]) > r

    for r in [1.,2.,3.,4.,5.]
      b = MetricBall(r, Chebyshev())
      m = metric(b)
      r = range(b)
      for i in 0.0:1.0:r, j in 0.0:1.0:r
        @test evaluate(m, T[0,0], T[i,j]) ≤ r
      end
    end

    # 2D ellipsoid rotated 45 degrees in GSLIB convention
    b = MetricBall(T[2,1], T[45], convention=GSLIB)
    m = metric(b)
    r = range(b)

    # tests along main semiaxes, slightly below threshold
    @test evaluate(m, T[0,0], T(1.9) * T[√2/2, √2/2]) ≤ r
    @test evaluate(m, T[0,0], T(0.9) * T[√2/2,-√2/2]) ≤ r

    # tests along main semiaxes, slightly above threshold
    @test evaluate(m, T[0,0], T(2.1) * T[√2/2, √2/2]) > r
    @test evaluate(m, T[0,0], T(1.1) * T[√2/2,-√2/2]) > r

    # 3D ellipsoid rotated (45, -45, 0) in GSLIB convention
    b = MetricBall(T[3,2,1], T[45,-45,0], convention=GSLIB)
    m = metric(b)
    r = range(b)

    # tests along main semiaxes, slightly below threshold
    @test evaluate(m, T[0,0,0], T(2.9) * T[0.5,0.5,-√2/2]) ≤ r
    @test evaluate(m, T[0,0,0], T(1.9) * T[√2/2,-√2/2,0.0]) ≤ r
    @test evaluate(m, T[0,0,0], T(0.9) * T[0.5,0.5,√2/2]) ≤ r

    # Tests along main semiaxes, slightly above threshold
    @test evaluate(m, T[0,0,0], T(3.1) * T[0.5,0.5,-√2/2]) > r
    @test evaluate(m, T[0,0,0], T(2.1) * T[√2/2,-√2/2,0.0]) > r
    @test evaluate(m, T[0,0,0], T(1.1) * T[0.5,0.5,√2/2]) > r

    # 2D simple test of default convention
    m₁ = metric(MetricBall([1.,1.], [0.]))
    m₂ = metric(MetricBall([1.,2.], [0.]))
    @test evaluate(m₁, [1.,0.], [0.,0.]) == evaluate(m₁, [0.,1.], [0.,0.])
    @test evaluate(m₂, [1.,0.], [0.,0.]) != evaluate(m₂, [0.,1.], [0.,0.])

    # 3D simple test of default convention
    m₃ = metric(MetricBall([1.,.5,.5], [π/4,0.,0.]))
    @test evaluate(m₃, [1.,1.,0.], [0.,0.,0.]) ≈ √2
    @test evaluate(m₃, [-1.,1.,0.], [0.,0.,0.]) ≈ √8

    # test of intrinsic conventions
    gslib = metric(MetricBall([50.,25.,5.], [30.,-30.,30.], convention=GSLIB))
    tait  = metric(MetricBall([25.,50.,5.], [-π/6,-π/6,π/6], convention=TaitBryanIntr))
    euler = metric(MetricBall([50.,25.,5.], deg2rad.([-78,-41,-50]), convention=EulerIntr))
    lpf   = metric(MetricBall([50.,25.,5.], [78.,41.,50.], convention=Leapfrog))
    dm    = metric(MetricBall([50.,25.,5.], [78.,41.,50.], convention=Datamine))

    @test evaluate(gslib, [1.,0.,0.], [0.,0.,0.]) ≈ 0.1325707358356285
    @test evaluate(gslib, [0.,1.,0.], [0.,0.,0.]) ≈ 0.039051248379533283
    @test evaluate(gslib, [0.,0.,1.], [0.,0.,0.]) ≈ 0.15132745950421558

    @test evaluate(gslib, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(tait, [0.,0.,1.], [0.,0.,0.])
    @test evaluate(euler, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(dm, [0.,0.,1.], [0.,0.,0.])
    @test evaluate(euler, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(lpf, [0.,0.,1.], [0.,0.,0.])
    @test evaluate(euler, [0.,0.,1.], [0.,0.,0.]) - evaluate(gslib, [0.,0.,1.], [0.,0.,0.]) < 10^-3

    # test for https://github.com/JuliaEarth/GeoStats.jl/issues/197
    gslib = metric(MetricBall([1.0,0.5,0.3], [100.,-10.,-20.], convention=GSLIB))

    @test evaluate(gslib, [1.,0.,0.], [0.,0.,0.]) ≈ 1.233956165693094
    @test evaluate(gslib, [0.,1.,0.], [0.,0.,0.]) ≈ 2.14219475359467
    @test evaluate(gslib, [0.,0.,1.], [0.,0.,0.]) ≈ 3.1621930570302914

    # test of extrinsic conventions
    xtait  = metric(MetricBall([50.,25.,5.], [π,0,π/2], convention=TaitBryanExtr))
    xeuler = metric(MetricBall([50.,25.,5.], [-π/2,-π/2,-π/2], convention=EulerExtr))

    @test evaluate(xtait, [1.,0.,0.], [0.,0.,0.]) ≈ 0.20
    @test evaluate(xtait, [0.,1.,0.], [0.,0.,0.]) ≈ 0.04
    @test evaluate(xtait, [0.,0.,1.], [0.,0.,0.]) ≈ 0.02

    @test evaluate(xtait, [1.,0.,0.], [0.,0.,0.]) ≈ evaluate(xeuler, [1.,0.,0.], [0.,0.,0.])
    @test evaluate(xtait, [0.,1.,0.], [0.,0.,0.]) ≈ evaluate(xeuler, [0.,1.,0.], [0.,0.,0.])
    @test evaluate(xtait, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(xeuler, [0.,0.,1.], [0.,0.,0.])
  end
end
