@testset "Neighborhoods" begin
  @testset "NormBall" begin
    # Euclidean metric
    ball = NormBall(T(1/2))
    @test evaluate(metric(ball), T[0], T[0]) ≤ radius(ball)
    @test evaluate(metric(ball), T[0], T[1]) > radius(ball)

    ball = NormBall(T(1))
    @test evaluate(metric(ball), T[0,0], T[0,0]) ≤ radius(ball)
    @test evaluate(metric(ball), T[0,0], T[1,0]) ≤ radius(ball)
    @test evaluate(metric(ball), T[0,0], T[0,1]) ≤ radius(ball)

    # Chebyshev metric
    ball = NormBall(T(1/2), Chebyshev())
    @test evaluate(metric(ball), T[0], T[0]) ≤ radius(ball)
    @test evaluate(metric(ball), T[0], T[1]) > radius(ball)

    for r in [1.,2.,3.,4.,5.]
      ball = NormBall(r, Chebyshev())
      for i in 0.0:1.0:r, j in 0.0:1.0:r
        @test evaluate(metric(ball), T[0,0], T[i,j]) ≤ radius(ball)
      end
    end

    ball = NormBall(T(1))
    if T == Float32
      @test sprint(show, ball) == "NormBall{Float32}(1.0, Euclidean(0.0))"
    else
      @test sprint(show, ball) == "NormBall{Float64}(1.0, Euclidean(0.0))"
    end
  end

  @testset "Ellipsoid" begin
    # 2D ellipsoid rotated 45 degrees in GSLIB convention
    ellipsoid = Ellipsoid(T[2,1], T[45], convention=GSLIB)

    # tests along main semiaxes, slightly below threshold
    @test evaluate(metric(ellipsoid), T[0,0], T(1.9) * T[√2/2, √2/2]) ≤ T(1)
    @test evaluate(metric(ellipsoid), T[0,0], T(0.9) * T[√2/2,-√2/2]) ≤ T(1)

    # tests along main semiaxes, slightly above threshold
    @test evaluate(metric(ellipsoid), T[0,0], T(2.1) * T[√2/2, √2/2]) > T(1)
    @test evaluate(metric(ellipsoid), T[0,0], T(1.1) * T[√2/2,-√2/2]) > T(1)

    # 3D ellipsoid rotated (45, -45, 0) in GSLIB convention
    ellipsoid = Ellipsoid(T[3,2,1], T[45,-45,0], convention=GSLIB)

    # tests along main semiaxes, slightly below threshold
    @test evaluate(metric(ellipsoid), T[0,0,0], T(2.9) * T[0.5,0.5,-√2/2]) ≤ T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(1.9) * T[√2/2,-√2/2,0.0]) ≤ T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(0.9) * T[0.5,0.5,√2/2]) ≤ T(1)

    # Tests along main semiaxes, slightly above threshold
    @test evaluate(metric(ellipsoid), T[0,0,0], T(3.1) * T[0.5,0.5,-√2/2]) > T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(2.1) * T[√2/2,-√2/2,0.0]) > T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(1.1) * T[0.5,0.5,√2/2]) > T(1)

    ellipsoid = Ellipsoid(T.((2,1)), T.((45,)), convention=GSLIB)
    if T == Float32
      @test sprint(show, ellipsoid) == "Ellipsoid{Float32}((2.0f0, 1.0f0), (45.0f0,), GSLIB)"
    else
      @test sprint(show, ellipsoid) == "Ellipsoid{Float64}((2.0, 1.0), (45.0,), GSLIB)"
    end
  end
  
  @testset "RotationConvention" begin
    # 2d simple test
    d₁ = metric(Ellipsoid([1.,1.], [0.]))
    d₂ = metric(Ellipsoid([1.,2.], [0.]))
    @test evaluate(d₁, [1.,0.], [0.,0.]) == evaluate(d₁, [0.,1.], [0.,0.])
    @test evaluate(d₂, [1.,0.], [0.,0.]) != evaluate(d₂, [0.,1.], [0.,0.])

    # 3d simple test
    d₃ = metric(Ellipsoid([1.,.5,.5], [π/4,0.,0.]))
    @test evaluate(d₃, [1.,1.,0.], [0.,0.,0.]) ≈ √2
    @test evaluate(d₃, [-1.,1.,0.], [0.,0.,0.]) ≈ √8

    # intrinsic conventions
    gslib = metric(Ellipsoid([50.,25.,5.], [30.,-30.,30.], convention=GSLIB))
    tait  = metric(Ellipsoid([25.,50.,5.], [-π/6,-π/6,π/6], convention=TaitBryanIntr))
    euler = metric(Ellipsoid([50.,25.,5.], deg2rad.([-78,-41,-50]), convention=EulerIntr))
    lpf   = metric(Ellipsoid([50.,25.,5.], [78.,41.,50.], convention=Leapfrog))
    dm    = metric(Ellipsoid([50.,25.,5.], [78.,41.,50.], convention=Datamine))

    @test evaluate(gslib, [1.,0.,0.], [0.,0.,0.]) ≈ 0.1325707358356285
    @test evaluate(gslib, [0.,1.,0.], [0.,0.,0.]) ≈ 0.039051248379533283
    @test evaluate(gslib, [0.,0.,1.], [0.,0.,0.]) ≈ 0.15132745950421558

    @test evaluate(gslib, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(tait, [0.,0.,1.], [0.,0.,0.])
    @test evaluate(euler, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(dm, [0.,0.,1.], [0.,0.,0.])
    @test evaluate(euler, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(lpf, [0.,0.,1.], [0.,0.,0.])
    @test evaluate(euler, [0.,0.,1.], [0.,0.,0.]) - evaluate(gslib, [0.,0.,1.], [0.,0.,0.]) < 10^-3

    # extrinsic conventions
    xtait  = metric(Ellipsoid([50.,25.,5.], [π,0,π/2], convention=TaitBryanExtr))
    xeuler = metric(Ellipsoid([50.,25.,5.], [-π/2,-π/2,-π/2], convention=EulerExtr))

    @test evaluate(xtait, [1.,0.,0.], [0.,0.,0.]) ≈ 0.20
    @test evaluate(xtait, [0.,1.,0.], [0.,0.,0.]) ≈ 0.04
    @test evaluate(xtait, [0.,0.,1.], [0.,0.,0.]) ≈ 0.02

    @test evaluate(xtait, [1.,0.,0.], [0.,0.,0.]) ≈ evaluate(xeuler, [1.,0.,0.], [0.,0.,0.])
    @test evaluate(xtait, [0.,1.,0.], [0.,0.,0.]) ≈ evaluate(xeuler, [0.,1.,0.], [0.,0.,0.])
    @test evaluate(xtait, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(xeuler, [0.,0.,1.], [0.,0.,0.])
  end
end
