@testset "Neighborhoods" begin
  @testset "IsotropicBall" begin
    # Euclidean metric
    ball = IsotropicBall(T(1/2))
    @test all(evaluate(metric(ball), T[0], T[0]) .≤ radii(ball))
    @test all(evaluate(metric(ball), T[0], T[1]) .> radii(ball))

    ball = IsotropicBall(T(1))
    @test all(evaluate(metric(ball), T[0,0], T[0,0]) .≤ radii(ball))
    @test all(evaluate(metric(ball), T[0,0], T[1,0]) .≤ radii(ball))
    @test all(evaluate(metric(ball), T[0,0], T[0,1]) .≤ radii(ball))

    # Chebyshev metric
    ball = IsotropicBall(T(1/2), Chebyshev())
    @test all(evaluate(metric(ball), T[0], T[0]) .≤ radii(ball))
    @test all(evaluate(metric(ball), T[0], T[1]) .> radii(ball))

    for r in [1.,2.,3.,4.,5.]
      ball = IsotropicBall(r, Chebyshev())
      for i in 0.0:1.0:r, j in 0.0:1.0:r
        @test all(evaluate(metric(ball), T[0,0], T[i,j]) .≤ radii(ball))
      end
    end

    ball = IsotropicBall(T(1))
    if T == Float32
      @test sprint(show, ball) == "IsotropicBall{Float32}(1.0, Euclidean(0.0))"
    else
      @test sprint(show, ball) == "IsotropicBall{Float64}(1.0, Euclidean(0.0))"
    end
  end

  @testset "AnisotropicBall" begin
    # 2D ellipsoid rotated 45 degrees in GSLIB convention
    ellipsoid = AnisotropicBall(T[2,1], T[45], convention=GSLIB)

    # tests along main semiaxes, slightly below threshold
    @test evaluate(metric(ellipsoid), T[0,0], T(1.9) * T[√2/2, √2/2]) ≤ T(1)
    @test evaluate(metric(ellipsoid), T[0,0], T(0.9) * T[√2/2,-√2/2]) ≤ T(1)

    # tests along main semiaxes, slightly above threshold
    @test evaluate(metric(ellipsoid), T[0,0], T(2.1) * T[√2/2, √2/2]) > T(1)
    @test evaluate(metric(ellipsoid), T[0,0], T(1.1) * T[√2/2,-√2/2]) > T(1)

    # 3D ellipsoid rotated (45, -45, 0) in GSLIB convention
    ellipsoid = AnisotropicBall(T[3,2,1], T[45,-45,0], convention=GSLIB)

    # tests along main semiaxes, slightly below threshold
    @test evaluate(metric(ellipsoid), T[0,0,0], T(2.9) * T[0.5,0.5,-√2/2]) ≤ T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(1.9) * T[√2/2,-√2/2,0.0]) ≤ T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(0.9) * T[0.5,0.5,√2/2]) ≤ T(1)

    # Tests along main semiaxes, slightly above threshold
    @test evaluate(metric(ellipsoid), T[0,0,0], T(3.1) * T[0.5,0.5,-√2/2]) > T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(2.1) * T[√2/2,-√2/2,0.0]) > T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(1.1) * T[0.5,0.5,√2/2]) > T(1)

    # 2D simple test of default convention
    d₁ = metric(AnisotropicBall([1.,1.], [0.]))
    d₂ = metric(AnisotropicBall([1.,2.], [0.]))
    @test evaluate(d₁, [1.,0.], [0.,0.]) == evaluate(d₁, [0.,1.], [0.,0.])
    @test evaluate(d₂, [1.,0.], [0.,0.]) != evaluate(d₂, [0.,1.], [0.,0.])

    # 3D simple test of default convention
    d₃ = metric(AnisotropicBall([1.,.5,.5], [π/4,0.,0.]))
    @test evaluate(d₃, [1.,1.,0.], [0.,0.,0.]) ≈ √2
    @test evaluate(d₃, [-1.,1.,0.], [0.,0.,0.]) ≈ √8

    # test of intrinsic conventions
    gslib = metric(AnisotropicBall([50.,25.,5.], [30.,-30.,30.], convention=GSLIB))
    tait  = metric(AnisotropicBall([25.,50.,5.], [-π/6,-π/6,π/6], convention=TaitBryanIntr))
    euler = metric(AnisotropicBall([50.,25.,5.], deg2rad.([-78,-41,-50]), convention=EulerIntr))
    lpf   = metric(AnisotropicBall([50.,25.,5.], [78.,41.,50.], convention=Leapfrog))
    dm    = metric(AnisotropicBall([50.,25.,5.], [78.,41.,50.], convention=Datamine))

    @test evaluate(gslib, [1.,0.,0.], [0.,0.,0.]) ≈ 0.1325707358356285
    @test evaluate(gslib, [0.,1.,0.], [0.,0.,0.]) ≈ 0.039051248379533283
    @test evaluate(gslib, [0.,0.,1.], [0.,0.,0.]) ≈ 0.15132745950421558

    @test evaluate(gslib, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(tait, [0.,0.,1.], [0.,0.,0.])
    @test evaluate(euler, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(dm, [0.,0.,1.], [0.,0.,0.])
    @test evaluate(euler, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(lpf, [0.,0.,1.], [0.,0.,0.])
    @test evaluate(euler, [0.,0.,1.], [0.,0.,0.]) - evaluate(gslib, [0.,0.,1.], [0.,0.,0.]) < 10^-3

    # test for https://github.com/JuliaEarth/GeoStats.jl/issues/197
    gslib = metric(AnisotropicBall([1.0, 0.5, 0.3], [100., -10., -20.], convention = GSLIB))

    @test evaluate(gslib, [1.,0.,0.], [0.,0.,0.]) ≈ 1.233956165693094
    @test evaluate(gslib, [0.,1.,0.], [0.,0.,0.]) ≈ 2.14219475359467
    @test evaluate(gslib, [0.,0.,1.], [0.,0.,0.]) ≈ 3.1621930570302914

    # test of extrinsic conventions
    xtait  = metric(AnisotropicBall([50.,25.,5.], [π,0,π/2], convention=TaitBryanExtr))
    xeuler = metric(AnisotropicBall([50.,25.,5.], [-π/2,-π/2,-π/2], convention=EulerExtr))

    @test evaluate(xtait, [1.,0.,0.], [0.,0.,0.]) ≈ 0.20
    @test evaluate(xtait, [0.,1.,0.], [0.,0.,0.]) ≈ 0.04
    @test evaluate(xtait, [0.,0.,1.], [0.,0.,0.]) ≈ 0.02

    @test evaluate(xtait, [1.,0.,0.], [0.,0.,0.]) ≈ evaluate(xeuler, [1.,0.,0.], [0.,0.,0.])
    @test evaluate(xtait, [0.,1.,0.], [0.,0.,0.]) ≈ evaluate(xeuler, [0.,1.,0.], [0.,0.,0.])
    @test evaluate(xtait, [0.,0.,1.], [0.,0.,0.]) ≈ evaluate(xeuler, [0.,0.,1.], [0.,0.,0.])

    ellipsoid = AnisotropicBall(T.((2,1)), T.((45,)), convention=GSLIB)
    if T == Float32
      @test sprint(show, ellipsoid) == "AnisotropicBall{Float32}((2.0f0, 1.0f0), (45.0f0,), GSLIB)"
    else
      @test sprint(show, ellipsoid) == "AnisotropicBall{Float64}((2.0, 1.0), (45.0,), GSLIB)"
    end
  end
end
