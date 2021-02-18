@testset "Conventions" begin
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
