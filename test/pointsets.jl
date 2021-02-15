@testset "PointSets" begin
  pset = PointSet(rand(P2, 100))
  @test embeddim(pset) == 2
  @test coordtype(pset) == T
  @test nelements(pset) == 100

  pset = PointSet(rand(P3, 100))
  @test embeddim(pset) == 3
  @test coordtype(pset) == T
  @test nelements(pset) == 100

  pset1 = PointSet([P3(1,2,3), P3(4,5,6)])
  pset2 = PointSet(P3(1,2,3), P3(4,5,6))
  pset3 = PointSet([T.((1,2,3)), T.((4,5,6))])
  pset4 = PointSet(T.((1,2,3)), T.((4,5,6)))
  pset5 = PointSet([T[1,2,3], T[4,5,6]])
  pset6 = PointSet(T[1,2,3], T[4,5,6])
  pset7 = PointSet(T[1 4; 2 5; 3 6])
  for pset in [pset1, pset2, pset3, pset4,
               pset5, pset6, pset7]
    @test embeddim(pset) == 3
    @test coordtype(pset) == T
    @test nelements(pset) == 2
    @test pset[1] == P3(1,2,3)
    @test pset[2] == P3(4,5,6)
  end

  pset = PointSet(P2[(1,0), (0,1)])
  @test coordinates(pset, 1) == T[1, 0]
  @test coordinates(pset, 2) == T[0, 1]
  @test coordinates(pset, 1:2) == T[1 0; 0 1]
  @test nelements(pset) == 2

  pset = PointSet(P2[(1,0), (0,1)])
  @test sprint(show, pset) == "2 PointSet{2,$T}"
  if T == Float32
    @test sprint(show, MIME"text/plain"(), pset) == "2 PointSet{2,Float32}\n  └─Point(1.0f0, 0.0f0)\n  └─Point(0.0f0, 1.0f0)"
  elseif T == Float64
    @test sprint(show, MIME"text/plain"(), pset) == "2 PointSet{2,Float64}\n  └─Point(1.0, 0.0)\n  └─Point(0.0, 1.0)"
  end

  if visualtests
    Random.seed!(2021)
    @test_ref_plot "data/pset-1D-$T.png" plot(PointSet(rand(P1,10)))
    @test_ref_plot "data/pset-2D-$T.png" plot(PointSet(rand(P2,10)))
    @test_ref_plot "data/pset-3D-$T.png" plot(PointSet(rand(P3,10)))
    @test_ref_plot "data/pset-1D-$T-data.png" plot(PointSet(rand(P1,10)),1:10)
    @test_ref_plot "data/pset-2D-$T-data.png" plot(PointSet(rand(P2,10)),1:10)
    @test_ref_plot "data/pset-3D-$T-data.png" plot(PointSet(rand(P3,10)),1:10)
  end
end
