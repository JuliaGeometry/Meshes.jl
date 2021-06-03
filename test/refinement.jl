@testset "Refinement" begin
  @testset "QuadRefinement" begin
    points = P2[(0,0), (1,0), (0,1), (1,1), (0.25,0.25), (0.75,0.25), (0.5,0.75)]
    connec = connect.([(1,2,6,5),(1,5,7,3),(2,4,7,6),(3,7,4)])
    mesh   = SimpleMesh(points, connec)
    ref1   = refine(mesh, QuadRefinement())
    ref2   = refine(ref1, QuadRefinement())
    ref3   = refine(ref2, QuadRefinement())

    if visualtests
      p1 = plot(ref1, fillcolor=false)
      p2 = plot(ref2, fillcolor=false)
      p3 = plot(ref3, fillcolor=false)
      p = plot(p1, p2, p3, layout=(1,3), size=(900,300))
      @test_reference "data/quadrefine-2-$T.png" p
    end
  end

  @testset "CatmullClark" begin
    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)])
    mesh   = SimpleMesh(points, connec)
    ref1   = refine(mesh, CatmullClark())
    ref2   = refine(ref1, CatmullClark())
    ref3   = refine(ref2, CatmullClark())

    if visualtests
      p1 = plot(ref1, fillcolor=false)
      p2 = plot(ref2, fillcolor=false)
      p3 = plot(ref3, fillcolor=false)
      p = plot(p1, p2, p3, layout=(1,3), size=(900,300))
      @test_reference "data/catmullclark-1-$T.png" p
    end

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.25,0.25), (0.75,0.25), (0.5,0.75)]
    connec = connect.([(1,2,6,5),(1,5,7,3),(2,4,7,6),(3,7,4)])
    mesh   = SimpleMesh(points, connec)
    ref1   = refine(mesh, CatmullClark())
    ref2   = refine(ref1, CatmullClark())
    ref3   = refine(ref2, CatmullClark())

    if visualtests
      p1 = plot(ref1, fillcolor=false)
      p2 = plot(ref2, fillcolor=false)
      p3 = plot(ref3, fillcolor=false)
      p = plot(p1, p2, p3, layout=(1,3), size=(900,300))
      @test_reference "data/catmullclark-2-$T.png" p
    end

    points = P3[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
    connec = connect.([(1,4,3,2),(5,6,7,8),(1,2,6,5),(3,4,8,7),(1,5,8,4),(2,3,7,6)])
    mesh   = SimpleMesh(points, connec)
    ref1   = refine(mesh, CatmullClark())
    ref2   = refine(ref1, CatmullClark())
    ref3   = refine(ref2, CatmullClark())

    if visualtests
      p1 = plot(ref1, fillcolor=false)
      p2 = plot(ref2, fillcolor=false)
      p3 = plot(ref3, fillcolor=false)
      p = plot(p1, p2, p3, layout=(1,3), size=(900,300))
      @test_reference "data/catmullclark-3-$T.png" p
    end
  end
end
