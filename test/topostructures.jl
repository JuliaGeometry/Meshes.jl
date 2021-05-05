@testset "Topostructures" begin
  function test_halfedge(connec, structure)
    @test nelements(structure) == length(connec)
    for e in 1:nelements(structure)
      inds = collect(indices(connec[e]))
      cvec = CircularVector(inds)
      segs = [connect((cvec[i], cvec[i+1]), Segment) for i in 1:length(cvec)]
      he = edgeonelem(structure, e)
      @test he.elem == e
      @test he.head ∈ inds
      @test boundary(connec[e], 0, s2) == inds
      @test boundary(connec[e], 1, s2) == segs
      for seg in segs
        @test boundary(seg, 0, structure) == collect(indices(seg))
      end
    end
  end

  # 2 triangles
  connec = connect.([(1,2,3),(4,3,2)], Triangle)
  s1 = FullStructure(connec)
  s2 = convert(HalfEdgeStructure, s1)
  s3 = convert(FullStructure, s2)
  @test collect(elements(s1)) == connec
  @test collect(elements(s2)) == connec
  @test s3 == s1
  test_halfedge(connec, s2)
  @test adjacency(1, s2) == [2,3]
  @test adjacency(2, s2) == [4,3,1]
  @test adjacency(3, s2) == [1,2,4]
  @test adjacency(4, s2) == [3,2]
  for v in 1:4
    @test coboundary(v, 1, s2) == connect.([(v, u) for u in adjacency(v, s2)], Segment)
  end

  # 2 triangles + 2 quadrangles
  connec = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(1,5,3)], Ngon)
  s1 = FullStructure(connec)
  s2 = convert(HalfEdgeStructure, s1)
  s3 = convert(FullStructure, s2)
  @test collect(elements(s1)) == connec
  @test collect(elements(s2)) == connec
  @test s3 == s1
  test_halfedge(connec, s2)
  @test adjacency(1, s2) == [2,5,3]
  @test adjacency(2, s2) == [4,6,1]
  @test adjacency(3, s2) == [1,5,4]
  @test adjacency(4, s2) == [3,6,2]
  @test adjacency(5, s2) == [1,6,3]
  @test adjacency(6, s2) == [5,2,4]
  for v in 1:6
    @test coboundary(v, 1, s2) == connect.([(v, u) for u in adjacency(v, s2)], Segment)
  end

  # 1 triangle + 3 quadrangles + 1 hole
  connec = connect.([(1,2,6,5),(2,4,7,6),(4,3,7),(3,1,5,7)], Ngon)
  s1 = FullStructure(connec)
  s2 = convert(HalfEdgeStructure, s1)
  s3 = convert(FullStructure, s2)
  @test collect(elements(s1)) == connec
  @test collect(elements(s2)) == connec
  @test s3 == s1
  test_halfedge(connec, s2)
  @test adjacency(1, s2) == [2,5,3]
  @test adjacency(2, s2) == [4,6,1]
  @test adjacency(3, s2) == [1,7,4]
  @test adjacency(4, s2) == [3,7,2]
  @test adjacency(5, s2) == [7,1,6]
  @test adjacency(6, s2) == [5,2,7]
  @test adjacency(7, s2) == [6,4,3,5]
  for v in 1:7
    @test coboundary(v, 1, s2) == connect.([(v, u) for u in adjacency(v, s2)], Segment)
  end
end
