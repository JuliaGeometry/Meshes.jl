@testset "Topostructures" begin
  connec = connect.([(1,2,3),(4,3,2)], Triangle)
  s1 = ElementListStructure(connec)
  s2 = convert(HalfEdgeStructure, s1)
  @test nelements(s2) == 2
  for e in 1:nelements(s2)
    inds = collect(indices(connec[e]))
    cvec = CircularVector(inds)
    segs = [connect((cvec[i], cvec[i+1]), Segment) for i in 1:length(cvec)]
    he = edgeonelem(s2, e)
    @test he.elem == e
    @test he.head ∈ inds
    @test boundary(connec[e], 0, s2) == inds
    @test boundary(connec[e], 1, s2) == segs
  end
  @test boundary(connect((1,2), Segment), 0, s2) == [1,2]
  @test boundary(connect((2,3), Segment), 0, s2) == [2,3]
  @test boundary(connect((3,1), Segment), 0, s2) == [3,1]
  @test boundary(connect((4,3), Segment), 0, s2) == [4,3]
  @test boundary(connect((3,2), Segment), 0, s2) == [3,2]
  @test boundary(connect((2,4), Segment), 0, s2) == [2,4]
  @test adjacency(1, s2) == [2,3]
  @test adjacency(2, s2) == [4,3,1]
  @test adjacency(3, s2) == [1,2,4]
  @test adjacency(4, s2) == [3,2]
  s3 = convert(ElementListStructure, s2)
  @test s3 == s1

  connec = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(1,5,3)], Ngon)
  s1 = ElementListStructure(connec)
  s2 = convert(HalfEdgeStructure, s1)
  @test nelements(s2) == 4
  for e in 1:nelements(s2)
    inds = collect(indices(connec[e]))
    cvec = CircularVector(inds)
    segs = [connect((cvec[i], cvec[i+1]), Segment) for i in 1:length(cvec)]
    he = edgeonelem(s2, e)
    @test he.elem == e
    @test he.head ∈ inds
    @test boundary(connec[e], 0, s2) == inds
    @test boundary(connec[e], 1, s2) == segs
  end
  s3 = convert(ElementListStructure, s2)
  @test s3 == s1
end
