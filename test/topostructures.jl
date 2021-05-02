@testset "Topostructures" begin
  connec = connect.([(1,2,3),(4,3,2)], Triangle)
  s1 = ElementListStructure(connec)
  s2 = convert(HalfEdgeStructure, s1)
  @test nelements(s2) == 2
  @test all(edgeonelem(s2, e).elem == e for e in 1:nelements(s2))
  @test edgeonelem(s2, 1).head ∈ 1:3
  @test edgeonelem(s2, 2).head ∈ 4:6
  @test boundary(connec[1], 1, s2) == connect.([(1,2),(2,3),(3,1)], Segment)
  @test boundary(connec[2], 1, s2) == connect.([(4,3),(3,2),(2,4)], Segment)
  @test boundary(connec[1], 0, s2) == [1,2,3]
  @test boundary(connec[2], 0, s2) == [4,3,2]
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
  @test all(edgeonelem(s2, e).elem == e for e in 1:nelements(s2))
  @test edgeonelem(s2, 1).head ∈ (1,2,6,5)
  @test edgeonelem(s2, 2).head ∈ (2,4,6)
  @test edgeonelem(s2, 3).head ∈ (4,3,5,6)
  @test edgeonelem(s2, 4).head ∈ (1,5,3)
  @test boundary(connec[1], 1, s2) == connect.([(1,2),(2,6),(6,5),(5,1)], Segment)
  @test boundary(connec[2], 1, s2) == connect.([(2,4),(4,6),(6,2)], Segment)
  @test boundary(connec[3], 1, s2) == connect.([(4,3),(3,5),(5,6),(6,4)], Segment)
  @test boundary(connec[4], 1, s2) == connect.([(1,5),(5,3),(3,1)], Segment)
  s3 = convert(ElementListStructure, s2)
  @test s3 == s1
end
