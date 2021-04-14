@testset "Topostructures" begin
  s1 = ExplicitStructure(connect.([(1,2,3),(4,3,2)], Triangle))
  @test connectivities(s1) == connect.([(1,2,3),(4,3,2)], Triangle) 

  s2 = convert(HalfEdgeStructure, s1)
  @test ncells(s2) == 2
  @test all(edgeoncell(s2, c).cell == c for c in 1:ncells(s2))
  @test edgeoncell(s2, 1).head ∈ 1:3
  @test edgeoncell(s2, 2).head ∈ 4:6

  s3 = convert(ExplicitStructure, s2)
  @test s3 == s1
end
