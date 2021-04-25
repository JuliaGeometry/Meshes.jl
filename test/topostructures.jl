@testset "Topostructures" begin
  s1 = ExplicitStructure(connect.([(1,2,3),(4,3,2)], Triangle))
  @test connectivities(s1) == connect.([(1,2,3),(4,3,2)], Triangle) 

  s2 = convert(HalfEdgeStructure, s1)
  @test nelements(s2) == 2
  @test all(edgeonelem(s2, e).elem == e for e in 1:nelements(s2))
  @test edgeonelem(s2, 1).head ∈ 1:3
  @test edgeonelem(s2, 2).head ∈ 4:6

  s3 = convert(ExplicitStructure, s2)
  @test s3 == s1
end
