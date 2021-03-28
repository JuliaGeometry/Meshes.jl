@testset "Simplification" begin
  c = Chain(P2[(0,0),(1,0),(1.5,0.5),(1,1),(0,1),(0,0)])
  s1 = simplify(c, DouglasPeucker(T(0.1)))
  s2 = simplify(c, DouglasPeucker(T(0.5)))
  @test s1 == Chain(P2[(0,0),(1,0),(1.5,0.5),(1,1),(0,1),(0,0)])
  @test s2 == Chain(P2[(0,0),(1.5,0.5),(0,1),(0,0)])
end
