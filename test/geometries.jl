@testset "Parametric dimension" begin
  @test paramdim(Segment) == 1
  @test paramdim(Triangle) == 2
  @test paramdim(Quadrangle) == 2
  @test paramdim(PolySurface) == 2
  @test paramdim(Hexahedron) == 3
  @test paramdim(Pyramid) == 3
  @test paramdim(Tetrahedron) == 3
end