@testset "Discretization" begin
  @testset "FIST" begin
    p = PolySurface((0,0), (1,0), (1,0), (1,1), (1,2), (0,2), (0,1), (0,1), (0,0))
    verts, perms = Meshes._fist_remove_duplicates(p)
    @test verts == [Point.([(0,0), (1,0), (1,1), (1,2), (0,2), (0,1)])]
    @test perms == [[1,6,5,2,3,4]]
  end
end
