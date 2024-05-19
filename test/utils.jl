@testset "Utilities" begin
  a, b, c = point(0, 0), point(1, 0), point(0, 1)
  @test signarea(a, b, c) == T(0.5) * u"m^2"
  a, b, c = point(0, 0), point(0, 1), point(1, 0)
  @test signarea(a, b, c) == T(-0.5) * u"m^2"

  normals = [
    vector(1, 0, 0),
    vector(0, 1, 0),
    vector(0, 0, 1),
    vector(-1, 0, 0),
    vector(0, -1, 0),
    vector(0, 0, -1),
    vector(ntuple(i -> rand() - 0.5, 3))
  ]
  for n in normals
    u, v = householderbasis(n)
    @test u isa Vec{3}
    @test v isa Vec{3}
    @test ustrip.(u × v) ≈ n ./ norm(n)
  end

  @test Meshes.mayberound(1.1, 1.0, 0.2) ≈ 1.0
  @test Meshes.mayberound(1.1, 1.0, 0.10000000000000001) ≈ 1.1
  @test Meshes.mayberound(1.1, 1.0, 0.05) ≈ 1.1

  # intersect parameters
  p1, p2 = point(0, 0), point(1, 1)
  p3, p4 = point(1, 0), point(0, 1)
  @inferred Meshes.intersectparameters(p1, p2, p3, p4)
  @inferred Meshes.intersectparameters(p1, p3, p2, p4)
  @inferred Meshes.intersectparameters(p1, p2, p1, p2)

  p1, p2 = point(0, 0, 0), point(1, 1, 1)
  p3, p4 = point(1, 0, 0), point(0, 1, 1)
  @inferred Meshes.intersectparameters(p1, p2, p3, p4)
  @inferred Meshes.intersectparameters(p1, p3, p2, p4)
  @inferred Meshes.intersectparameters(p1, p2, p1, p2)
end
