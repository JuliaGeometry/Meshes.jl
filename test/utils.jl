@testset "Utilities" begin
  a, b, c = P2(0, 0), P2(1, 0), P2(0, 1)
  @test signarea(a, b, c) == T(0.5)
  a, b, c = P2(0, 0), P2(0, 1), P2(1, 0)
  @test signarea(a, b, c) == T(-0.5)

  normals = [V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1), V3(-1, 0, 0), V3(0, -1, 0), V3(0, 0, -1), V3(rand(3) .- 0.5)]
  for n in normals
    u, v = householderbasis(n)
    @test u × v ≈ n ./ norm(n)
  end

  @test Meshes.mayberound(1.1, 1.0, 0.2) ≈ 1.0
  @test Meshes.mayberound(1.1, 1.0, 0.10000000000000001) ≈ 1.1
  @test Meshes.mayberound(1.1, 1.0, 0.05) ≈ 1.1

  # intersect parameters
  p1, p2 = P2(0, 0), P2(1, 1)
  p3, p4 = P2(1, 0), P2(0, 1)
  @inferred Meshes.intersectparameters(p1, p2, p3, p4)
  @inferred Meshes.intersectparameters(p1, p3, p2, p4)
  @inferred Meshes.intersectparameters(p1, p2, p1, p2)

  p1, p2 = P3(0, 0, 0), P3(1, 1, 1)
  p3, p4 = P3(1, 0, 0), P3(0, 1, 1)
  @inferred Meshes.intersectparameters(p1, p2, p3, p4)
  @inferred Meshes.intersectparameters(p1, p3, p2, p4)
  @inferred Meshes.intersectparameters(p1, p2, p1, p2)
end
