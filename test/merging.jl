@testset "Merging" begin
  s = Sphere(point(0, 0, 0), T(1))
  c = CylinderSurface(T(1))
  m = merge(s, c)
  @test m isa Multi
  @test eltype(parent(m)) <: Primitive

  s = Sphere(point(0, 0, 0), T(1))
  b = Box(point(0, 0, 0), point(1, 1, 1))
  ms = Multi([s])
  mb = Multi([b])
  @test merge(ms, b) == merge(ms, mb) == merge(s, mb)
  m = merge(ms, mb)
  @test m isa Multi
  @test eltype(parent(m)) <: Primitive

  m1 = SimpleMesh(randpoint3(3), [connect((1, 2, 3))])
  m2 = SimpleMesh(randpoint3(4), [connect((1, 2, 3, 4))])
  m = merge(m1, m2)
  @test m isa Mesh
  @test eltype(m) <: Ngon
end
