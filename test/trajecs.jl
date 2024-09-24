@testitem "CylindricalTrajectory" setup = [Setup] begin
  s = Segment(cart(0, 0, 0), cart(0, 0, 1))
  c = [s(t) for t in range(T(0), stop=T(1), length=10)]
  t = CylindricalTrajectory(c)
  @test crs(t) <: Cartesian{NoDatum}
  @test Meshes.lentype(t) == ℳ
  @test eltype(t) <: Cylinder
  @test nelements(t) == 10
  @test radius(t) == T(1) * u"m"
  @test topology(t) == GridTopology(10)
  @test centroid(t, 1) == cart(0, 0, 0)
  @test centroid(t, 10) == cart(0, 0, 1)

  b = BezierCurve([cart(0, 0, 0), cart(3, 3, 0), cart(3, 0, 7)])
  c = [b(t) for t in range(T(0), stop=T(1), length=20)]
  t = CylindricalTrajectory(c, T(2))
  @test crs(t) <: Cartesian{NoDatum}
  @test Meshes.lentype(t) == ℳ
  @test eltype(t) <: Cylinder
  @test nelements(t) == 20
  @test radius(t) == T(2) * u"m"
  @test topology(t) == GridTopology(20)
  @test centroid(t, 1) == cart(0, 0, 0)
  @test centroid(t, 20) == cart(3, 0, 7)

  # trajectory with single cylinder
  t = CylindricalTrajectory([cart(0, 0, 0)], T(1))
  @test crs(t) <: Cartesian{NoDatum}
  @test Meshes.lentype(t) == ℳ
  @test eltype(t) <: Cylinder
  @test nelements(t) == 1
  @test radius(t) == T(1) * u"m"
  @test topology(t) == GridTopology(1)
  @test t[1] == Cylinder(cart(0, 0, -0.5), cart(0, 0, 0.5), T(1))
end
