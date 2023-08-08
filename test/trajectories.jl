@testset "Trajectories" begin
  @testset "CylindricalTrajectory" begin
    s = Segment(P3(0, 0, 0), P3(0, 0, 1))
    c = [s(t) for t in range(T(0), stop=T(1), length=10)]
    t = CylindricalTrajectory(c)
    @test eltype(t) <: Cylinder
    @test nelements(t) == 10
    @test radius(t) == T(1)

    b = BezierCurve([P3(0, 0, 0), P3(3, 3, 0), P3(3, 0, 7)])
    c = [b(t) for t in range(T(0), stop=T(1), length=20)]
    t = CylindricalTrajectory(c, T(2))
    @test eltype(t) <: Cylinder
    @test nelements(t) == 20
    @test radius(t) == T(2)
  end
end
