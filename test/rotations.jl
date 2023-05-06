@testset "Rotations" begin
  # ClockwiseAngle
  R = Angle2d(T(-π/2))
  @test eltype(R) == T
  @test R ≈ T[0 1; -1 0]

  # CounterClockwiseAngle
  R = Angle2d(T(π/2))
  @test eltype(R) == T
  @test R ≈ T[0 -1; 1 0]
end