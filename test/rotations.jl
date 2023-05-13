@testset "Rotations" begin
  cw = ClockwiseAngle(T(π / 2))
  R = convert(DCM, cw)
  @test eltype(R) == T
  @test R ≈ T[0 1; -1 0]

  ccw = CounterClockwiseAngle(T(π / 2))
  R = convert(DCM, ccw)
  @test eltype(R) == T
  @test R ≈ T[0 -1; 1 0]
end
