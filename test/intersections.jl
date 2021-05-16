@testset "Intersections" begin
  @testset "Segments" begin
    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0.5,0.0), P2(2,0))
    @test s1 ∩ s2 == Segment(P2(0.5,0.0), P2(1,0))
    @test s2 ∩ s1 == Segment(P2(0.5,0.0), P2(1,0))

    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0,0), P2(0,1))
    @test s1 ∩ s2 == P2(0,0)
    @test s2 ∩ s1 == P2(0,0)

    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0,0), P2(-1,0))
    @test s1 ∩ s2 == P2(0,0)
    @test s2 ∩ s1 == P2(0,0)

    s1 = Segment(P2(0,0), P2(0,1))
    s2 = Segment(P2(0,0), P2(0,-1))
    @test s1 ∩ s2 == P2(0,0)
    @test s2 ∩ s1 == P2(0,0)

    s1 = Segment(P2(1,1), P2(1,2))
    s2 = Segment(P2(1,1), P2(1,0))
    @test s1 ∩ s2 == P2(1,1)
    @test s2 ∩ s1 == P2(1,1)

    s1 = Segment(P2(1,1), P2(2,1))
    s2 = Segment(P2(1,0), P2(3,0))
    @test s1 ∩ s2 === nothing
    @test s2 ∩ s1 === nothing

    s1 = Segment(P2(0.181429364026879, 0.546811355144474),
                  P2(0.38282226144778, 0.107781953228536))
    s2 = Segment(P2(0.412498700935005, 0.212081819871479),
                  P2(0.395936725690311, 0.252041094122474))
    @test s1 ∩ s2 === nothing
    @test s2 ∩ s1 === nothing

    s1 = Segment(P2(1,2), P2(1,0))
    s2 = Segment(P2(1,0), P2(1,1))
    @test s1 ∩ s2 == Segment(P2(1,0), P2(1,1))

    s1 = Segment(P2(0,0), P2(2,0))
    s2 = Segment(P2(-2,0), P2(-1,0))
    s3 = Segment(P2(-1,0), P2(-2,0))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing
    @test s1 ∩ s3 === s3 ∩ s1 === nothing

    s1 = Segment(P2(-1,0), P2(0,0))
    s2 = Segment(P2(0,0), P2(2,0))
    @test s1 ∩ s2 == s2 ∩ s1 == P2(0,0)

    s1 = Segment(P2(-1,0), P2(1,0))
    s2 = Segment(P2(0,0), P2(3,0))
    @test s1 ∩ s2 == s2 ∩ s1 == Segment(P2(0,0), P2(1,0))

    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0,0), P2(2,0))
    @test s1 ∩ s2 == s2 ∩ s1 == Segment(P2(0,0), P2(1,0))

    s1 = Segment(P2(0,0), P2(3,0))
    s2 = Segment(P2(1,0), P2(2,0))
    @test s1 ∩ s2 == s2 ∩ s1 == s2

    s1 = Segment(P2(0,0), P2(2,0))
    s2 = Segment(P2(1,0), P2(2,0))
    @test s1 ∩ s2 == s2 ∩ s1 == s2

    s1 = Segment(P2(0,0), P2(2,0))
    s2 = Segment(P2(1,0), P2(3,0))
    @test s1 ∩ s2 == s2 ∩ s1 == Segment(P2(1,0), P2(2,0))

    s1 = Segment(P2(0,0), P2(2,0))
    s2 = Segment(P2(2,0), P2(3,0))
    @test s1 ∩ s2 == s2 ∩ s1 == P2(2,0)

    s1 = Segment(P2(0,0), P2(2,0))
    s2 = Segment(P2(3,0), P2(4,0))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing

    s1 = Segment(P2(2,1), P2(1,2))
    s2 = Segment(P2(1,0), P2(1,1))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing
    
    s1 = Segment(P2(1.5,1.5), P2(3.0,1.5))
    s2 = Segment(P2(3.0,1.0), P2(2.0,2.0))
    @test s1 ∩ s2 == s2 ∩ s1 == P2(2.5,1.5)

    s1 = Segment(P2(0.94495744, 0.53224397), P2(0.94798386, 0.5344541))
    s2 = Segment(P2(0.94798386, 0.5344541), P2(0.9472896, 0.5340202))
    @test s1 ∩ s2 == s2 ∩ s1 == P2(0.94798386, 0.5344541) 
  end

  @testset "Lines" begin
    l1 = Line(P2(0,0), P2(1,0))
    l2 = Line(P2(-1,-1), P2(-1,1))
    @test l1 ∩ l2 == l2 ∩ l1 == P2(-1,0)

    l1 = Line(P2(0,0), P2(1,0))
    l2 = Line(P2(0,1), P2(1,1))
    @test l1 ∩ l2 === l2 ∩ l1 === nothing

    l1 = Line(P2(0,0), P2(1,0))
    l2 = Line(P2(1,0), P2(2,0))
    @test l1 == l2
    @test l1 ∩ l2 == l2 ∩ l1 == l1
  end

  @testset "Boxes" begin
    b1 = Box(P2(0,0), P2(1,1))
    b2 = Box(P2(0.5,0.5), P2(2,2))
    b3 = Box(P2(2,2), P2(3,3))
    b4 = Box(P2(1,1), P2(2,2))
    b5 = Box(P2(1.0,0.5), P2(2,2))
    @test intersecttype(b1, b2) isa OverlappingBoxes
    @test b1 ∩ b2 == Box(P2(0.5,0.5), P2(1,1))
    @test intersecttype(b1, b3) isa NoIntersection
    @test b1 ∩ b3 === nothing
    @test intersecttype(b1, b4) isa CornerTouchingBoxes
    @test b1 ∩ b4 == P2(1,1)
    @test intersecttype(b1, b5) isa FaceTouchingBoxes
    @test b1 ∩ b5 == Box(P2(1.0,0.5), P2(1,1))
  end
end
