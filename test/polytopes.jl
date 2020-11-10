@testset "Polytopes" begin
  @testset "Segment" begin
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
  end
end
