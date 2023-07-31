@testset "Predicates" begin
  @testset "issubset" begin
    point = P2(0.5, 0.5)
    box = Box(P2(0, 0), P2(1, 1))
    ball = Ball(P2(0, 0))
    tri = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    quad = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    @test point ⊆ box
    @test point ⊆ ball
    @test point ⊆ tri
    @test point ⊆ quad

    s1 = Segment(P2(0, 0), P2(1, 1))
    s2 = Segment(P2(0.5, 0.5), P2(1, 1))
    s3 = Segment(P2(0, 0), P2(0.5, 0.5))
    @test s2 ⊆ s1
    @test s3 ⊆ s1

    seg = Segment(P2(0, 0), P2(1, 1))
    box = Box(P2(0, 0), P2(1, 1))
    ball = Ball(P2(0, 0))
    @test seg ⊆ box
    @test !(seg ⊆ ball)

    t1 = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    t2 = Triangle(P2(0, 0), P2(1, 0), P2(0.8, 0.8))
    t3 = Triangle(P2(0, 0), P2(1, 0), P2(1.1, 1.1))
    @test t2 ⊆ t1
    @test !(t3 ⊆ t1)

    tri = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    box = Box(P2(0, 0), P2(1, 1))
    ball = Ball(P2(0, 0))
    quad = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    pent = Pentagon(P2(0, 0), P2(1, 0), P2(1, 1), P2(0.5, 1.5), P2(0, 1))
    poly = PolyArea(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    @test tri ⊆ quad
    @test !(quad ⊆ tri)
    @test tri ⊆ box
    @test !(box ⊆ tri)
    @test !(tri ⊆ ball)
    @test !(ball ⊆ tri)
    @test tri ⊆ pent
    @test !(pent ⊆ tri)
    @test quad ⊆ pent
    @test !(pent ⊆ quad)
    @test tri ⊆ poly
    @test !(poly ⊆ tri)
    @test quad ⊆ poly
    @test poly ⊆ quad
  end
end
