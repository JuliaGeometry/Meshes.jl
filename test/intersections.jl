using Base: hash_integer
@testset "Intersections" begin
  @testset "Segments" begin
    # segments in 2D
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

    # segments in 3D
    s1 = Segment(P3(0.0, 0.0, 0.0), P3(1.0, 0.0, 0.0))
    s2 = Segment(P3(0.5, 1.0, 0.0), P3(0.5, -1.0, 0.0))
    s3 = Segment(P3(0.5, 0.0, 0.0), P3(1.5, 0.0, 0.0))
    s4 = Segment(P3(0.0, 1.0, 0.0), P3(0.0, -2.0, 0.0))
    s5 = Segment(P3(-1.0, 1.0, 0.0), P3(2.0, -2.0, 0.0))
    s6 = Segment(P3(0.0, 0.0, 0.0), P3(0.0, 1.0, 0.0))
    s7 = Segment(P3(-1.0, 1.0, 0.0), P3(-1.0, -1.0, 0.0))
    s8 = Segment(P3(-1.0, 1.0, 1.0), P3(-1.0, -1.0, 1.0))
    s9 = Segment(P3(0.5, 1.0, 1.0), P3(0.5, -1.0, 1.0))
    s10 = Segment(P3(0.0, 1.0, 0.0), P3(1.0, 1.0, 0.0))
    s11 = Segment(P3(1.5, 0.0, 0.0), P3(2.5, 0.0, 0.0))
    s12 = Segment(P3(1.0, 0.0, 0.0), P3(2.0, 0.0, 0.0))

    @test isa(intersecttype(s1, s2), CrossingSegments)        # CrossingSegments
    @test s1 ∩ s2 == P3(0.5, 0.0, 0.0)
    @test isa(intersecttype(s1, s3), OverlappingSegments)     # OverlappingSegments
    @test s1 ∩ s3 == Segment(P3(0.5, 0.0, 0.0), P3(1.0, 0.0, 0.0))
    @test isa(intersecttype(s1, s4), MidTouchingSegments)     # MidTouchingSegments (perpendicular)
    @test s1 ∩ s4 == P3(0.0, 0.0, 0.0)
    @test isa(intersecttype(s1, s5), MidTouchingSegments)     # MidTouchingSegments (obtuse)
    @test s1 ∩ s5 == P3(0.0, 0.0, 0.0)
    @test isa(intersecttype(s1, s6), CornerTouchingSegments)  # CornerTouchingSegments
    @test s1 ∩ s6 == P3(0.0, 0.0, 0.0)
    @test isa(intersecttype(s1, s7), NoIntersection)          # NoIntersection (but coplanar)
    @test isnothing(s1 ∩ s7)
    @test isa(intersecttype(s1, s8), NoIntersection)          # NoIntersection (non-coplanar)
    @test isnothing(s1 ∩ s8)
    @test isa(intersecttype(s1, s9), NoIntersection)          # NoIntersection (non-coplanar)
    @test isnothing(s1 ∩ s9)
    @test isa(intersecttype(s1, s10), NoIntersection)         # NoIntersection (parallel)
    @test isnothing(s1 ∩ s10)
    @test isa(intersecttype(s1, s11), NoIntersection)         # NoIntersection (colinear, not-overlapping)
    @test isnothing(s1 ∩ s11)
    @test isa(intersecttype(s1, s12), CornerTouchingSegments) # CornerTouchingSegments (colinear)
    @test s1 ∩ s12 == P3(1.0, 0.0, 0.0)
  end
    
  @testset "Triangles" begin
    # Utility to reverse segments, to more fully test branches in the intersection algorithm
    reverse_segment(s) = Segment(vertices(s)[2], vertices(s)[1])

    # segments and triangles in 3D
    t₁ = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))  # Triangle lying in XY plane
    t₂ = Triangle(P3(0, 0, 0), P3(2, 0, 0), P3(0, 2, 2))  # Inclined triangle
    
    ## Intersections with t₁
    # Intersects through t₁
    s₁ = Segment(P3(0.2, 0.2, 1.0), P3(0.2, 0.2, -1.0))
    @test isa(intersecttype(s₁, t₁), IntersectingSegmentTriangle)
    @test s₁ ∩ t₁ == P3(0.2, 0.2, 0.0)
    
    # Intersects at a vertex of t₁
    s₂ = Segment(P3(0.0, 0.0, 1.0), P3(0.0, 0.0, -1.0))
    @test isa(intersecttype(s₂, t₁), IntersectingSegmentTriangle)
    @test s₂ ∩ t₁ == P3(0.0, 0.0, 0.0)
    
    # Normal to, doesn't intersect with t₁
    s₃ = Segment(P3(0.9, 0.9, 1.0), P3(0.9, 0.9, -1.0))
    @test isa(intersecttype(s₃, t₁), NoIntersection)
    @test isnothing(s₃ ∩ t₁)
    
    # Co-planar, intersects with t₁ (but should return NoIntersection)
    s₄ = Segment(P3(-0.2, 0.2, 0.0), P3(1.2, 0.2, 0.0))
    @test isa(intersecttype(s₄, t₁), NoIntersection)
    @test isnothing(s₄ ∩ t₁)
    
    # Co-planar, doesn't intersect with t₁
    s₅ = Segment(P3(-0.2, -0.2, 0.0), P3(1.2, -0.2, 0.0))
    @test isa(intersecttype(s₅, t₁), NoIntersection)
    @test isnothing(s₅ ∩ t₁)
    
    # Parallel, above, doesn't intersect with t₁  
    s₆ = Segment(P3(-0.2, 0.2, 1.0), P3(1.2, 0.2, 1.0))
    @test isa(intersecttype(s₆, t₁), NoIntersection)
    @test isnothing(s₆ ∩ t₁)
    
    # Parallel, below, doesn't intersect with t₁  
    s₇ = Segment(P3(-0.2, 0.2, -1.0), P3(1.2, 0.2, -1.0))
    @test isa(intersecttype(s₇, t₁), NoIntersection)
    @test isnothing(s₇ ∩ t₁)
    
    # Segment colinear with edge of t₁ (but should return NoIntersection)
    s₈ = Segment(P3(-1.0, 0.0, 0.0), P3(1.0, 0.0, 0.0))
    @test isa(intersecttype(s₈, t₁), NoIntersection)
    @test isnothing(s₈ ∩ t₁)
    
    # Co-planar, within bounding box of t₁, no intersection 
    s₉ = Segment(P3(0.7, 0.8, 0.0), P3(0.8, 0.7, 0.0))
    @test isa(intersecttype(s₉, t₁), NoIntersection)
    @test isnothing(s₉ ∩ t₁)
    
    # Segment above and to right of t₁, no intersection
    s₁₀ = Segment(P3(1.0, 1.0, 0.0), P3(1.0, 1.0, 1.0))
    @test isa(intersecttype(s₁₀, t₁), NoIntersection)
    @test isnothing(s₁₀ ∩ t₁)
    
    # Segment below t₁, no intersection
    s₁₁ = Segment(P3(0.5, -1.0, 0.0), P3(0.5, -1.0, 1.0))
    @test isa(intersecttype(s₁₁, t₁), NoIntersection)
    @test isnothing(s₁₁ ∩ t₁)
    
    # Segment left of t₁, no intersection
    s₁₂ = Segment(P3(-1.0, 0.5, 0.0), P3(-1.0, 0.5, 1.0))
    @test isa(intersecttype(s₁₂, t₁), NoIntersection)
    @test isnothing(s₁₂ ∩ t₁)
    
    # Segment above and to right of t₁, no intersection
    s₁₃ = Segment(P3(1.0, 1.0, 0.0), P3(1.0, 1.0, -1.0))
    @test isa(intersecttype(s₁₃, t₁), NoIntersection)
    @test isnothing(s₁₃ ∩ t₁)
    @test isa(intersecttype(reverse_segment(s₁₃), t₁), NoIntersection)
    @test isnothing(reverse_segment(s₁₃) ∩ t₁)
    
    # Segment below t₁, no intersection
    s₁₄ = Segment(P3(0.5, -1.0, 0.0), P3(0.5, -1.0, -1.0))
    @test isa(intersecttype(s₁₄, t₁), NoIntersection)
    @test isnothing(s₁₄ ∩ t₁)
    @test isa(intersecttype(reverse_segment(s₁₄), t₁), NoIntersection)
    @test isnothing(reverse_segment(s₁₄) ∩ t₁)
    
    # Segment left of t₁, no intersection
    s₁₅ = Segment(P3(-1.0, 0.5, 0.0), P3(-1.0, 0.5, -1.0))
    @test isa(intersecttype(s₁₅, t₁), NoIntersection)
    @test isnothing(s₁₅ ∩ t₁)
    @test isa(intersecttype(reverse_segment(s₁₅), t₁), NoIntersection)
    @test isnothing(reverse_segment(s₁₅) ∩ t₁)

    # Segment above and to right of t₁, no intersection
    s₁₆ = Segment(P3(1.0, 1.0, 1.0), P3(1.0, 1.0, 0.0))
    @test isa(intersecttype(s₁₆, t₁), NoIntersection)
    @test isnothing(s₁₆ ∩ t₁)

    # Segment below t₁, no intersection
    s₁₇ = Segment(P3(0.5, -1.0, 1.0), P3(0.5, -1.0, 0.0))
    @test isa(intersecttype(s₁₇, t₁), NoIntersection)
    @test isnothing(s₁₇ ∩ t₁)

    # Segment left of t₁, no intersection
    s₁₈ = Segment(P3(-1.0, 0.5, 1.0), P3(-1.0, 0.5, 0.0))
    @test isa(intersecttype(s₁₈, t₁), NoIntersection)
    @test isnothing(s₁₈ ∩ t₁)

    ## Intersections with inclined triangle t₂
    # Doesn't reach t₂, no intersection
    s₁₉ = Segment(P3(0.5, 0.5, 1.9), P3(0.5, 0.5, 1.8))
    @test isa(intersecttype(s₁₉, t₂), NoIntersection)
    @test isnothing(s₁₉ ∩ t₂)

    # Parallel, offset from t₂, no intersection
    s₂₀ = Segment(P3(0.0, 0.5, 1.0), P3(1.0, 0.5, 1.0))
    @test isa(intersecttype(s₂₀, t₂), NoIntersection)
    @test isnothing(s₂₀ ∩ t₂)
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

  @testset "Misc" begin
    t = Triangle(P2[(0,0),(1,0),(0,1)])
    q = Quadrangle(P2[(1,1),(2,1),(2,2),(1,2)])
    @test hasintersect(t, t)
    @test hasintersect(q, q)
    @test !hasintersect(t, q)
    @test !hasintersect(q, t)

    t = Triangle(P2[(1,0),(2,0),(1,1)])
    q = Quadrangle(P2[(1.3,0.5),(2.3,0.5),(2.3,1.5),(1.3,1.5)])
    @test hasintersect(t, t)
    @test hasintersect(q, q)
    @test hasintersect(t, q)
    @test hasintersect(q, t)

    t = Triangle(P2[(0,0),(1,0),(0,1)])
    b = Ball(P2(0,0), T(1))
    @test hasintersect(t, t)
    @test hasintersect(b, b)
    @test hasintersect(t, b)
    @test hasintersect(b, t)

    t = Triangle(P2[(1,0),(2,0),(1,1)])
    b = Ball(P2(0,0), T(1))
    @test hasintersect(t, t)
    @test hasintersect(b, b)
    @test hasintersect(t, b)
    @test hasintersect(b, t)

    t = Triangle(P2[(1,0),(2,0),(1,1)])
    b = Ball(P2(-0.01,0), T(1))
    @test hasintersect(t, t)
    @test hasintersect(b, b)
    @test !hasintersect(t, b)
    @test !hasintersect(b, t)

    outer = P2[(0,0),(1,0),(1,1),(0,1),(0,0)]
    hole1 = P2[(0.2,0.2),(0.4,0.2),(0.4,0.4),(0.2,0.4),(0.2,0.2)]
    hole2 = P2[(0.6,0.2),(0.8,0.2),(0.8,0.4),(0.6,0.4),(0.6,0.2)]
    poly1 = PolyArea(outer)
    poly2 = PolyArea(outer, [hole1, hole2])
    ball1 = Ball(P2(0.5,0.5), T(0.05))
    ball2 = Ball(P2(0.3,0.3), T(0.05))
    ball3 = Ball(P2(0.7,0.3), T(0.05))
    ball4 = Ball(P2(0.3,0.3), T(0.15))
    @test hasintersect(poly1, poly1)
    @test hasintersect(poly2, poly2)
    @test hasintersect(poly1, poly2)
    @test hasintersect(poly2, poly1)
    @test hasintersect(poly1, ball1)
    @test hasintersect(poly2, ball1)
    @test hasintersect(poly1, ball2)
    @test !hasintersect(poly2, ball2)
    @test hasintersect(poly1, ball3)
    @test !hasintersect(poly2, ball3)
    @test hasintersect(poly1, ball4)
    @test hasintersect(poly2, ball4)
    mesh1 = discretize(poly1, Dehn1899())
    mesh2 = discretize(poly2, Dehn1899())
    @test hasintersect(mesh1, mesh1)
    @test hasintersect(mesh2, mesh2)
    @test hasintersect(mesh1, mesh2)
    @test hasintersect(mesh2, mesh1)
  end
end
