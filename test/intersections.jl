@testset "Intersections" begin
  # helper function for type stability tests
  function someornone(g1, g2)
    intersection(g1, g2) do I
      if type(I) == NoIntersection
        "None"
      else
        "Some"
      end
    end
  end

  @testset "Segments" begin
    # segments in 2D
    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0.5,0.0), P2(2,0))
    @test s1 ∩ s2 == Segment(P2(0.5,0.0), P2(1,0))
    @test s2 ∩ s1 == Segment(P2(0.5,0.0), P2(1,0))

    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0,0), P2(0,1))
    @test s1 ∩ s2 ≈ P2(0,0)
    @test s2 ∩ s1 ≈ P2(0,0)

    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0,0), P2(-1,0))
    @test s1 ∩ s2 ≈ P2(0,0)
    @test s2 ∩ s1 ≈ P2(0,0)

    s1 = Segment(P2(0,0), P2(0,1))
    s2 = Segment(P2(0,0), P2(0,-1))
    @test s1 ∩ s2 ≈ P2(0,0)
    @test s2 ∩ s1 ≈ P2(0,0)

    s1 = Segment(P2(1,1), P2(1,2))
    s2 = Segment(P2(1,1), P2(1,0))
    @test s1 ∩ s2 ≈ P2(1,1)
    @test s2 ∩ s1 ≈ P2(1,1)

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
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ P2(0,0)

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
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ P2(2,0)

    s1 = Segment(P2(0,0), P2(2,0))
    s2 = Segment(P2(3,0), P2(4,0))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing

    s1 = Segment(P2(2,1), P2(1,2))
    s2 = Segment(P2(1,0), P2(1,1))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing

    s1 = Segment(P2(1.5,1.5), P2(3.0,1.5))
    s2 = Segment(P2(3.0,1.0), P2(2.0,2.0))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ P2(2.5,1.5)

    s1 = Segment(P2(0.94495744, 0.53224397), P2(0.94798386, 0.5344541))
    s2 = Segment(P2(0.94798386, 0.5344541), P2(0.9472896, 0.5340202))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ P2(0.94798386, 0.5344541)

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

    @test intersection(s1, s2) |> type == CrossingSegments
    @test s1 ∩ s2 ≈ P3(0.5, 0.0, 0.0)
    @test intersection(s1, s3) |> type == OverlappingSegments
    @test s1 ∩ s3 == Segment(P3(0.5, 0.0, 0.0), P3(1.0, 0.0, 0.0))
    @test intersection(s1, s4) |> type == MidTouchingSegments
    @test s1 ∩ s4 ≈ P3(0.0, 0.0, 0.0)
    @test intersection(s1, s5) |> type == MidTouchingSegments
    @test s1 ∩ s5 ≈ P3(0.0, 0.0, 0.0)
    @test intersection(s1, s6) |> type == CornerTouchingSegments
    @test s1 ∩ s6 ≈ P3(0.0, 0.0, 0.0)
    @test intersection(s1, s7) |> type == NoIntersection
    @test isnothing(s1 ∩ s7)
    @test intersection(s1, s8) |> type == NoIntersection
    @test isnothing(s1 ∩ s8)
    @test intersection(s1, s9) |> type == NoIntersection
    @test isnothing(s1 ∩ s9)
    @test intersection(s1, s10) |> type == NoIntersection
    @test isnothing(s1 ∩ s10)
    @test intersection(s1, s11) |> type == NoIntersection
    @test isnothing(s1 ∩ s11)
    @test intersection(s1, s12) |> type == CornerTouchingSegments
    @test s1 ∩ s12 ≈ P3(1.0, 0.0, 0.0)

    # type stability tests
    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0.5,0.0), P2(2,0))
    @inferred someornone(s1, s2)

    s1 = Segment(P3(0.0, 0.0, 0.0), P3(1.0, 0.0, 0.0))
    s2 = Segment(P3(0.5, 1.0, 0.0), P3(0.5, -1.0, 0.0))
    @inferred someornone(s1, s2)
  end

  @testset "Triangles" begin
    # utility to reverse segments, to more fully
    # test branches in the intersection algorithm
    reverse_segment(s) = Segment(vertices(s)[2], vertices(s)[1])

    # intersections with triangle lying in XY plane
    t = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))

    # intersects through t
    s = Segment(P3(0.2, 0.2, 1.0), P3(0.2, 0.2, -1.0))
    @test intersection(s, t) |> type == IntersectingSegmentTriangle
    @test s ∩ t == P3(0.2, 0.2, 0.0)

    # intersects at a vertex of t
    s = Segment(P3(0.0, 0.0, 1.0), P3(0.0, 0.0, -1.0))
    @test intersection(s, t) |> type == IntersectingSegmentTriangle
    @test s ∩ t == P3(0.0, 0.0, 0.0)

    # normal to, doesn't intersect with t
    s = Segment(P3(0.9, 0.9, 1.0), P3(0.9, 0.9, -1.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # coplanar, intersects with t (but should return NoIntersection)
    s = Segment(P3(-0.2, 0.2, 0.0), P3(1.2, 0.2, 0.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # coplanar, doesn't intersect with t
    s = Segment(P3(-0.2, -0.2, 0.0), P3(1.2, -0.2, 0.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # parallel, above, doesn't intersect with t
    s = Segment(P3(-0.2, 0.2, 1.0), P3(1.2, 0.2, 1.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # parallel, below, doesn't intersect with t
    s = Segment(P3(-0.2, 0.2, -1.0), P3(1.2, 0.2, -1.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # segment colinear with edge of t (but should return NoIntersection)
    s = Segment(P3(-1.0, 0.0, 0.0), P3(1.0, 0.0, 0.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # coplanar, within bounding box of t, no intersection
    s = Segment(P3(0.7, 0.8, 0.0), P3(0.8, 0.7, 0.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # segment above and to right of t, no intersection
    s = Segment(P3(1.0, 1.0, 0.0), P3(1.0, 1.0, 1.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # segment below t, no intersection
    s = Segment(P3(0.5, -1.0, 0.0), P3(0.5, -1.0, 1.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # segment left of t, no intersection
    s = Segment(P3(-1.0, 0.5, 0.0), P3(-1.0, 0.5, 1.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # segment above and to right of t, no intersection
    s = Segment(P3(1.0, 1.0, 0.0), P3(1.0, 1.0, -1.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)
    @test intersection(reverse_segment(s), t) |> type == NoIntersection
    @test isnothing(reverse_segment(s) ∩ t)

    # segment below t, no intersection
    s = Segment(P3(0.5, -1.0, 0.0), P3(0.5, -1.0, -1.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)
    @test intersection(reverse_segment(s), t) |> type == NoIntersection
    @test isnothing(reverse_segment(s) ∩ t)

    # segment left of t, no intersection
    s = Segment(P3(-1.0, 0.5, 0.0), P3(-1.0, 0.5, -1.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)
    @test intersection(reverse_segment(s), t) |> type == NoIntersection
    @test isnothing(reverse_segment(s) ∩ t)

    # segment above and to right of t, no intersection
    s = Segment(P3(1.0, 1.0, 1.0), P3(1.0, 1.0, 0.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # segment below t, no intersection
    s = Segment(P3(0.5, -1.0, 1.0), P3(0.5, -1.0, 0.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # segment left of t, no intersection
    s = Segment(P3(-1.0, 0.5, 1.0), P3(-1.0, 0.5, 0.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # intersections with an inclined inclined triangle t
    t = Triangle(P3(0, 0, 0), P3(2, 0, 0), P3(0, 2, 2))

    # doesn't reach t, no intersection
    s = Segment(P3(0.5, 0.5, 1.9), P3(0.5, 0.5, 1.8))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # parallel, offset from t, no intersection
    s = Segment(P3(0.0, 0.5, 1.0), P3(1.0, 0.5, 1.0))
    @test intersection(s, t) |> type == NoIntersection
    @test isnothing(s ∩ t)

    # triangle as first argument
    t = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))
    s = Segment(P3(0.2, 0.2, 1.0), P3(0.2, 0.2, -1.0))
    @test intersection(t, s) |> type == IntersectingSegmentTriangle
    @test s ∩ t == t ∩ s == P3(0.2, 0.2, 0.0)

    # type stability tests
    s = Segment(P3(0.2, 0.2, 1.0), P3(0.2, 0.2, -1.0))
    t = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))
    @inferred someornone(s, t)

    # Intersection for a triangle and a ray
    t = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))

    # intersects through t
    r = Ray(P3(0.2, 0.2, 1.0), V3(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == IntersectingRayTriangle
    @test r ∩ t == P3(0.2, 0.2, 0.0)
    # Special case: the direction vector is not length enough to cross triangle
    r = Ray(P3(0.2, 0.2, 1.0), V3(0.0, 0.0, -0.00001))
    @test intersection(r, t) |> type == IntersectingRayTriangle
    @test r ∩ t ≈ P3(0.2, 0.2, 0.0)
    # Special case: reverse direction vector should not hit the triangle
    r = Ray(P3(0.2, 0.2, 1.0), V3(0.0, 0.0, 1.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # intersects at a vertex of t
    r = Ray(P3(0.0, 0.0, 1.0), V3(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == IntersectingRayTriangle
    @test r ∩ t ≈ P3(0.0, 0.0, 0.0)

    # normal to, doesn't intersect with t
    r = Ray(P3(0.9, 0.9, 1.0), V3(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # coplanar, intersects with t (but should return NoIntersection)
    r = Ray(P3(-0.2, 0.2, 0.0), V3(1.0, 0.0, 0.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # coplanar, doesn't intersect with t
    r = Ray(P3(-0.2, -0.2, 0.0), V3(1.0, 0.0, 0.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # parallel, above, doesn't intersect with t
    r = Ray(P3(-0.2, 0.2, 1.0), V3(1.0, 0.0, 0.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # parallel, below, doesn't intersect with t
    r = Ray(P3(-0.2, 0.2, -1.0), V3(1.0, 0.0, 0.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # ray colinear with edge of t (but should return NoIntersection)
    r = Ray(P3(-1.0, 0.0, 0.0), V3(1.0, 0.0, 0.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # coplanar, within bounding box of t, no intersection
    r = Ray(P3(0.7, 0.8, 0.0), V3(1.0, -1.0, 0.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # ray above and to right of t, no intersection
    r = Ray(P3(1.0, 1.0, 0.0), V3(0.0, 0.0, 1.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # ray below t, no intersection
    r = Ray(P3(0.5, -1.0, 0.0), V3(0.0, 0.0, 1.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # ray left of t, no intersection
    r = Ray(P3(-1.0, 0.5, 0.0), V3(0.0, 0.0, 1.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # ray above and to right of t, no intersection
    r = Ray(P3(1.0, 1.0, 0.0), V3(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # ray below t, no intersection
    r = Ray(P3(0.5, -1.0, 0.0), V3(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # ray left of t, no intersection
    r = Ray(P3(-1.0, 0.5, 0.0), V3(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # ray above and to right of t, no intersection
    r = Ray(P3(1.0, 1.0, 1.0), V3(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # ray below t, no intersection
    r = Ray(P3(0.5, -1.0, 1.0), V3(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # ray left of t, no intersection
    r = Ray(P3(-1.0, 0.5, 1.0), V3(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)

    # intersections with an inclined inclined triangle t
    t = Triangle(P3(0, 0, 0), P3(2, 0, 0), P3(0, 2, 2))

    # doesn't reach t, but a ray can hit the triangle
    r = Ray(P3(0.5, 0.5, 1.9), V3(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == IntersectingRayTriangle
    @test r ∩ t ≈ P3(0.5, 0.5, 0.5)

    # parallel, offset from t, no intersection
    r = Ray(P3(0.0, 0.5, 1.0), V3(1.0, 0.0, 0.0))
    @test intersection(r, t) |> type == NoIntersection
    @test isnothing(r ∩ t)
  end

  @testset "Ngons" begin
    o = Octagon([P3(0.0,0.0,1.0), P3(0.5,-0.5,0.0), P3(1.0,0.0,0.0), P3(1.5,0.5,-0.5),
                 P3(1.0,1.0,0.0), P3(0.5,1.5,0.0), P3(0.0,1.0,0.0), P3(-0.5,0.5,0.0)])

    r = Ray(P3(-1.0, -1.0, -1.0), V3(1.0, 1.0, 1.0))
    @test intersection(r, o) |> type == IntersectingRayTriangle
    @test r ∩ o ≈ P3(0.0, 0.0, 0.0)

    r = Ray(P3(-1.0, -1.0, -1.0), V3(-1.0, -1.0, -1.0))
    @test intersection(r, o) |> type == NoIntersection
    @test isnothing(r ∩ o)
  end

  @testset "Planes" begin
    p = Plane(P3(0, 0, 1), V3(1, 0, 0), V3(0, 1, 0))

    # intersecting segment and plane
    s = Segment(P3(0, 0, 0), P3(0, 2, 2))
    @test intersection(s, p) |> type == CrossingSegmentPlane
    @test s ∩ p == P3(0, 1, 1)

    # intersecting segment and plane with λ ≈ 0
    s = Segment(P3(0, 0, 1), P3(0, 2, 2))
    @test intersection(s, p) |> type == TouchingSegmentPlane
    @test s ∩ p == P3(0, 0, 1)

    # intersecting segment and plane with λ ≈ 1
    s = Segment(P3(0, 0, 2), P3(0, 2, 1))
    @test intersection(s, p) |> type == TouchingSegmentPlane
    @test s ∩ p == P3(0, 2, 1)

    # segment contained within plane
    s = Segment(P3(0, 0, 1), P3(0, -2, 1))
    @test intersection(s, p) |> type == OverlappingSegmentPlane
    @test s ∩ p == s

    # segment below plane, non-intersecting
    s = Segment(P3(0, 0, 0), P3(0, -2, -2))
    @test intersection(s, p) |> type == NoIntersection
    @test isnothing(s ∩ p)

    # segment parallel to plane, offset, non-intersecting
    s = Segment(P3(0, 0, -1), P3(0, -2, -1))
    @test intersection(s, p) |> type == NoIntersection
    @test isnothing(s ∩ p)

    # plane as first argument
    p = Plane(P3(0, 0, 1), V3(1, 0, 0), V3(0, 1, 0))
    s = Segment(P3(0, 0, 0), P3(0, 2, 2))
    @test intersection(p, s) |> type == CrossingSegmentPlane
    @test s ∩ p == p ∩ s == P3(0, 1, 1)

    # type stability tests
    s = Segment(P3(0, 0, 0), P3(0, 2, 2))
    p = Plane(P3(0, 0, 1), V3(1, 0, 0), V3(0, 1, 0))
    @inferred someornone(s, p)
  end

  @testset "Boxes" begin
    b1 = Box(P2(0,0), P2(1,1))
    b2 = Box(P2(0.5,0.5), P2(2,2))
    b3 = Box(P2(2,2), P2(3,3))
    b4 = Box(P2(1,1), P2(2,2))
    b5 = Box(P2(1.0,0.5), P2(2,2))
    @test intersection(b1, b2) |> type == OverlappingBoxes
    @test b1 ∩ b2 == Box(P2(0.5,0.5), P2(1,1))
    @test intersection(b1, b3) |> type == NoIntersection
    @test b1 ∩ b3 === nothing
    @test intersection(b1, b4) |> type == CornerTouchingBoxes
    @test b1 ∩ b4 == P2(1,1)
    @test intersection(b1, b5) |> type == FaceTouchingBoxes
    @test b1 ∩ b5 == Box(P2(1.0,0.5), P2(1,1))

    # type stability tests
    b1 = Box(P2(0,0), P2(1,1))
    b2 = Box(P2(0.5,0.5), P2(2,2))
    @inferred someornone(b1, b2)

    # Ray-Box intersection
    b = Box(P3(0,0,0), P3(1,1,1))

    r = Ray(P3(0,0,0), V3(1,1,1))
    @test intersection(r, b) |> type == CrossingRayBox
    @test r ∩ b == Segment(P3(0,0,0), P3(1,1,1))

    r = Ray(P3(-0.5,0,0), V3(1.0,1.0,1.0))
    @test intersection(r, b) |> type == CrossingRayBox
    @test r ∩ b == Segment(P3(0.0,0.5,0.5), P3(0.5,1.0,1.0))

    r = Ray(P3(3.0,0.0,0.5), V3(-1.0,1.0,0.0))
    @test intersection(r, b) |> type == NoIntersection

    r = Ray(P3(2.0,0.0,0.5), V3(-1.0,1.0,0.0))
    @test intersection(r, b) |> type == TouchingRayBox
    @test r ∩ b == P3(1.0,1.0,0.5)

    # the ray on a face of the box, got NaN in calculation
    r = Ray(P3(1.5,0.0,0.0), V3(-1.0,1.0,0.0))
    @test intersection(r, b) |> type == CrossingRayBox
    @test r ∩ b == Segment(P3(1.0,0.5,0.0), P3(0.5,1.0,0.0))
  end

  @testset "Lines" begin
    # lines in 2D
    l1 = Line(P2(0,0), P2(1,0))
    l2 = Line(P2(-1,-1), P2(-1,1))
    @test l1 ∩ l2 ≈ l2 ∩ l1 ≈ P2(-1,0)

    l1 = Line(P2(0,0), P2(1,0))
    l2 = Line(P2(0,1), P2(1,1))
    @test l1 ∩ l2 === l2 ∩ l1 === nothing

    l1 = Line(P2(0,0), P2(1,0))
    l2 = Line(P2(1,0), P2(2,0))
    @test l1 == l2
    @test l1 ∩ l2 == l2 ∩ l1 == l1

    # lines in 3D
    # not in same plane
    l1 = Line(P3(0,0,0), P3(1,0,0))
    l2 = Line(P3(1,1,1), P3(1,2,1))
    @test l1 ∩ l2 == l2 ∩ l1 === nothing
        
    # in same plane but parallel
    l1 = Line(P3(0,0,0), P3(1,0,0))
    l2 = Line(P3(0,1,1), P3(1,1,1))
    @test l1 ∩ l2 == l2 ∩ l1 === nothing

    # in same plane and colinear
    l1 = Line(P3(0,0,0), P3(1,0,0))
    l2 = Line(P3(2,0,0), P3(3,0,0))
    @test l1 ∩ l2 == l2 ∩ l1 == l1

    # crossing in one point
    l1 = Line(P3(1,2,3), P3(2,1,0))
    l2 = Line(P3(1,2,3), P3(1,1,1))
    @test l1 ∩ l2 ≈ l2 ∩ l1 ≈ P3(1,2,3)

    # type stability tests
    l1 = Line(P2(0,0), P2(1,0))
    l2 = Line(P2(-1,-1), P2(-1,1))
    @inferred someornone(l1, l2)

    p1, p2 = P2(0,0), P2(1,1)
    p3, p4 = P2(1,0), P2(0,1)
    @inferred Meshes.intersectparameters(p1, p2, p3, p4)
    @inferred Meshes.intersectparameters(p1, p3, p2, p4)
    @inferred Meshes.intersectparameters(p1, p2, p1, p2)

    p1, p2 = P3(0,0,0), P3(1,1,1)
    p3, p4 = P3(1,0,0), P3(0,1,1)
    @inferred Meshes.intersectparameters(p1, p2, p3, p4)
    @inferred Meshes.intersectparameters(p1, p3, p2, p4)
    @inferred Meshes.intersectparameters(p1, p2, p1, p2)
  end

  @testset "missing" begin
    # https://github.com/JuliaGeometry/Meshes.jl/issues/246
    b = Ball(P2(0, 0), T(1))
    t = Triangle(P2(0, 0), P2(2, 0), P2(0, 2))
    @test_throws ErrorException("not implemented") b ∩ t
    @test_throws ErrorException("not implemented") t ∩ b
  end

  @testset "hasintersect" begin
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

    # https://github.com/JuliaGeometry/Meshes.jl/issues/250
    t1 = Triangle(P3[(0,0,0), (2,0,0), (1,2,0)])
    t2 = Triangle(P3[(1,0,0), (3,0,0), (2,2,0)])
    t3 = Triangle(P3[(3,0,0), (5,0,0), (4,2,0)])
    @test hasintersect(t1, t2)
    @test hasintersect(t2, t3)
    @test !hasintersect(t1, t3)

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
