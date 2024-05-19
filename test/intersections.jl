@testset "Intersections" begin
  # helper function for type stability tests
  function someornone(g1, g2)
    intersection(g1, g2) do I
      if type(I) == NotIntersecting
        "None"
      else
        "Some"
      end
    end
  end

  @testset "Points" begin
    p = point(0, 0)
    q = point(-1, -1)
    b = Box(point(0, 0), point(1, 1))
    @test p ∩ p == p
    @test q ∩ q == q
    @test p ∩ b == b ∩ p == p
    @test isnothing(p ∩ q)
    @test isnothing(q ∩ b)
  end

  @testset "Segments" begin
    # segments in 2D
    s1 = Segment(point(0, 0), point(1, 0))
    s2 = Segment(point(0.5, 0.0), point(2, 0))
    @test s1 ∩ s2 ≈ Segment(point(0.5, 0.0), point(1, 0))
    @test s2 ∩ s1 ≈ Segment(point(0.5, 0.0), point(1, 0))

    s1 = Segment(point(0, 0), point(1, -1))
    s2 = Segment(point(0.5, -0.5), point(1.5, -1.5))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ Segment(point(0.5, -0.5), point(1, -1))

    s1 = Segment(point(0, 0), point(1, 0))
    s2 = Segment(point(0, 0), point(0, 1))
    @test s1 ∩ s2 ≈ point(0, 0)
    @test s2 ∩ s1 ≈ point(0, 0)

    s1 = Segment(point(0, 0), point(1, 0))
    s2 = Segment(point(0, 0), point(-1, 0))
    @test s1 ∩ s2 ≈ point(0, 0)
    @test s2 ∩ s1 ≈ point(0, 0)

    s1 = Segment(point(0, 0), point(0, 1))
    s2 = Segment(point(0, 0), point(0, -1))
    @test s1 ∩ s2 ≈ point(0, 0)
    @test s2 ∩ s1 ≈ point(0, 0)

    s1 = Segment(point(1, 1), point(1, 2))
    s2 = Segment(point(1, 1), point(1, 0))
    @test s1 ∩ s2 ≈ point(1, 1)
    @test s2 ∩ s1 ≈ point(1, 1)

    s1 = Segment(point(1, 1), point(2, 1))
    s2 = Segment(point(1, 0), point(3, 0))
    @test s1 ∩ s2 === nothing
    @test s2 ∩ s1 === nothing

    s1 = Segment(point(0.181429364026879, 0.546811355144474), point(0.38282226144778, 0.107781953228536))
    s2 = Segment(point(0.412498700935005, 0.212081819871479), point(0.395936725690311, 0.252041094122474))
    @test s1 ∩ s2 === nothing
    @test s2 ∩ s1 === nothing

    s1 = Segment(point(1, 2), point(1, 0))
    s2 = Segment(point(1, 0), point(1, 1))
    @test s1 ∩ s2 ≈ Segment(point(1, 1), point(1, 0))
    @test s2 ∩ s1 ≈ Segment(point(1, 0), point(1, 1))

    s1 = Segment(point(0, 0), point(2, 0))
    s2 = Segment(point(-2, 0), point(-1, 0))
    s3 = Segment(point(-1, 0), point(-2, 0))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing
    @test s1 ∩ s3 === s3 ∩ s1 === nothing

    s1 = Segment(point(-1, 0), point(0, 0))
    s2 = Segment(point(0, 0), point(2, 0))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ point(0, 0)

    s1 = Segment(point(-1, 0), point(1, 0))
    s2 = Segment(point(0, 0), point(3, 0))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ Segment(point(0, 0), point(1, 0))

    s1 = Segment(point(0, 0), point(1, 0))
    s2 = Segment(point(0, 0), point(2, 0))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ Segment(point(0, 0), point(1, 0))

    s1 = Segment(point(0, 0), point(3, 0))
    s2 = Segment(point(1, 0), point(2, 0))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ s2

    s1 = Segment(point(0, 0), point(2, 0))
    s2 = Segment(point(1, 0), point(2, 0))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ s2

    s1 = Segment(point(0, 0), point(2, 0))
    s2 = Segment(point(1, 0), point(3, 0))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ Segment(point(1, 0), point(2, 0))

    s1 = Segment(point(0, 0), point(2, 0))
    s2 = Segment(point(2, 0), point(3, 0))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ point(2, 0)

    s1 = Segment(point(0, 0), point(2, 0))
    s2 = Segment(point(3, 0), point(4, 0))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing

    s1 = Segment(point(2, 1), point(1, 2))
    s2 = Segment(point(1, 0), point(1, 1))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing

    s1 = Segment(point(1.5, 1.5), point(3.0, 1.5))
    s2 = Segment(point(3.0, 1.0), point(2.0, 2.0))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ point(2.5, 1.5)

    s1 = Segment(point(0.94495744, 0.53224397), point(0.94798386, 0.5344541))
    s2 = Segment(point(0.94798386, 0.5344541), point(0.9472896, 0.5340202))
    @test s1 ∩ s2 ≈ s2 ∩ s1 ≈ point(0.94798386, 0.5344541)

    s₁ = Segment(point(0, 0), point(3, 4))
    s₂ = Segment(point(1, 2), point(3, -2))
    s₃ = Segment(point(2, 0), point(-2, 0))
    s₄ = Segment(point(0, 0), point(1, 2))
    s₅ = Segment(point(1, 2), point(3, 4))
    s₆ = Segment(point(-1, -4 / 3), point(0, 0))
    s₇ = Segment(point(1, 2), point(0, 4))
    s₈ = Segment(point(4, 16 / 3), point(3, 4))

    s₉ = Segment(point(-1, 5), point(1, 4))
    s₁₀ = Segment(point(1, 4), point(-1, 5))
    s₁₁ = Segment(point(-2, 5.5), point(-0.8, 4.9))
    s₁₂ = Segment(point(-0.8, 4.9), point(-2, 5.5))
    s₁₃ = Segment(point(-0.5, 4.75), point(0.2, 4.4))
    s₁₄ = Segment(point(0.2, 4.4), point(-0.5, 4.75))
    s₁₅ = Segment(point(0.5, 4.25), point(1, 4))
    s₁₆ = Segment(point(1, 4), point(0.5, 4.25))
    s₁₇ = Segment(point(2, 3.5), point(1.5, 3.75))
    s₁₈ = Segment(point(1.5, 3.75), point(2, 3.5))

    @test s₁ ∩ s₂ ≈ s₂ ∩ s₁ ≈ point(1.2, 1.6) # CASE 1: Crossing Segments
    @test intersection(s₁, s₂) |> type == Crossing
    @test intersection(s₂, s₁) |> type == Crossing

    @test s₁ ∩ s₃ ≈ s₃ ∩ s₁ ≈ point(0, 0) # CASE 2: EdgeTouching (s₁(0))
    @test intersection(s₁, s₃) |> type == EdgeTouching
    @test intersection(s₃, s₁) |> type == EdgeTouching

    @test s₂ ∩ s₃ ≈ s₃ ∩ s₂ ≈ point(2, 0) # CASE 2: EdgeTouching (s₃(1))
    @test intersection(s₂, s₃) |> type == EdgeTouching
    @test intersection(s₃, s₂) |> type == EdgeTouching

    @test s₁ ∩ s₄ ≈ s₄ ∩ s₁ ≈ point(0, 0) # CASE 3: CornerTouching (s₁(0), s₄(0))
    @test intersection(s₁, s₄) |> type == CornerTouching
    @test intersection(s₄, s₁) |> type == CornerTouching

    @test s₂ ∩ s₄ ≈ s₄ ∩ s₂ ≈ point(1, 2) # CASE 3: CornerTouching (s₂(0), s₄(1))
    @test intersection(s₂, s₄) |> type == CornerTouching
    @test intersection(s₄, s₂) |> type == CornerTouching

    @test s₁ ∩ s₅ ≈ s₅ ∩ s₁ ≈ point(3, 4) # CASE 3: CornerTouching (s₁(1), s₅(1))
    @test intersection(s₂, s₄) |> type == CornerTouching
    @test intersection(s₄, s₂) |> type == CornerTouching

    @test s₁ ∩ s₆ ≈ s₆ ∩ s₁ ≈ point(0, 0) # CASE 3: CornerTouching (s₁(0), s₆(1)), collinear
    @test intersection(s₁, s₆) |> type == CornerTouching
    @test intersection(s₆, s₁) |> type == CornerTouching

    @test s₂ ∩ s₇ ≈ s₇ ∩ s₂ ≈ point(1, 2) # CASE 3: CornerTouching (s₂(0), s₇(0)), collinear
    @test intersection(s₂, s₇) |> type == CornerTouching
    @test intersection(s₇, s₂) |> type == CornerTouching

    @test s₁ ∩ s₈ ≈ s₈ ∩ s₁ ≈ point(3, 4) # CASE 3: CornerTouching (s₁(1), s₈(1)), collinear
    @test intersection(s₁, s₈) |> type == CornerTouching
    @test intersection(s₈, s₁) |> type == CornerTouching

    @test s₉ ∩ s₉ ≈ s₉ # CASE 4: Overlapping (same segment)
    @test intersection(s₉, s₉) |> type == Overlapping

    @test s₉ ∩ s₁₀ ≈ s₉ # CASE 4: Overlapping (same segment, flipped points)
    @test s₁₀ ∩ s₉ ≈ s₁₀
    @test intersection(s₉, s₁₀) |> type == Overlapping
    @test intersection(s₁₀, s₉) |> type == Overlapping

    @test s₉ ∩ s₁₁ ≈ s₁₁ ∩ s₉ ≈ Segment(point(-1, 5), point(-0.8, 4.9)) # CASE 4: Overlapping (same alignment)
    @test intersection(s₉, s₁₁) |> type == Overlapping
    @test intersection(s₁₁, s₉) |> type == Overlapping

    @test s₉ ∩ s₁₂ ≈ Segment(point(-1, 5), point(-0.8, 4.9)) # CASE 4: Overlapping (opposite alignment, λ = 0 involved)
    @test s₁₂ ∩ s₉ ≈ Segment(point(-0.8, 4.9), point(-1, 5)) # flipped Points in Segment
    @test intersection(s₉, s₁₂) |> type == Overlapping
    @test intersection(s₁₂, s₉) |> type == Overlapping

    @test s₁₀ ∩ s₁₁ ≈ Segment(point(-0.8, 4.9), point(-1, 5)) # CASE 4: Overlapping (opposite alignment, λ = 1 involved)
    @test s₁₁ ∩ s₁₀ ≈ Segment(point(-1, 5), point(-0.8, 4.9)) # flipped Points in Segment
    @test intersection(s₁₀, s₁₁) |> type == Overlapping
    @test intersection(s₁₁, s₁₀) |> type == Overlapping

    @test s₉ ∩ s₁₃ ≈ s₁₃ ∩ s₉ ≈ s₁₃ # CASE 4: Overlapping (same alignment)
    @test intersection(s₉, s₁₃) |> type == Overlapping
    @test intersection(s₁₃, s₉) |> type == Overlapping

    @test s₁₄ ∩ s₉ ≈ s₁₄ # CASE 4: Overlapping (opposite alignment)
    @test s₉ ∩ s₁₄ ≈ s₁₃ # flipped Points in Segment
    @test intersection(s₉, s₁₄) |> type == Overlapping
    @test intersection(s₁₄, s₉) |> type == Overlapping

    @test s₉ ∩ s₁₅ ≈ s₁₅ ∩ s₉ ≈ s₁₅ # CASE 4: Overlapping (same alignment, corner case)
    @test intersection(s₉, s₁₅) |> type == Overlapping
    @test intersection(s₁₅, s₉) |> type == Overlapping

    @test s₁₅ ∩ s₁₀ ≈ s₁₅ # CASE 4: Overlapping (same alignment, corner case)
    @test s₁₀ ∩ s₁₅ ≈ s₁₆ # flipped Points in Segment
    @test intersection(s₁₀, s₁₅) |> type == Overlapping
    @test intersection(s₁₅, s₁₀) |> type == Overlapping

    @test s₁₆ ∩ s₉ ≈ s₁₆ # CASE 4: Overlapping (opposite alignment, corner case)
    @test s₉ ∩ s₁₆ ≈ s₁₅ # flipped Points in Segment
    @test intersection(s₉, s₁₆) |> type == Overlapping
    @test intersection(s₁₆, s₉) |> type == Overlapping

    @test s₁₀ ∩ s₁₆ ≈ s₁₆ ∩ s₁₀ ≈ s₁₆ # CASE 4: Overlapping (same alignment, corner case)
    @test intersection(s₁₀, s₁₆) |> type == Overlapping
    @test intersection(s₁₆, s₁₀) |> type == Overlapping

    @test s₉ ∩ s₁₇ === s₁₇ ∩ s₉ === nothing # CASE 5: NotIntersecting (collinear, same alignment)
    @test intersection(s₉, s₁₇) |> type == NotIntersecting
    @test intersection(s₁₇, s₉) |> type == NotIntersecting

    @test s₁₀ ∩ s₁₇ === s₁₇ ∩ s₁₀ === nothing # CASE 5: NotIntersecting (collinear, opposite alignment)
    @test intersection(s₁₀, s₁₇) |> type == NotIntersecting
    @test intersection(s₁₇, s₁₀) |> type == NotIntersecting

    @test s₉ ∩ s₁₈ === s₁₈ ∩ s₉ === nothing # CASE 5: NotIntersecting (collinear, opposite alignment)
    @test intersection(s₉, s₁₈) |> type == NotIntersecting
    @test intersection(s₁₈, s₉) |> type == NotIntersecting

    @test s₁ ∩ s₉ === s₉ ∩ s₁ === nothing # CASE 5: NotIntersecting, one λ in range
    @test intersection(s₉, s₁) |> type == NotIntersecting
    @test intersection(s₁, s₉) |> type == NotIntersecting

    @test s₁ ∩ s₁₀ === s₁₀ ∩ s₁ === nothing # CASE 5: NotIntersecting, one λ in range
    @test intersection(s₁₀, s₁) |> type == NotIntersecting
    @test intersection(s₁, s₁₀) |> type == NotIntersecting

    @test s₃ ∩ s₉ === s₉ ∩ s₃ === nothing # CASE 5: NotIntersecting
    @test intersection(s₉, s₁) |> type == NotIntersecting
    @test intersection(s₁, s₉) |> type == NotIntersecting

    @test s₃ ∩ s₁₀ === s₁₀ ∩ s₃ === nothing # CASE 5: NotIntersecting
    @test intersection(s₁₀, s₃) |> type == NotIntersecting
    @test intersection(s₃, s₁₀) |> type == NotIntersecting

    # segments in 3D
    s1 = Segment(point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0))
    s2 = Segment(point(0.5, 1.0, 0.0), point(0.5, -1.0, 0.0))
    s3 = Segment(point(0.5, 0.0, 0.0), point(1.5, 0.0, 0.0))
    s4 = Segment(point(0.0, 1.0, 0.0), point(0.0, -2.0, 0.0))
    s5 = Segment(point(-1.0, 1.0, 0.0), point(2.0, -2.0, 0.0))
    s6 = Segment(point(0.0, 0.0, 0.0), point(0.0, 1.0, 0.0))
    s7 = Segment(point(-1.0, 1.0, 0.0), point(-1.0, -1.0, 0.0))
    s8 = Segment(point(-1.0, 1.0, 1.0), point(-1.0, -1.0, 1.0))
    s9 = Segment(point(0.5, 1.0, 1.0), point(0.5, -1.0, 1.0))
    s10 = Segment(point(0.0, 1.0, 0.0), point(1.0, 1.0, 0.0))
    s11 = Segment(point(1.5, 0.0, 0.0), point(2.5, 0.0, 0.0))
    s12 = Segment(point(1.0, 0.0, 0.0), point(2.0, 0.0, 0.0))

    @test intersection(s1, s2) |> type == Crossing
    @test s1 ∩ s2 ≈ point(0.5, 0.0, 0.0)
    @test intersection(s1, s3) |> type == Overlapping
    @test s1 ∩ s3 ≈ Segment(point(0.5, 0.0, 0.0), point(1.0, 0.0, 0.0))
    @test intersection(s1, s4) |> type == EdgeTouching
    @test s1 ∩ s4 ≈ point(0.0, 0.0, 0.0)
    @test intersection(s1, s5) |> type == EdgeTouching
    @test s1 ∩ s5 ≈ point(0.0, 0.0, 0.0)
    @test intersection(s1, s6) |> type == CornerTouching
    @test s1 ∩ s6 ≈ point(0.0, 0.0, 0.0)
    @test intersection(s1, s7) |> type == NotIntersecting
    @test isnothing(s1 ∩ s7)
    @test intersection(s1, s8) |> type == NotIntersecting
    @test isnothing(s1 ∩ s8)
    @test intersection(s1, s9) |> type == NotIntersecting
    @test isnothing(s1 ∩ s9)
    @test intersection(s1, s10) |> type == NotIntersecting
    @test isnothing(s1 ∩ s10)
    @test intersection(s1, s11) |> type == NotIntersecting
    @test isnothing(s1 ∩ s11)
    @test intersection(s1, s12) |> type == CornerTouching
    @test s1 ∩ s12 ≈ point(1.0, 0.0, 0.0)

    # precision test
    s1 = Segment(point(2.0, 2.0), point(3.0, 1.0))
    s2 = Segment(point(2.12505, 1.87503), point(50000.0, 30000.0))
    s3 = Segment(point(2.125005, 1.875003), point(50000.0, 30000.0))
    s4 = Segment(point(2.125005, 1.875003), point(50002.125005, 30001.875003))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing
    @test s1 ∩ s3 === s3 ∩ s1 === ((T == Float32) ? point(2.125005, 1.875003) : nothing)
    @test s1 ∩ s4 === s4 ∩ s1 === ((T == Float32) ? point(2.125005, 1.875003) : nothing)

    # type stability tests
    s1 = Segment(point(0, 0), point(1, 0))
    s2 = Segment(point(0.5, 0.0), point(2, 0))
    @inferred someornone(s1, s2)

    s1 = Segment(point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0))
    s2 = Segment(point(0.5, 1.0, 0.0), point(0.5, -1.0, 0.0))
    @inferred someornone(s1, s2)

    # rays and segments in 2D
    r₁ = Ray(point(1, 0), vector(2, 1))
    s₁ = Segment(point(0, 2), point(2, -1)) # Crossing
    s₂ = Segment(point(0, 2), point(1, 0.5)) # NotIntersecting
    s₃ = Segment(point(0, 2), point(0.5, -0.5)) # NotIntersecting
    s₄ = Segment(point(0.5, 1), point(1.5, -1)) # EdgeTouching
    s₅ = Segment(point(1.5, 0.25), point(1.5, 2)) # EdgeTouching
    s₆ = Segment(point(1, 0), point(1, -1)) # CornerTouching
    s₇ = Segment(point(0.5, -1), point(1, 0)) # CornerTouching

    @test intersection(r₁, s₁) |> type == Crossing #CASE 1
    @test r₁ ∩ s₁ ≈ s₁ ∩ r₁ ≈ point(1.25, 0.125)
    @test intersection(r₁, s₂) |> type == NotIntersecting # CASE 5
    @test r₁ ∩ s₂ === s₂ ∩ r₁ === nothing
    @test intersection(r₁, s₃) |> type == NotIntersecting # CASE 5
    @test r₁ ∩ s₃ === s₃ ∩ r₁ === nothing
    @test intersection(r₁, s₄) |> type == EdgeTouching # CASE 2
    @test r₁ ∩ s₄ ≈ s₄ ∩ r₁ ≈ r₁(0)
    @test intersection(r₁, s₅) |> type == EdgeTouching # CASE 2
    @test r₁ ∩ s₅ ≈ s₅ ∩ r₁ ≈ point(1.5, 0.25)
    @test intersection(r₁, s₆) |> type == CornerTouching # CASE 3
    @test r₁ ∩ s₆ ≈ s₆ ∩ r₁ ≈ r₁(0)
    @test intersection(r₁, s₇) |> type == CornerTouching # CASE 3
    @test r₁ ∩ s₇ ≈ s₇ ∩ r₁ ≈ r₁(0)

    r₂ = Ray(point(3, 2), vector(1, 1))
    s₈ = Segment(point(4, 3), point(5, 4)) # Overlapping
    s₉ = Segment(point(2.5, 1.5), point(3.3, 2.3)) # Overlapping s(1)
    s₁₀ = Segment(point(3.6, 2.6), point(2.6, 1.6)) # Overlapping s(0)
    s₁₁ = Segment(point(2.2, 1.2), point(3, 2)) # CornerTouching, colinear, s(1)
    s₁₂ = Segment(point(3, 2), point(2.4, 1.4)) # CornerTouching, colinear, s(0)
    s₁₃ = Segment(point(3, 2), point(3.1, 2.1)) # Overlapping s(0) = r(0)
    s₁₄ = Segment(point(3.2, 2.2), point(3, 2)) # Overlapping s(1) = r(0)
    s₁₅ = Segment(point(2, 1), point(1.6, 0.6)) # No Intersection, colinear
    s₁₆ = Segment(point(3, 1), point(4, 2)) # No Intersection, parallel
    @test intersection(r₂, s₈) |> type == Overlapping # CASE 4
    @test r₂ ∩ s₈ === s₈ ∩ r₂ === s₈
    @test intersection(r₂, s₉) |> type == Overlapping # CASE 4
    @test r₂ ∩ s₉ == s₉ ∩ r₂ == Segment(r₂(0), s₉(1))
    @test intersection(r₂, s₁₀) |> type == Overlapping # CASE 4
    @test r₂ ∩ s₁₀ == s₁₀ ∩ r₂ == Segment(r₂(0), s₁₀(0))
    @test intersection(r₂, s₁₁) |> type == CornerTouching # CASE 3
    @test r₂ ∩ s₁₁ ≈ s₁₁ ∩ r₂ ≈ r₂(0)
    @test intersection(r₂, s₁₂) |> type == CornerTouching # CASE 3
    @test r₂ ∩ s₁₂ ≈ s₁₂ ∩ r₂ ≈ r₂(0)
    @test intersection(r₂, s₁₃) |> type == Overlapping # CASE 4
    @test r₂ ∩ s₁₃ === s₁₃ ∩ r₂ === s₁₃
    @test intersection(r₂, s₁₄) |> type == Overlapping # CASE 4
    @test r₂ ∩ s₁₄ === s₁₄ ∩ r₂ === s₁₄
    @test intersection(r₂, s₁₅) |> type == NotIntersecting # CASE 5
    @test r₂ ∩ s₁₅ === s₁₅ ∩ r₂ === nothing
    @test intersection(r₂, s₁₆) |> type == NotIntersecting # CASE 5
    @test r₂ ∩ s₁₆ === s₁₆ ∩ r₂ === nothing

    # type stability tests
    r₁ = Ray(point(0, 0), vector(1, 0))
    s₁ = Segment(point(-1, -1), point(-1, 1))
    @inferred someornone(r₁, s₁)

    # 3D test
    r₁ = Ray(point(1, 2, 3), vector(1, 2, 3))
    s₁ = Segment(point(1, 3, 5), point(3, 5, 7))
    @test intersection(r₁, s₁) |> type === Crossing # CASE 1
    @test r₁ ∩ s₁ ≈ s₁ ∩ r₁ ≈ point(2, 4, 6)

    s₂ = Segment(point(0, 1, 2), point(2, 3, 4))
    @test intersection(r₁, s₂) |> type === EdgeTouching # CASE 2
    @test r₁ ∩ s₂ == s₂ ∩ r₁ == r₁(0)

    s₃ = Segment(point(0.23, 1, 2.3), point(1, 2, 3))
    @test intersection(r₁, s₃) |> type === CornerTouching # CASE 3
    @test r₁ ∩ s₃ == s₃ ∩ r₁ == r₁(0)

    s₄ = Segment(point(0, 0, 0), point(2, 4, 6))
    @test intersection(r₁, s₄) |> type === Overlapping # CASE 4
    @test r₁ ∩ s₄ == s₄ ∩ r₁ == Segment(point(1, 2, 3), point(2, 4, 6))

    s₅ = Segment(point(0, 0, 0), point(0.5, 1, 1.5))
    @test intersection(r₁, s₅) |> type === NotIntersecting # CASE 5
    @test r₁ ∩ s₅ === s₅ ∩ r₁ === nothing

    l₁ = Line(point(1, 0), point(3, 1))
    s₁ = Segment(point(0, 2), point(2, -1)) # Crossing
    s₂ = Segment(point(0.5, 1), point(0, 0)) # NotIntersecting
    s₃ = Segment(point(0, 2), point(-2, 1)) # NotIntersecting
    s₄ = Segment(point(0.5, -1), point(1, 0)) # Touching
    s₅ = Segment(point(1.5, 0.25), point(1.5, 2)) # Touching
    s₆ = Segment(point(-3, -2), point(4, 1.5)) # Overlapping

    @test intersection(l₁, s₁) |> type == Crossing #CASE 1
    @test l₁ ∩ s₁ ≈ s₁ ∩ l₁ ≈ point(1.25, 0.125)
    @test intersection(l₁, s₂) |> type == NotIntersecting # CASE 4
    @test l₁ ∩ s₂ === s₂ ∩ l₁ === nothing
    @test intersection(l₁, s₃) |> type == NotIntersecting # CASE 4
    @test l₁ ∩ s₃ === s₃ ∩ l₁ === nothing
    @test intersection(l₁, s₄) |> type == Touching # CASE 2
    @test l₁ ∩ s₄ ≈ s₄ ∩ l₁ ≈ s₄(1)
    @test intersection(l₁, s₅) |> type == Touching # CASE 2
    @test l₁ ∩ s₅ ≈ s₅ ∩ l₁ ≈ s₅(0)
    @test intersection(l₁, s₆) |> type == Overlapping # CASE 3
    @test l₁ ∩ s₆ ≈ s₆ ∩ l₁ ≈ s₆

    # type stability tests
    @inferred someornone(l₁, s₁)
    @inferred someornone(l₁, s₂)

    # 3d tests
    l₁ = Line(point(1, 0, 1), point(3, 1, 1))
    s₁ = Segment(point(0, 2, 1), point(2, -1, 1)) # Crossing
    s₂ = Segment(point(0.5, 1, 1), point(0, 0, 1)) # NotIntersecting
    s₃ = Segment(point(0, 2, 1), point(-2, 1, 1)) # NotIntersecting
    s₄ = Segment(point(0.5, -1, 1), point(1, 0, 1)) # Touching
    s₅ = Segment(point(1.5, 0.25, 1), point(1.5, 2, 1)) # Touching
    s₆ = Segment(point(-3, -2, 1), point(4, 1.5, 1)) # Overlapping
    s₇ = Segment(point(0, 2, 1), point(2, -1, 1.1)) # NotIntersecting

    @test intersection(l₁, s₁) |> type == Crossing #CASE 1
    @test l₁ ∩ s₁ ≈ s₁ ∩ l₁ ≈ point(1.25, 0.125, 1)
    @test intersection(l₁, s₂) |> type == NotIntersecting # CASE 4
    @test l₁ ∩ s₂ === s₂ ∩ l₁ === nothing
    @test intersection(l₁, s₃) |> type == NotIntersecting # CASE 4
    @test l₁ ∩ s₃ === s₃ ∩ l₁ === nothing
    @test intersection(l₁, s₄) |> type == Touching # CASE 2
    @test l₁ ∩ s₄ ≈ s₄ ∩ l₁ ≈ s₄(1)
    @test intersection(l₁, s₅) |> type == Touching # CASE 2
    @test l₁ ∩ s₅ ≈ s₅ ∩ l₁ ≈ s₅(0)
    @test intersection(l₁, s₆) |> type == Overlapping # CASE 3
    @test l₁ ∩ s₆ ≈ s₆ ∩ l₁ ≈ s₆
    @test intersection(l₁, s₇) |> type == NotIntersecting # CASE 4
    @test l₁ ∩ s₇ === s₇ ∩ l₁ === nothing

    # degenerate segments
    A = point(0.0, 0.0)
    B = point(0.5, 0.0)
    C = point(1.0, 0.0)
    s₀ = Segment(A, C)
    s₁ = Segment(A, A)
    s₂ = Segment(B, B)
    s₃ = Segment(C, C)
    @test s₀ ∩ s₁ ≈ s₁ ∩ s₀ ≈ A
    @test s₀ ∩ s₂ ≈ s₂ ∩ s₀ ≈ B
    @test s₀ ∩ s₃ ≈ s₃ ∩ s₀ ≈ C
    @test intersection(s₀, s₁) |> type == CornerTouching
    @test intersection(s₀, s₂) |> type == EdgeTouching
    @test intersection(s₀, s₃) |> type == CornerTouching
    @test s₁ ∩ s₂ === s₂ ∩ s₁ === nothing
    @test s₁ ∩ s₃ === s₃ ∩ s₁ === nothing
    @test s₂ ∩ s₃ === s₃ ∩ s₂ === nothing
    @test intersection(s₁, s₂) |> type == NotIntersecting
    @test intersection(s₁, s₃) |> type == NotIntersecting
    @test intersection(s₂, s₃) |> type == NotIntersecting
    @test s₁ ∩ s₁ ≈ A
    @test s₂ ∩ s₂ ≈ B
    @test s₃ ∩ s₃ ≈ C
    @test intersection(s₁, s₁) |> type == CornerTouching
    @test intersection(s₂, s₂) |> type == CornerTouching
    @test intersection(s₃, s₃) |> type == CornerTouching

    # utils
    @test Meshes._sort4vals(2.5, 1.4, 1.1, 2.0) == (1.4, 2.0)
    @test Meshes._sort4vals(2.0, 1.1, 1.4, 2.5) == (1.4, 2.0)
    @test Meshes._sort4vals(2.0, 2.5, 1.1, 1.4) == (1.4, 2.0)
  end

  @testset "Rays" begin
    # rays in 2D
    r₁ = Ray(point(1, 0), vector(2, 1))
    r₂ = Ray(point(0, 2), vector(2, -3))
    r₃ = Ray(point(0.5, 1), vector(1, -2))
    r₄ = Ray(point(0, 2), vector(1, -3))
    r₅ = Ray(point(4, 1.5), vector(4, 2))
    r₆ = Ray(point(2, 0.5), vector(-0.5, -0.25))
    r₇ = Ray(point(4, 0), vector(0, 1))
    @test intersection(r₁, r₂) |> type == Crossing #CASE 1
    @test r₁ ∩ r₂ ≈ point(1.25, 0.125)
    @test r₁ ∩ r₇ ≈ point(4, 1.5)
    @test intersection(r₁, r₃) |> type == EdgeTouching #CASE 2
    @test r₁ ∩ r₃ ≈ r₁(0) # origin of first ray
    @test r₅ ∩ r₇ ≈ r₅(0)
    @test intersection(r₃, r₁) |> type == EdgeTouching
    @test r₃ ∩ r₁ ≈ r₁(0) # origin of second ray
    @test r₇ ∩ r₅ ≈ r₅(0)
    @test intersection(r₂, r₄) |> type == CornerTouching #CASE 3
    @test r₂ ∩ r₄ ≈ r₂(0) ≈ r₄(0)
    @test intersection(r₅, r₁) |> type == PosOverlapping #CASE 4
    @test r₅ ∩ r₁ == r₅ # first ray
    @test intersection(r₁, r₅) |> type == PosOverlapping #CASE 4
    @test r₁ ∩ r₅ == r₅ # second ray
    @test intersection(r₁, r₆) |> type == NegOverlapping #CASE 5
    @test r₁ ∩ r₆ == Segment(r₁(0), r₆(0))
    @test intersection(r₁, r₄) |> type == NotIntersecting #CASE 6
    @test r₁ ∩ r₄ === r₄ ∩ r₁ === nothing

    # lines and rays in 2D
    l₁ = Line(point(0, 0), point(4, 5))
    r₁ = Ray(point(3, 4), vector(1, -2)) # crossing ray
    r₂ = Ray(point(1, 1.25), vector(1, 0.3)) # touching ray
    r₃ = Ray(point(-1, -1.25), vector(-1, -1.25)) # overlapping ray
    r₄ = Ray(point(1, 3), vector(1, 1.25)) # parallel ray
    r₅ = Ray(point(1, 1), vector(1, -1)) # no Intersection

    @test l₁ ∩ r₁ ≈ r₁ ∩ l₁ ≈ point(3.0769230769230766, 3.846153846153846) # CASE 1
    @test intersection(l₁, r₁) |> type === Crossing

    @test l₁ ∩ r₂ == r₂ ∩ l₁ == r₂(0) # CASE 2
    @test intersection(l₁, r₂) |> type === Touching

    @test l₁ ∩ r₃ == r₃ ∩ l₁ == r₃ # CASE 3
    @test intersection(l₁, r₃) |> type === Overlapping

    @test l₁ ∩ r₄ == r₄ ∩ l₁ === nothing # CASE 4 parallel
    @test intersection(l₁, r₄) |> type === NotIntersecting

    @test l₁ ∩ r₅ == r₅ ∩ l₁ === nothing # CASE 4 no intersection
    @test intersection(l₁, r₅) |> type === NotIntersecting

    # type stability tests
    @inferred someornone(l₁, r₁)
    @inferred someornone(l₁, r₅)

    # 3D tests
    # lines and rays in 3D
    l₁ = Line(point(0, 0, 0.1), point(4, 5, 0.1))
    r₁ = Ray(point(3, 4, 0.1), vector(1, -2, 0)) # crossing ray
    r₂ = Ray(point(1, 1.25, 0.1), vector(1, 0.3, 0)) # touching ray
    r₃ = Ray(point(-1, -1.25, 0.1), vector(-1, -1.25, 0)) # overlapping ray
    r₄ = Ray(point(1, 3, 0.1), vector(1, 1.25, 0)) # parallel ray
    r₅ = Ray(point(1, 1, 0.1), vector(1, -1, 0)) # no Intersection
    r₆ = Ray(point(3, 4, 0), vector(1, -2, 1)) # crossing ray

    @test l₁ ∩ r₁ ≈ r₁ ∩ l₁ ≈ point(3.0769230769230766, 3.846153846153846, 0.1) # CASE 1
    @test intersection(l₁, r₁) |> type === Crossing

    @test l₁ ∩ r₂ == r₂ ∩ l₁ == r₂(0) # CASE 2
    @test intersection(l₁, r₂) |> type === Touching

    @test l₁ ∩ r₃ == r₃ ∩ l₁ == r₃ # CASE 3
    @test intersection(l₁, r₃) |> type === Overlapping

    @test l₁ ∩ r₄ == r₄ ∩ l₁ === nothing # CASE 4 parallel
    @test intersection(l₁, r₄) |> type === NotIntersecting

    @test l₁ ∩ r₅ == r₅ ∩ l₁ === nothing # CASE 4 no intersection
    @test intersection(l₁, r₅) |> type === NotIntersecting

    @test l₁ ∩ r₆ == r₆ ∩ l₁ === nothing # CASE 4 no intersection
    @test intersection(l₁, r₆) |> type === NotIntersecting
  end

  @testset "Lines" begin
    # lines in 2D
    l1 = Line(point(0, 0), point(1, 0))
    l2 = Line(point(-1, -1), point(-1, 1))
    @test l1 ∩ l2 ≈ l2 ∩ l1 ≈ point(-1, 0)

    l1 = Line(point(0, 0), point(1, 0))
    l2 = Line(point(0, 1), point(1, 1))
    @test l1 ∩ l2 === l2 ∩ l1 === nothing

    l1 = Line(point(0, 0), point(1, 0))
    l2 = Line(point(1, 0), point(2, 0))
    @test l1 == l2
    @test l1 ∩ l2 == l2 ∩ l1 == l1

    # rounding errors
    l1 = Line(point(3.0, 1.0), point(2.0, 2.0))
    for k in 1:1000
      Δ = k * atol(T)
      l2 = Line(point(1.5, 1.5 + Δ), point(3.0, 1.5 + Δ))
      p = point(2.5 - Δ, 1.5 + Δ)
      @test l1 ∩ l2 ≈ l2 ∩ l1 ≈ p
    end

    # lines in 3D
    # not in same plane
    l1 = Line(point(0, 0, 0), point(1, 0, 0))
    l2 = Line(point(1, 1, 1), point(1, 2, 1))
    @test l1 ∩ l2 == l2 ∩ l1 === nothing

    # in same plane but parallel
    l1 = Line(point(0, 0, 0), point(1, 0, 0))
    l2 = Line(point(0, 1, 1), point(1, 1, 1))
    @test l1 ∩ l2 == l2 ∩ l1 === nothing

    # in same plane and colinear
    l1 = Line(point(0, 0, 0), point(1, 0, 0))
    l2 = Line(point(2, 0, 0), point(3, 0, 0))
    @test l1 ∩ l2 == l2 ∩ l1 == l1

    # crossing in one point
    l1 = Line(point(1, 2, 3), point(2, 1, 0))
    l2 = Line(point(1, 2, 3), point(1, 1, 1))
    @test l1 ∩ l2 ≈ l2 ∩ l1 ≈ point(1, 2, 3)

    # type stability tests
    l1 = Line(point(0, 0), point(1, 0))
    l2 = Line(point(-1, -1), point(-1, 1))
    @inferred someornone(l1, l2)
  end

  @testset "Chains" begin
    # https://github.com/JuliaGeometry/Meshes.jl/issues/644
    r = Rope(point(0, 0), point(1, 1))
    @test r ∩ r == GeometrySet([Segment(point(0, 0), point(1, 1))])
    @inferred someornone(r, r)
  end

  @testset "Planes" begin
    # ---------
    # SEGMENTS
    # ---------

    p = Plane(point(0, 0, 1), vector(1, 0, 0), vector(0, 1, 0))

    # intersecting segment and plane
    s = Segment(point(0, 0, 0), point(0, 2, 2))
    @test intersection(s, p) |> type == Crossing
    @test s ∩ p == point(0, 1, 1)

    # intersecting segment and plane with λ ≈ 0
    s = Segment(point(0, 0, 1), point(0, 2, 2))
    @test intersection(s, p) |> type == Touching
    @test s ∩ p == point(0, 0, 1)

    # intersecting segment and plane with λ ≈ 1
    s = Segment(point(0, 0, 2), point(0, 2, 1))
    @test intersection(s, p) |> type == Touching
    @test s ∩ p == point(0, 2, 1)

    # segment contained within plane
    s = Segment(point(0, 0, 1), point(0, -2, 1))
    @test intersection(s, p) |> type == Overlapping
    @test s ∩ p == s

    # segment below plane, non-intersecting
    s = Segment(point(0, 0, 0), point(0, -2, -2))
    @test intersection(s, p) |> type == NotIntersecting
    @test isnothing(s ∩ p)

    # segment parallel to plane, offset, non-intersecting
    s = Segment(point(0, 0, -1), point(0, -2, -1))
    @test intersection(s, p) |> type == NotIntersecting
    @test isnothing(s ∩ p)

    # plane as first argument
    p = Plane(point(0, 0, 1), vector(1, 0, 0), vector(0, 1, 0))
    s = Segment(point(0, 0, 0), point(0, 2, 2))
    @test intersection(p, s) |> type == Crossing
    @test s ∩ p == p ∩ s == point(0, 1, 1)

    # type stability tests
    s = Segment(point(0, 0, 0), point(0, 2, 2))
    p = Plane(point(0, 0, 1), vector(1, 0, 0), vector(0, 1, 0))
    @inferred someornone(s, p)

    # -----
    # RAYS
    # -----

    p = Plane(point(0, 0, 1), vector(1, 0, 0), vector(0, 1, 0))

    # intersecting ray and plane
    r = Ray(point(0, 0, 0), vector(0, 2, 2))
    @test intersection(r, p) |> type == Crossing
    @test r ∩ p == point(0, 1, 1)

    # intersecting ray and plane with λ ≈ 0
    r = Ray(point(0, 0, 1), vector(0, 2, 1))
    @test intersection(r, p) |> type == Touching
    @test r ∩ p == point(0, 0, 1)

    # intersecting ray and plane with λ ≈ 1 (only case where Ray different to Segment)
    r = Ray(point(0, 0, 2), vector(0, 2, -1))
    @test intersection(r, p) |> type == Crossing
    @test r ∩ p == point(0, 2, 1)

    # ray contained within plane
    r = Ray(point(0, 0, 1), vector(0, -2, 0))
    @test intersection(r, p) |> type == Overlapping
    @test r ∩ p == r

    # ray below plane, non-intersecting
    r = Ray(point(0, 0, 0), vector(0, -2, -2))
    @test intersection(r, p) |> type == NotIntersecting
    @test isnothing(r ∩ p)

    # ray parallel to plane, offset, non-intersecting
    r = Ray(point(0, 0, -1), vector(0, -2, 0))
    @test intersection(r, p) |> type == NotIntersecting
    @test isnothing(r ∩ p)

    # plane as first argument
    p = Plane(point(0, 0, 1), vector(1, 0, 0), vector(0, 1, 0))
    r = Ray(point(0, 0, 0), vector(0, 2, 2))
    @test intersection(p, r) |> type == Crossing
    @test r ∩ p == p ∩ r == point(0, 1, 1)

    # ------
    # LINES
    # ------

    p = Plane(point(0, 0, 1), vector(1, 0, 0), vector(0, 1, 0))

    # intersecting line and plane
    l = Line(point(0, 0, 0), point(0, 2, 2))
    @test intersection(l, p) |> type == Crossing
    @test l ∩ p == point(0, 1, 1)

    # intersecting line and plane with λ ≈ 0
    l = Line(point(0, 0, 1), point(0, 2, 2))
    @test intersection(l, p) |> type == Crossing
    @test l ∩ p == point(0, 0, 1)

    # intersecting line and plane with λ ≈ 1
    l = Line(point(0, 0, 2), point(0, 2, 1))
    @test intersection(l, p) |> type == Crossing
    @test l ∩ p == point(0, 2, 1)

    # line contained within plane
    l = Line(point(0, 0, 1), point(0, -2, 1))
    @test intersection(l, p) |> type == Overlapping
    @test l ∩ p == l

    # line below plane, non-intersecting
    l = Line(point(0, 0, 0), point(0, -2, -2))
    @test intersection(l, p) |> type == Crossing
    @test l ∩ p == point(0, 1, 1)

    # line parallel to plane, offset, non-intersecting
    l = Line(point(0, 0, -1), point(0, -2, -1))
    @test intersection(l, p) |> type == NotIntersecting
    @test isnothing(l ∩ p)

    # plane as first argument
    p = Plane(point(0, 0, 1), vector(1, 0, 0), vector(0, 1, 0))
    l = Line(point(0, 0, 0), point(0, 2, 2))
    @test intersection(p, l) |> type == Crossing
    @test l ∩ p == p ∩ l == point(0, 1, 1)

    # type stability tests
    l = Line(point(0, 0, 0), point(0, 2, 2))
    p = Plane(point(0, 0, 1), vector(1, 0, 0), vector(0, 1, 0))
    @inferred someornone(l, p)

    # ------
    # PLANES
    # ------

    p1 = Plane(point(0, 0, 0), vector(0, 0, 1))

    # p1 parallel to p2
    p2 = Plane(point(0, 0, 1), vector(0, 0, 1))
    @test intersection(p1, p2) |> type == NotIntersecting
    @test isnothing(p1 ∩ p2)

    # p1 intersects p2
    p2 = Plane(point(0, 0, 1), vector(1 / sqrt(2), 0, 1 / sqrt(2)))
    @test intersection(p1, p2) |> type == Intersecting
    @test p1 ∩ p2 == Line(point(1, 0, 0), point(1, 1, 0))
  end

  @testset "Boxes" begin
    b1 = Box(point(0, 0), point(1, 1))
    b2 = Box(point(0.5, 0.5), point(2, 2))
    b3 = Box(point(2, 2), point(3, 3))
    b4 = Box(point(1, 1), point(2, 2))
    b5 = Box(point(1.0, 0.5), point(2, 2))
    b6 = Box(point(0, 2), point(1, 3))
    b7 = Box(point(0, 1), point(1, 2))
    b8 = Box(point(0, -1), point(1, 0))
    b9 = Box(point(1, 0), point(2, 1))
    b10 = Box(point(-1, 0), point(0, 1))
    @test intersection(b1, b2) |> type == Overlapping
    @test b1 ∩ b2 == Box(point(0.5, 0.5), point(1, 1))
    @test intersection(b1, b3) |> type == NotIntersecting
    @test isnothing(b1 ∩ b3)
    @test intersection(b1, b4) |> type == CornerTouching
    @test b1 ∩ b4 == point(1, 1)
    @test intersection(b1, b5) |> type == Touching
    @test b1 ∩ b5 == Box(point(1.0, 0.5), point(1, 1))
    @test intersection(b1, b6) |> type == NotIntersecting
    @test isnothing(b1 ∩ b6)
    @test intersection(b1, b7) |> type == Touching
    @test b1 ∩ b7 == Box(point(0, 1), point(1, 1))
    @test intersection(b1, b8) |> type == Touching
    @test b1 ∩ b8 == Box(point(0, 0), point(1, 0))
    @test intersection(b1, b9) |> type == Touching
    @test b1 ∩ b9 == Box(point(1, 0), point(1, 1))
    @test intersection(b1, b10) |> type == Touching
    @test b1 ∩ b10 == Box(point(0, 0), point(0, 1))

    # more touching examples
    b1 = Box(point(0, 0), point(1, 1))
    b2 = Box(point(1.0, 0.5), point(2, 1))
    b3 = Box(point(-1, 0), point(0.0, 0.5))
    b4 = Box(point(0, 1), point(0.5, 2.0))
    b5 = Box(point(0.5, -1.0), point(1, 0))
    @test intersection(b1, b2) |> type == Touching
    @test b1 ∩ b2 == Box(point(1.0, 0.5), point(1, 1))
    @test intersection(b1, b3) |> type == Touching
    @test b1 ∩ b3 == Box(point(0.0, 0.0), point(0.0, 0.5))
    @test intersection(b1, b4) |> type == Touching
    @test b1 ∩ b4 == Box(point(0.0, 1.0), point(0.5, 1.0))
    @test intersection(b1, b5) |> type == Touching
    @test b1 ∩ b5 == Box(point(0.5, 0.0), point(1.0, 0.0))

    # tricky examples with degenerate boxes
    b1 = Box(point(0, 0, 0), point(2, 2, 0))
    b2 = Box(point(3, 0, 0), point(5, 2, 0))
    b3 = Box(point(1, 0, 0), point(3, 2, 0))
    @test intersection(b1, b2) |> type == NotIntersecting
    @test isnothing(b1 ∩ b2)
    @test intersection(b1, b3) |> type == Touching
    @test b1 ∩ b3 == Box(point(1, 0, 0), point(2, 2, 0))

    # type stability tests
    b1 = Box(point(0, 0), point(1, 1))
    b2 = Box(point(0.5, 0.5), point(2, 2))
    @inferred someornone(b1, b2)

    # Ray-Box intersection
    b = Box(point(0, 0, 0), point(1, 1, 1))

    r = Ray(point(0, 0, 0), vector(1, 1, 1))
    @test intersection(r, b) |> type == Crossing
    @test r ∩ b == Segment(point(0, 0, 0), point(1, 1, 1))

    r = Ray(point(-0.5, 0, 0), vector(1.0, 1.0, 1.0))
    @test intersection(r, b) |> type == Crossing
    @test r ∩ b == Segment(point(0.0, 0.5, 0.5), point(0.5, 1.0, 1.0))

    r = Ray(point(3.0, 0.0, 0.5), vector(-1.0, 1.0, 0.0))
    @test intersection(r, b) |> type == NotIntersecting

    r = Ray(point(2.0, 0.0, 0.5), vector(-1.0, 1.0, 0.0))
    @test intersection(r, b) |> type == Touching
    @test r ∩ b == point(1.0, 1.0, 0.5)

    # the ray on a face of the box, got NaN in calculation
    r = Ray(point(1.5, 0.0, 0.0), vector(-1.0, 1.0, 0.0))
    @test intersection(r, b) |> type == Crossing
    @test r ∩ b == Segment(point(1.0, 0.5, 0.0), point(0.5, 1.0, 0.0))
  end

  @testset "Triangles" begin
    # utility to reverse segments, to more fully
    # test branches in the intersection algorithm
    reverse_segment(s) = Segment(vertices(s)[2], vertices(s)[1])

    # intersections with triangle lying in XY plane
    t = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 0))

    # intersects through t
    s = Segment(point(0.2, 0.2, 1.0), point(0.2, 0.2, -1.0))
    @test intersection(s, t) |> type == Intersecting
    @test s ∩ t == point(0.2, 0.2, 0.0)

    # intersects at a vertex of t
    s = Segment(point(0.0, 0.0, 1.0), point(0.0, 0.0, -1.0))
    @test intersection(s, t) |> type == Intersecting
    @test s ∩ t == point(0.0, 0.0, 0.0)

    # normal to, doesn't intersect with t
    s = Segment(point(0.9, 0.9, 1.0), point(0.9, 0.9, -1.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # coplanar, doesn't intersect with t
    s = Segment(point(-0.2, -0.2, 0.0), point(1.2, -0.2, 0.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # parallel, above, doesn't intersect with t
    s = Segment(point(-0.2, 0.2, 1.0), point(1.2, 0.2, 1.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # parallel, below, doesn't intersect with t
    s = Segment(point(-0.2, 0.2, -1.0), point(1.2, 0.2, -1.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # coplanar, within bounding box of t, no intersection
    s = Segment(point(0.7, 0.8, 0.0), point(0.8, 0.7, 0.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # segment above and to right of t, no intersection
    s = Segment(point(1.0, 1.0, 0.0), point(1.0, 1.0, 1.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # segment below t, no intersection
    s = Segment(point(0.5, -1.0, 0.0), point(0.5, -1.0, 1.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # segment left of t, no intersection
    s = Segment(point(-1.0, 0.5, 0.0), point(-1.0, 0.5, 1.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # segment above and to right of t, no intersection
    s = Segment(point(1.0, 1.0, 0.0), point(1.0, 1.0, -1.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)
    @test intersection(reverse_segment(s), t) |> type == NotIntersecting
    @test isnothing(reverse_segment(s) ∩ t)

    # segment below t, no intersection
    s = Segment(point(0.5, -1.0, 0.0), point(0.5, -1.0, -1.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)
    @test intersection(reverse_segment(s), t) |> type == NotIntersecting
    @test isnothing(reverse_segment(s) ∩ t)

    # segment left of t, no intersection
    s = Segment(point(-1.0, 0.5, 0.0), point(-1.0, 0.5, -1.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)
    @test intersection(reverse_segment(s), t) |> type == NotIntersecting
    @test isnothing(reverse_segment(s) ∩ t)

    # segment above and to right of t, no intersection
    s = Segment(point(1.0, 1.0, 1.0), point(1.0, 1.0, 0.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # segment below t, no intersection
    s = Segment(point(0.5, -1.0, 1.0), point(0.5, -1.0, 0.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # segment left of t, no intersection
    s = Segment(point(-1.0, 0.5, 1.0), point(-1.0, 0.5, 0.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # intersections with an inclined inclined triangle t
    t = Triangle(point(0, 0, 0), point(2, 0, 0), point(0, 2, 2))

    # doesn't reach t, no intersection
    s = Segment(point(0.5, 0.5, 1.9), point(0.5, 0.5, 1.8))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # parallel, offset from t, no intersection
    s = Segment(point(0.0, 0.5, 1.0), point(1.0, 0.5, 1.0))
    @test intersection(s, t) |> type == NotIntersecting
    @test isnothing(s ∩ t)

    # triangle as first argument
    t = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 0))
    s = Segment(point(0.2, 0.2, 1.0), point(0.2, 0.2, -1.0))
    @test intersection(t, s) |> type == Intersecting
    @test s ∩ t == t ∩ s == point(0.2, 0.2, 0.0)

    # type stability tests
    s = Segment(point(0.2, 0.2, 1.0), point(0.2, 0.2, -1.0))
    t = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 0))
    @inferred someornone(s, t)

    # https://github.com/JuliaGeometry/Meshes.jl/issues/728
    s = Segment(point(0.5, 0.5, 0.0), point(0.5, 0.5, 2.0))
    t = Triangle(point(1.0, 0.0, 0.0), point(0.0, 1.0, 0.0), point(0.0, 0.0, 1.0))
    @test intersection(s, t) |> type == Intersecting
    @test s ∩ t == t ∩ s == point(0.5, 0.5, 0.0)
    s = Segment(point(0.5, 0.5, 2.0), point(0.5, 0.5, 0.0))
    @test intersection(s, t) |> type == Intersecting
    @test s ∩ t == t ∩ s == point(0.5, 0.5, 0.0)

    # Intersection for a triangle and a ray
    t = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 0))

    # intersects through t
    r = Ray(point(0.2, 0.2, 1.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == Crossing
    @test r ∩ t == point(0.2, 0.2, 0.0)
    # origin of ray intersects with middle of triangle
    r = Ray(point(0.2, 0.2, 0.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == Touching
    @test r ∩ t == point(0.2, 0.2, 0.0)
    # Special case: the direction vector is not length enough to cross triangle
    r = Ray(point(0.2, 0.2, 1.0), vector(0.0, 0.0, -0.00001))
    @test intersection(r, t) |> type == Crossing
    if T == Float64
      @test r ∩ t ≈ point(0.2, 0.2, 0.0)
    end
    # Special case: reverse direction vector should not hit the triangle
    r = Ray(point(0.2, 0.2, 1.0), vector(0.0, 0.0, 1.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # intersects at a vertex of t
    r = Ray(point(0.0, 0.0, 1.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == CornerCrossing
    @test r ∩ t ≈ point(0.0, 0.0, 0.0)

    # normal to, doesn't intersect with t
    r = Ray(point(0.9, 0.9, 1.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # coplanar, doesn't intersect with t
    r = Ray(point(-0.2, -0.2, 0.0), vector(1.0, 0.0, 0.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # parallel, above, doesn't intersect with t
    r = Ray(point(-0.2, 0.2, 1.0), vector(1.0, 0.0, 0.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # parallel, below, doesn't intersect with t
    r = Ray(point(-0.2, 0.2, -1.0), vector(1.0, 0.0, 0.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # coplanar, within bounding box of t, no intersection
    r = Ray(point(0.7, 0.8, 0.0), vector(1.0, -1.0, 0.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # ray above and to right of t, no intersection
    r = Ray(point(1.0, 1.0, 0.0), vector(0.0, 0.0, 1.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # ray below t, no intersection
    r = Ray(point(0.5, -1.0, 0.0), vector(0.0, 0.0, 1.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # ray left of t, no intersection
    r = Ray(point(-1.0, 0.5, 0.0), vector(0.0, 0.0, 1.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # ray above and to right of t, no intersection
    r = Ray(point(1.0, 1.0, 0.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # ray below t, no intersection
    r = Ray(point(0.5, -1.0, 0.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # ray left of t, no intersection
    r = Ray(point(-1.0, 0.5, 0.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # ray above and to right of t, no intersection
    r = Ray(point(1.0, 1.0, 1.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # ray below t, no intersection
    r = Ray(point(0.5, -1.0, 1.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # ray left of t, no intersection
    r = Ray(point(-1.0, 0.5, 1.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # intersections with an inclined inclined triangle t
    t = Triangle(point(0, 0, 0), point(2, 0, 0), point(0, 2, 2))

    # doesn't reach t, but a ray can hit the triangle
    r = Ray(point(0.5, 0.5, 1.9), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == Crossing
    @test r ∩ t ≈ point(0.5, 0.5, 0.5)

    # parallel, offset from t, no intersection
    r = Ray(point(0.0, 0.5, 1.0), vector(1.0, 0.0, 0.0))
    @test intersection(r, t) |> type == NotIntersecting
    @test isnothing(r ∩ t)

    # origin of ray intersects with vertex of triangle
    r = Ray(point(0.0, 0.0, 0.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == CornerTouching
    @test r ∩ t ≈ point(0.0, 0.0, 0.0)

    # origin of ray intersects with edge of triangle
    r = Ray(point(0.5, 0.0, 0.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == EdgeTouching
    @test r ∩ t ≈ point(0.5, 0.0, 0.0)

    # ray intersects with edge of triangle
    r = Ray(point(0.5, 0.0, 1.0), vector(0.0, 0.0, -1.0))
    @test intersection(r, t) |> type == EdgeCrossing
    @test r ∩ t ≈ point(0.5, 0.0, 0.0)
  end

  @testset "Ngons" begin
    o = Octagon(
      point(0.0, 0.0, 1.0),
      point(0.5, -0.5, 0.0),
      point(1.0, 0.0, 0.0),
      point(1.5, 0.5, -0.5),
      point(1.0, 1.0, 0.0),
      point(0.5, 1.5, 0.0),
      point(0.0, 1.0, 0.0),
      point(-0.5, 0.5, 0.0)
    )

    r = Ray(point(-1.0, -1.0, -1.0), vector(1.0, 1.0, 1.0))
    @test intersection(r, o) |> type == Intersecting
    @test r ∩ o == PointSet([point(0.0, 0.0, 0.0)])

    r = Ray(point(-1.0, -1.0, -1.0), vector(-1.0, -1.0, -1.0))
    @test intersection(r, o) |> type == NotIntersecting
    @test isnothing(r ∩ o)
  end

  @testset "Polygons" begin
    # triangle
    poly = Triangle(point(6, 2), point(3, 5), point(0, 2))
    other = Quadrangle(point(5, 0), point(5, 4), point(0, 4), point(0, 0))
    @test intersection(poly, other) |> type == Intersecting
    @test all(vertices(poly ∩ other) .≈ [point(5, 3), point(4, 4), point(2, 4), point(0, 2), point(5, 2)])

    # octagon
    poly =
      Octagon(point(8, -2), point(8, 5), point(2, 5), point(4, 3), point(6, 3), point(4, 1), point(2, 1), point(2, -2))
    other = Quadrangle(point(5, 0), point(5, 4), point(0, 4), point(0, 0))
    @test intersection(poly, other) |> type == Intersecting
    @test all(
      vertices(poly ∩ other) .≈ [
        point(3, 4),
        point(4, 3),
        point(5, 3),
        point(5, 2),
        point(4, 1),
        point(2, 1),
        point(2, 0),
        point(5, 0),
        point(5, 4)
      ]
    )

    # inside
    poly = Quadrangle(point(1, 0), point(1, 1), point(0, 1), point(0, 0))
    other = Quadrangle(point(5, 0), point(5, 4), point(0, 4), point(0, 0))
    @test intersection(poly, other) |> type == Intersecting
    @test all(vertices(poly ∩ other) .≈ vertices(poly))

    # outside
    poly = Quadrangle(point(7, 6), point(7, 7), point(6, 7), point(6, 6))
    other = Quadrangle(point(5, 0), point(5, 4), point(0, 4), point(0, 0))
    @test intersection(poly, other) |> type == NotIntersecting
    @test isnothing(poly ∩ other)

    # convex and non-convex polygons
    quad = Quadrangle(point(0, 0), point(0.1, 0.0), point(0.1, 0.1), point(0.0, 0.1))
    poly = PolyArea(point(0, 0), point(2, 0), point(1, 1), point(1, 0.5))
    @test intersection(quad, poly) |> type == Intersecting
    @test all(vertices(quad ∩ poly) .≈ [point(0, 0), point(0.1, 0), point(0.1, 0.05)])
  end

  @testset "Domains" begin
    grid = cartgrid(4, 4)
    pset = PointSet(centroid.(grid))
    ball = Ball(point(0, 0), T(1))
    @test pset ∩ pset == pset
    @test pset ∩ grid == grid ∩ pset == pset
    @test pset ∩ ball == ball ∩ pset == PointSet(point(0.5, 0.5))
  end
end
