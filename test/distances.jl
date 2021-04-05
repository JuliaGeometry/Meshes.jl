@testset "Distances" begin
    @testset "mindistance" begin
        p = P2(0, 1)
        l = Line(P2(0, 0), P2(1, 0))
        for metric in (Euclidean(), SqEuclidean())
            @test (@test_deprecated evaluate(metric, p, l)) == T(1)
            @test (@test_deprecated evaluate(metric, l, p)) == T(1)
            @test mindistance(metric, l, p) == mindistance(metric, p, l) == T(1)
        end

        p1, p2 = P2(1, 0), P2(0, 1)
        @test (@test_deprecated evaluate(Chebyshev(), p1, p2)) == T(1)
        @test mindistance(Chebyshev(), p1, p2) == T(1)

        p = P2(68, 259)
        l = Line(P2(68, 260), P2(69, 261))
        @test (@test_deprecated evaluate(Euclidean(), p, l)) ≤ T(0.8)
        @test mindistance(Euclidean(), l, p) == mindistance(Euclidean(), p, l) ≤ T(0.8)

        s = Segment(P2(0, 0), P2(1, 0))
        for metric in (Euclidean(), SqEuclidean())
            @test mindistance(metric, s, P2(0, 1)) == mindistance(metric, P2(0, 1), s) == T(1)
            @test mindistance(metric, s, P2(2, 0)) == mindistance(metric, P2(2, 0), s) == T(1)
            @test mindistance(metric, s, first(vertices(s))) ==
                  mindistance(metric, first(vertices(s)), s) ==
                  T(0)
        end

        c = Chain(P2(0, 0), P2(1, 0), P2(1, 1))
        for metric in (Euclidean(), SqEuclidean())
            @test mindistance(metric, c, P2(0, 1)) == mindistance(metric, P2(0, 1), c) == T(1)
            @test mindistance(metric, c, P2(2, 0)) == mindistance(metric, P2(2, 0), c) == T(1)
            @test mindistance(metric, c, first(vertices(c))) ==
                  mindistance(metric, first(vertices(c)), c) ==
                  T(0)
        end

    end

    @testset "closest_point" begin
        p = P2(0, 1)
        l = Line(P2(0, 0), P2(1, 0))

        for metric in (Euclidean(), SqEuclidean())
            @test closest_point(metric, l, p) == P2(0, 0)
            @test closest_point(metric, l, l.a) == l.a
            @test closest_point(metric, l, l.b) == l.b
            @test_throws MethodError closest_point(metric, p, l)
        end

        s = Segment(P2(0, 0), P2(1, 0))
        for metric in (Euclidean(), SqEuclidean())
            p_center = P2(((coordinates(s.vertices[1]) + coordinates(s.vertices[2])) / 2)...)
            @test closest_point(metric, s, P2(0, 1)) == s.vertices[1]
            @test closest_point(metric, s, P2(2, 0)) == s.vertices[2]
            @test closest_point(metric, s, s.vertices[1]) == s.vertices[1]
            @test closest_point(metric, s, p_center) == p_center
        end

        c = Chain(P2(0, 0), P2(1, 0), P2(1, 1))
        for metric in (Euclidean(), SqEuclidean())
            @test closest_point(metric, c, P2(0, 1)) in (c.vertices[1], c.vertices[3])
            @test closest_point(metric, c, P2(2, 0)) == c.vertices[2]
            @test closest_point(metric, c, c.vertices[1]) == c.vertices[1]
        end
    end
end
