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
    end
end
