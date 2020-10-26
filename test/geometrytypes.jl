using Test, Meshes

@testset "Cylinder" begin
    @testset "constructors" begin
        c = Cylinder(Point3(1,2,3), Point3(4,5,6), 5.0)
        @test ndims(c) == 3
        @test coordtype(c) == Float64
        @test radius(c) == 5.0
        @test height(c) == √27
        @test volume(c) == π*5.0^2*√27
    end

    @testset "decompose" begin
        c = Cylinder(Point3(1,2,3), Point3(4,5,6), 5.0)
        points = decompose(Point3, c)
    end
end

@testset "Boxes" begin
    a = Box(Point(0, 0), Point(1, 1))
    pt_expa = Point[(0, 0), (1, 0), (0, 1), (1, 1)]
    @test decompose(Point{2,Int}, a) == pt_expa
    mesh = Meshes.mesh(a)
    @test decompose(Point2f, mesh) == convert.(Point2f, pt_expa)

    b = Box(Point(1,1,1), Point(2,2,2))
    pt_expb = Point[(1, 1, 1), (1, 1, 2), (1, 2, 2), (1, 2, 1), (1, 1, 1),
                    (2, 1, 1), (2, 1, 2), (1, 1, 2), (1, 1, 1), (1, 2, 1),
                    (2, 2, 1), (2, 1, 1), (2, 2, 2), (1, 2, 2), (1, 1, 2),
                    (2, 1, 2), (2, 2, 2), (2, 1, 2), (2, 1, 1), (2, 2, 1),
                    (2, 2, 2), (2, 2, 1), (1, 2, 1), (1, 2, 2)]
    @test decompose(Point{3,Int}, b) == pt_expb
end

NFace = NgonFace

@testset "Faces" begin
    @test convert_simplex(TriangleFace, QuadFace((1, 2, 3, 4))) ==
          (TriangleFace((1, 2, 3)), TriangleFace((1, 3, 4)))
    @test convert_simplex(LineFace, QuadFace((1, 2, 3, 4))) ==
          (LineFace((1, 2)), LineFace((2, 3)), LineFace((3, 4)), LineFace((4, 1)))
end

@testset "Sphere" begin
    sphere = Sphere(Point3f(0,0,0), 1.0f0)
    points = decompose(Point3f, sphere)
    @test ndims(sphere) == 3
    @test coordtype(sphere) == Float32
    @test center(sphere) == Point3f(0,0,0)
    @test radius(sphere) == 1.0f0

    circle = Sphere(Point2f(0, 0), 1.0f0)
    points = decompose(Point2f, circle)
end

@testset "Boxes" begin
    p = Point(1.0, 1.0)
    r = Box(Point(0.0, 0.0), Point(1.0, 1.0))
    @test p ∈ r

    h1 = Box(Point(0.0, 0.0), Point(1.0, 1.0))
    h2 = Box(Point(1.0, 1.0), Point(3.0, 3.0))
    @test union(h1, h2) isa Box{2,Float64}
    @test Meshes.intersect(h1, h2) isa Box{2,Float64}

    rect1 = Box(Point(0.0, 0.0), Point(1.0, 1.0))
    rect2 = Box(Point(3.0, 1.0), Point(7.0, 3.0))
    @test !before(rect1, rect2)
    rect1 = Box(Point(0.0, 0.0), Point(1.0, 1.0))
    rect2 = Box(Point(3.0, 2.0), Point(7.0, 4.0))
    @test before(rect1, rect2)

    rect1 = Box(Point(1.0, 1.0), Point(3.0, 3.0))
    rect2 = Box(Point(0.0, 0.0), Point(2.0, 1.0))
    @test !overlaps(rect1, rect2)
    rect1 = Box(Point(1.0, 1.0), Point(3.0, 3.0))
    rect2 = Box(Point(1.5, 1.5), Point(3.5, 3.5))
    @test overlaps(rect1, rect2)

    rect1 = Box(Point(1.0, 1.0), Point(3.0, 3.0))
    rect2 = Box(Point(0.0, 0.0), Point(2.0, 1.0))
    @test !Meshes.starts(rect1, rect2)
    rect2 = Box(Point(1.0, 1.0), Point(2.5, 2.5))
    @test !Meshes.starts(rect1, rect2)
    rect2 = Box(Point(1.0, 1.0), Point(4.0, 4.0))
    @test Meshes.starts(rect1, rect2)

    rect1 = Box(Point(1.0, 1.0), Point(3.0, 3.0))
    rect2 = Box(Point(0.0, 0.0), Point(4.0, 4.0))
    @test during(rect1, rect2)
    rect1 = Box(Point(0.0, 0.0), Point(2.0, 3.0))
    rect2 = Box(Point(1.0, 1.0), Point(5.0, 3.0))
    @test !during(rect1, rect2)

    rect1 = Box(Point(1.0, 1.0), Point(3.0, 3.0))
    rect2 = Box(Point(0.0, 0.0), Point(4.0, 4.0))
    @test !finishes(rect1, rect2)
    rect1 = Box(Point(1.0, 0.0), Point(2.0, 1.0))
    rect2 = Box(Point(0.0, 0.0), Point(2.0, 1.0))
    @test !finishes(rect1, rect2)
    rect1 = Box(Point(1.0, 1.0), Point(2.0, 3.0))
    rect2 = Box(Point(0.0, 0.0), Point(2.0, 3.0))
    @test finishes(rect1, rect2)
end
