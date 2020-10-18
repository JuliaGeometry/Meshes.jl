using Test, Meshes

@testset "algorithms.jl" begin
    cube = Box(Point(-0.5,-0.5,-0.5), Vec(1.0,1.0,1.0))
    cube_faces = decompose(TriangleFace, faces(cube))
    cube_vertices = decompose(Point3f, cube)
    @test area(cube_vertices, cube_faces) == 6
    mesh = Mesh(cube_vertices, cube_faces)
    @test Meshes.volume(mesh) ≈ 1
end

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
        target = Point3[(4.535533905932738, -1.5355339059327373, 3.0),
                        (7.535533905932738, 1.4644660940672627, 6.0),
                        (3.0412414523193148, 4.041241452319315, -1.0824829046386295),
                        (6.041241452319315, 7.041241452319315, 1.9175170953613705),
                        (-2.535533905932737, 5.535533905932738, 2.9999999999999996),
                        (0.46446609406726314, 8.535533905932738, 6.0),
                        (-1.0412414523193152, -0.04124145231931431, 7.0824829046386295),
                        (1.9587585476806848, 2.9587585476806857, 10.08248290463863),
                        (1, 2, 3),
                        (4, 5, 6)]
        points = decompose(Point3, Tesselation(c, 8))
        @test coordinates.(points) ≈ coordinates.(target)
        faces = TriangleFace.([(3, 2, 1), (4, 2, 3), (5, 4, 3), (6, 4, 5), (7, 6, 5),
                               (8, 6, 7), (1, 8, 7), (2, 8, 1), (3, 1, 9), (2, 4, 10),
                               (5, 3, 9), (4, 6, 10), (7, 5, 9), (6, 8, 10), (1, 7, 9),
                               (8, 2, 10)])
        @test faces == decompose(TriangleFace, Tesselation(c, 8))

        m = triangle_mesh(Tesselation(c, 8))

        @test Meshes.faces(m) == faces
        points = coordinates(m)
        @test coordinates.(points) ≈ coordinates.(target)
    end
end

@testset "Boxes" begin
    a = Box(Point(0, 0), Vec(1, 1))
    pt_expa = Point[(0, 0), (1, 0), (0, 1), (1, 1)]
    @test decompose(Point{2,Int}, a) == pt_expa
    mesh = triangle_mesh(a)
    @test decompose(Point2f, mesh) == convert.(Point2f, pt_expa)

    b = Box(Point(1,1,1), Vec(1,1,1))
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
    @test ndims(sphere) == 3
    @test coordtype(sphere) == Float32
    @test center(sphere) == Point3f(0,0,0)
    @test radius(sphere) == 1.0f0

    points = decompose(Point3f, Tesselation(sphere, 3))
    target = Point3f[(0.0, 0.0, 1.0),
                     (1.0, 0.0, 6.12323e-17),
                     (1.22465e-16, 0.0, -1.0),
                     (-0.0, 0.0, 1.0),
                     (-1.0, 1.22465e-16, 6.12323e-17),
                     (-1.22465e-16, 1.49976e-32, -1.0),
                     (0.0, -0.0, 1.0),
                     (1.0, -2.44929e-16, 6.12323e-17),
                     (1.22465e-16, -2.99952e-32, -1.0)]
    @test coordinates.(points) ≈ coordinates.(target)

    f = decompose(TriangleFace, Tesselation(sphere, 3))
    face_target = TriangleFace.([(1, 2, 5), (1, 5, 4), (2, 3, 6), (2, 6, 5),
                                 (4, 5, 8), (4, 8, 7), (5, 6, 9), (5, 9, 8)])
    @test f == face_target
    circle = Sphere(Point2f(0, 0), 1.0f0)
    points = decompose(Point2f, Tesselation(circle, 20))
    @test length(points) == 20
    tess = Tesselation(circle, 32)
    mesh = triangle_mesh(tess)
    mpoints = decompose(Point2f, mesh)
    tpoints = decompose(Point2f, tess)
    @test coordinates.(mpoints) ≈ coordinates.(tpoints)
end

@testset "Boxes" begin
    rect = Box(Point2(0, 0), Vec2(1, 2))
    @test rect isa Box{2,Float64}

    split1, split2 = Meshes.split(rect, 2, 1)
    @test widths(split1) == widths(split2)
    @test origin(split1) == Point(0.0, 0.0)
    @test origin(split2) == Point(0.0, 1.0)
    @test split1 ∈ rect
    @test rect ∉ split1

    prim = Box(Point2(0, 0), Vec2(1, 1))
    @test length(prim) == 2

    p = Point(1.0, 1.0)
    r = Box(Point(0.0, 0.0), Vec(1.0, 1.0))
    @test p ∈ r

    h1 = Box(Point2(0.0, 0.0), Vec2(1.0, 1.0))
    h2 = Box(Point2(1.0, 1.0), Vec2(2.0, 2.0))
    @test union(h1, h2) isa Box{2,Float64}
    @test Meshes.intersect(h1, h2) isa Box{2,Float64}

    rect1 = Box(Point2(0.0, 0.0), Vec2(1.0, 1.0))
    rect2 = Box(Point2(3.0, 1.0), Vec2(4.0, 2.0))
    @test !before(rect1, rect2)
    rect1 = Box(Point2(0.0, 0.0), Vec2(1.0, 1.0))
    rect2 = Box(Point2(3.0, 2.0), Vec2(4.0, 2.0))
    @test before(rect1, rect2)

    rect1 = Box(Point2(1.0, 1.0), Vec2(2.0, 2.0))
    rect2 = Box(Point2(0.0, 0.0), Vec2(2.0, 1.0))
    @test !overlaps(rect1, rect2)
    rect1 = Box(Point2(1.0, 1.0), Vec2(2.0, 2.0))
    rect2 = Box(Point2(1.5, 1.5), Vec2(2.0, 2.0))
    @test overlaps(rect1, rect2)

    rect1 = Box(Point2(1.0, 1.0), Vec2(2.0, 2.0))
    rect2 = Box(Point2(0.0, 0.0), Vec2(2.0, 1.0))
    @test !Meshes.starts(rect1, rect2)
    rect2 = Box(Point2(1.0, 1.0), Vec2(1.5, 1.5))
    @test !Meshes.starts(rect1, rect2)
    rect2 = Box(Point2(1.0, 1.0), Vec2(3.0, 3.0))
    @test Meshes.starts(rect1, rect2)

    rect1 = Box(Point2(1.0, 1.0), Vec2(2.0, 2.0))
    rect2 = Box(Point2(0.0, 0.0), Vec2(4.0, 4.0))
    @test during(rect1, rect2)
    rect1 = Box(Point2(0.0, 0.0), Vec2(2.0, 3.0))
    rect2 = Box(Point2(1.0, 1.0), Vec2(4.0, 2.0))
    @test !during(rect1, rect2)

    rect1 = Box(Point2(1.0, 1.0), Vec2(2.0, 2.0))
    rect2 = Box(Point2(0.0, 0.0), Vec2(4.0, 4.0))
    @test !finishes(rect1, rect2)
    rect1 = Box(Point2(1.0, 0.0), Vec2(1.0, 1.0))
    rect2 = Box(Point2(0.0, 0.0), Vec2(2.0, 1.0))
    @test !finishes(rect1, rect2)
    rect1 = Box(Point2(1.0, 1.0), Vec2(1.0, 2.0))
    rect2 = Box(Point2(0.0, 0.0), Vec2(2.0, 3.0))
    @test finishes(rect1, rect2)
end
