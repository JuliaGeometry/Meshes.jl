using Meshes
using StaticArrays
using LinearAlgebra
using Test, Random

@testset "Meshes" begin

@testset "view" begin
    @testset "TupleView" begin
        x = [1, 2, 3, 4, 5, 6]
        y = TupleView{2, 1}(x)
        @test y == [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]

        y = TupleView{2}(x)
        @test y == [(1, 2), (3, 4), (5, 6)]

        y = TupleView{2, 3}(x)
        @test y == [(1, 2), (4, 5)]

        y = TupleView{3, 1}(x)
        @test y == [(1, 2, 3), (2, 3, 4), (3, 4, 5), (4, 5, 6)]

        y = TupleView{2, 1}(x, connect = true)
        @test y == [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]

    end

    @testset "connected views" begin
        points = Point{2,Int}[(1, 2), (3, 4), (5, 6)]
        lines = connect(points, Segment{2,Int}, 1)
        @test lines == [Segment(Point(1, 2), Point(3, 4)), Segment(Point(3, 4), Point(5, 6))]
        triangles = connect(points, Triangle{2,Int})
        @test triangles == [Triangle(Point(1, 2), Point(3, 4), Point(5, 6))]
    end

    @testset "face views" begin
        points = Point.([(1, 2), (3, 4), (5, 6)])
        faces = connect([1, 2, 3], TriangleFace)
        triangles = connect(points, faces)
        @test triangles == [Triangle(Point(1, 2), Point(3, 4), Point(5, 6))]

        x = Point3(1,1,1)
        triangles = connect([x], [TriangleFace((1, 1, 1))])
        @test triangles == [Triangle(x, x, x)]

        points = Point3[(1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12)]
        faces = connect([1, 2, 3, 4], TetrahedronFace)
        tetrahedron = connect(points, faces)[1]
        @test ndims(tetrahedron) == 3
        @test length(tetrahedron) == 4
    end

    @testset "reinterpret" begin
        numbers = collect(reshape(1:6, 2, 3))
        points = reinterpret(Point{2, Int}, numbers)
        @test points[1] === Point(1, 2)
        @test points[2] === Point(3, 4)
        numbers[4] = 0
        @test points[2] === Point(3, 0)
    end
end

@testset "constructors" begin
    @testset "Mesh" begin
        points = Point2[(1, 2), (3, 4), (5, 6)]
        mesh = Mesh(points, [1,2,3])
        @test elements(mesh)[1] == Triangle(points...)

        x = Point3(1,1,1)
        mesh = Mesh([x], [TriangleFace((1, 1, 1))])
        @test elements(mesh)[1] == Triangle(x, x, x)

        points = rand(Point3f, 8)
        tfaces = [TriangleFace((1, 2, 3)), TriangleFace((5, 6, 7))]
        mesh = Mesh(points, tfaces)
        @test mesh isa Mesh

        points = Point3[(1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12)]
        sfaces = connect([1, 2, 3, 4], TetrahedronFace)
        mesh = Mesh(points, sfaces)
        @test elements(mesh)[1] == Tetrahedron(points...)

        b = Box(Point2f(0,0), Point2f(2,2))
        m = Meshes.mesh(b, facetype=QuadFace)
    end
end

@testset "decompose/triangulation" begin
    primitive = Sphere(Point3f(0,0,0), 1.0f0)
    @test ndims(primitive) === 3
    mesh = Meshes.mesh(primitive)
    @test decompose(Point3f, mesh) isa Vector{Point3f}
    @test decompose(Point3f, primitive) isa Vector{Point3f}

    primitive = Box(Point2(0, 0), Point2(1, 1))
    mesh = Meshes.mesh(primitive)
    @test decompose(Point2f, mesh) isa Vector{Point2f}
    @test decompose(Point{2,Int}, primitive) isa Vector{Point{2,Int}}

    primitive = Box(Point3(0,0,0), Point3(1,1,1))
    mesh = Meshes.mesh(primitive)
    @test decompose(Point3, mesh) isa Vector{Point3}
end

@testset "mesh conversion" begin
    s = Sphere(Point3(0,0,0), 1.0)
    m = Meshes.mesh(s)
    @test m isa Mesh{3,Float64}
    points1 = coordinates(m)
    points2 = decompose(Point3, m)
    @test coordinates.(points1) ≈ coordinates.(points2)

    s = Sphere(Point3f(0,0,0), 1.0f0)
    m = Meshes.mesh(s)
    @test m isa Mesh{3,Float32}
end

@testset "lines intersects" begin
    a = Segment(Point(0.0, 0.0), Point(4.0, 1.0))
    b = Segment(Point(0.0, 0.25), Point(3.0, 0.25))
    c = Segment(Point(0.0, 0.25), Point(0.5, 0.25))
    d = Segment(Point(0.0, 0.0), Point(0.0, 4.0))
    e = Segment(Point(1.0, 0.0), Point(0.0, 4.0))
    f = Segment(Point(5.0, 0.0), Point(6.0, 0.0))

    @test a ∩ b === (true, Point(1.0, 0.25))
    @test a ∩ c === (false, Point(0.0, 0.0))
    @test d ∩ d === (false, Point(0.0, 0.0))
    found, point = d ∩ e
    @test found && coordinates(point) ≈ [0.0, 4.0]
    @test a ∩ f === (false, Point(0.0, 0.0))
end

@testset "Tests from GeometryTypes" begin
    include("geometrytypes.jl")
end

end
