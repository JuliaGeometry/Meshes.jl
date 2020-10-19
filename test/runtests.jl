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
        numbers = [1, 2, 3, 4, 5, 6]
        x = connect(numbers, Point2)
        @test x == Point{2,Int}[(1, 2), (3, 4), (5, 6)]

        line = connect(x, Line, 1)
        @test line == [Line(Point(1, 2), Point(3, 4)), Line(Point(3, 4), Point(5, 6))]

        triangles = connect(x, Triangle)
        @test triangles == [Triangle(Point(1, 2), Point(3, 4), Point(5, 6))]
    end

    @testset "face views" begin
        numbers = [1, 2, 3, 4, 5, 6]
        points = connect(numbers, Point{2})
        faces = connect([1, 2, 3], TriangleFace)
        triangles = connect(points, faces)
        @test triangles == [Triangle(Point(1, 2), Point(3, 4), Point(5, 6))]

        x = Point3(1,1,1)
        triangles = connect([x], [TriangleFace((1, 1, 1))])
        @test triangles == [Triangle(x, x, x)]

        points = connect([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], Point3)
        faces = connect([1, 2, 3, 4], TetrahedronFace{Int})
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
    @testset "LineFace" begin

        points = connect([1, 2, 3, 4, 5, 6], Point{2})
        linestring = LineString(points)
        @test linestring == [Line(points[1], points[2]), Line(points[2], points[3])]

        points = rand(Point2, 4)
        linestring = LineString(points, 2)
        @test linestring == [Line(points[1], points[2]), Line(points[3], points[4])]

        linestring = LineString([points[1] => points[2], points[2] => points[3]])
        @test linestring == [Line(points[1], points[2]), Line(points[2], points[3])]

        faces = [1, 2, 3]
        linestring = LineString(points, faces)
        @test linestring == LineString([points[1] => points[2], points[2] => points[3]])
        a, b, c, d = Point(1, 2), Point(3, 4), Point(5, 6), Point(7, 8)
        points = [a, b, c, d]; faces = [1, 2, 3, 4]
        linestring = LineString(points, faces, 2)
        @test linestring == LineString([a => b, c => d])

        faces = [LineFace((1, 2)), LineFace((3, 4))]
        linestring = LineString(points, faces)
        @test linestring == LineString([a => b, c => d])
    end

    @testset "Polygon" begin

        points = connect([1, 2, 3, 4, 5, 6], Point2)
        polygon = Polygon(points)
        @test polygon == Polygon(LineString(points))

        points = rand(Point2, 4)
        linestring = LineString(points, 2)
        @test Polygon(points, 2) == Polygon(linestring)

        faces = [1, 2, 3]
        polygon = Polygon(points, faces)
        @test polygon == Polygon(LineString(points, faces))

        a, b, c, d = Point(1, 2), Point(3, 4), Point(5, 6), Point(7, 8)
        points = [a, b, c, d]; faces = [1, 2, 3, 4]
        polygon = Polygon(points, faces, 2)
        @test polygon == Polygon(LineString(points, faces, 2))

        faces = [LineFace((1, 2)), LineFace((3, 4))]
        polygon = Polygon(points, faces)
        @test polygon == Polygon(LineString(points, faces))
        @test ndims(polygon) === 2
    end

    @testset "Mesh" begin

        numbers = [1, 2, 3, 4, 5, 6]
        points = connect(numbers, Point2)
        mesh = Mesh(points, [1,2,3])
        @test mesh == [Triangle(points...)]

        x = Point3(1,1,1)
        mesh = Mesh([x], [TriangleFace((1, 1, 1))])
        @test mesh == [Triangle(x, x, x)]

        points = rand(Point3f, 8)
        tfaces = [TriangleFace((1, 2, 3)), TriangleFace((5, 6, 7))]
        mesh = Mesh(points, tfaces)
        @test mesh isa Mesh

        points = connect([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], Point3)
        sfaces = connect([1, 2, 3, 4], TetrahedronFace{Int})
        mesh = Mesh(points, sfaces)
        @test mesh == [Meshes.NNgon{4}(points...)]

        t = Tesselation(Box(Point2f(0,0), Point2f(2,2)), (30, 30))
        m = Meshes.mesh(t, facetype=QuadFace)
        @test Meshes.faces(m) isa Vector{QuadFace}
        @test Meshes.coordinates(m) isa Vector{Point2f}
    end

    @testset "Multi geometries" begin
        # coordinates from https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry#Geometric_objects
        points = Point{2, Int}[(10, 40), (40, 30), (20, 20), (30, 10)]
        multipoint = MultiPoint(points)
        @test size(multipoint) === size(points)
        @test multipoint[3] === points[3]

        linestring1 = LineString(Point{2, Int}[(10, 10), (20, 20), (10, 40)])
        linestring2 = LineString(Point{2, Int}[(40, 40), (30, 30), (40, 20), (30, 10)])
        multilinestring = MultiLineString([linestring1, linestring2])
        @test size(multilinestring) === (2,)
        @test multilinestring[1] === linestring1
        @test multilinestring[2] === linestring2

        polygon11 = Polygon(Point{2, Int}[(30, 20), (45, 40), (10, 40), (30, 20)])
        polygon12 = Polygon(Point{2, Int}[(15, 5), (40, 10), (10, 20), (5, 10), (15, 5)])
        multipolygon1 = MultiPolygon([polygon11, polygon12])
        @test size(multipolygon1) === (2,)
        @test multipolygon1[1] === polygon11
        @test multipolygon1[2] === polygon12

        polygon21 = Polygon(Point{2, Int}[(40, 40), (20, 45), (45, 30), (40, 40)])
        polygon22 = Polygon(LineString(Point{2, Int}[(20, 35), (10, 30), (10, 10), (30, 5), (45, 20), (20, 35)]),
            [LineString(Point{2, Int}[(30, 20), (20, 15), (20, 25), (30, 20)])])
        multipolygon2 = MultiPolygon([polygon21, polygon22])
        @test size(multipolygon2) === (2,)
        @test multipolygon2[1] === polygon21
        @test multipolygon2[2] === polygon22
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

    points = decompose(Point2f, Sphere(Point2f(0, 0), 1.0f0))
    mesh = Meshes.mesh(points)
    @test coordinates(mesh) == points

    linestring = LineString(Point{2, Int}[(10, 10), (20, 20), (10, 40)])
    pts = Point{2, Int}[(10, 10), (20, 20), (10, 40)]
    linestring = LineString(pts)
    pts_decomp = decompose(Point{2, Int}, linestring)
    @test pts == pts_decomp

    pts_ext = Point{2, Int}[(5, 1), (3, 3), (4, 8), (1, 2), (5, 1)]
    ls_ext = LineString(pts_ext)
    pts_int1 = Point{2, Int}[(2, 2), (3, 8),(5, 6), (3, 4), (2, 2)]
    ls_int1 = LineString(pts_int1)
    pts_int2 = Point{2, Int}[(3, 2), (4, 5),(6, 1), (1, 4), (3, 2)]
    ls_int2 =  LineString(pts_int2)
    poly_ext = Polygon(ls_ext)
    poly_ext_int = Polygon(ls_ext, [ls_int1, ls_int2])
    @test decompose(Point{2, Int}, poly_ext) == pts_ext
    @test decompose(Point{2, Int}, poly_ext_int) == [pts_ext..., pts_int1..., pts_int2...]
end

@testset "mesh conversion" begin
    s = Sphere(Point3(0,0,0), 1.0)
    m = Meshes.mesh(s)
    @test m isa Mesh{3,Float64}
    @test coordinates(m) isa Vector{Point3}
    @test Meshes.faces(m) isa Vector{TriangleFace}
    points1 = coordinates(m)
    points2 = decompose(Point3, m)
    @test coordinates.(points1) ≈ coordinates.(points2)

    s = Sphere(Point3f(0,0,0), 1.0f0)
    m = Meshes.mesh(s)
    tmesh = Meshes.mesh(m)
    points1 = coordinates(tmesh)
    points2 = decompose(Point3f, tmesh)
    @test tmesh isa Mesh
    @test coordinates.(points1) ≈ coordinates.(points2)
    @test m isa Mesh{3,Float32}
    @test coordinates(m) isa Vector{Point3f}
    @test Meshes.faces(m) isa Vector{TriangleFace}
end

@testset "lines intersects" begin
    a = Line(Point(0.0, 0.0), Point(4.0, 1.0))
    b = Line(Point(0.0, 0.25), Point(3.0, 0.25))
    c = Line(Point(0.0, 0.25), Point(0.5, 0.25))
    d = Line(Point(0.0, 0.0), Point(0.0, 4.0))
    e = Line(Point(1.0, 0.0), Point(0.0, 4.0))
    f = Line(Point(5.0, 0.0), Point(6.0, 0.0))

    @test intersects(a, b) === (true, Point(1.0, 0.25))
    @test intersects(a, c) === (false, Point(0.0, 0.0))
    @test intersects(d, d) === (false, Point(0.0, 0.0))
    found, point = intersects(d, e)
    @test found && coordinates(point) ≈ [0.0, 4.0]
    @test intersects(a, f) === (false, Point(0.0, 0.0))
end

@testset "Tests from GeometryTypes" begin
    include("geometrytypes.jl")
end

end
