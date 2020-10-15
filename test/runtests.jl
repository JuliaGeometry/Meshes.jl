using Meshes
using Meshes: attributes
using StructArrays
using Tables
using StaticArrays
using LinearAlgebra
using Test, Random

@testset "Meshes" begin

@testset "embedding metadata" begin
    @testset "Meshes" begin

        @testset "per vertex attributes" begin
            points = rand(Point3, 8)
            tfaces = TetrahedronFace{Int}[(1, 2, 3, 4), (5, 6, 7, 8)]
            normals = rand(Vec3, 8)
            stress = LinRange(0, 1, 8)
            mesh = Mesh(meta(points, normals = normals, stress = stress), tfaces)

            @test hasproperty(coordinates(mesh), :stress)
            @test hasproperty(coordinates(mesh), :normals)
            @test coordinates(mesh).stress === stress
            @test coordinates(mesh).normals === normals
            @test coordinates(mesh).normals === normals
            @test Meshes.faces(mesh) === tfaces
            @test propertynames(coordinates(mesh)) == (:position, :normals, :stress)

        end

        @testset "per face attributes" begin

            # Construct a cube out of Quads
            points = Point{3, Float64}[
                (0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                (2.0, 2.0, 0.0), (0.0, 2.0, 0.0),
                (0.0, 0.0, 12.0), (2.0, 0.0, 12.0),
                (2.0, 2.0, 12.0), (0.0, 2.0, 12.0)
            ]

            facets = QuadFace[
                1:4,
                5:8,
                [1,5,6,2],
                [2,6,7,3],
                [3, 7, 8, 4],
                [4, 8, 5, 1]
            ]

            markers = [-1, -2, 0, 0, 0, 0]
            # attach some additional information to our faces!
            mesh = Mesh(points, meta(facets, markers = markers))
            @test hasproperty(Meshes.faces(mesh), :markers)
            # test with === to assert we're not doing any copies
            @test Meshes.faces(mesh).markers === markers
            @test coordinates(mesh) === points
            @test metafree(Meshes.faces(mesh)) === facets

        end

    end

    @testset "polygon with metadata" begin
        polys = [Polygon(rand(Point2f, 20)) for i in 1:10]
        pnames = [randstring(4) for i in 1:10]
        numbers = LinRange(0.0, 1.0, 10)
        bin = rand(Bool, 10)
        # create a polygon
        poly = PolygonMeta(polys[1], name = pnames[1], value = numbers[1], category = bin[1])
        # create a MultiPolygon with the right type & meta information!
        multipoly = MultiPolygonMeta(polys, name = pnames, value = numbers, category = bin)
        @test multipoly isa AbstractVector
        @test poly isa Meshes.AbstractPolygon

        @test Meshes.getcolumn(poly, :name) == pnames[1]
        @test Meshes.MetaFree(PolygonMeta) == Polygon

        @test Meshes.getcolumn(multipoly, :name) == pnames
        @test Meshes.MetaFree(MultiPolygonMeta) == MultiPolygon

        meta_p = meta(polys[1], boundingbox=Rectangle(Point(0,0), Vec(2,2)))
        @test meta_p.boundingbox === Rectangle(Point(0,0), Vec(2,2))
        @test metafree(meta_p) === polys[1]
        attributes(meta_p) == Dict{Symbol, Any}(:boundingbox => meta_p.boundingbox,
                                                :polygon => polys[1])
    end
    @testset "point with metadata" begin
        p = Point(1.1, 2.2)
        pm = Meshes.PointMeta(1.1, 2.2; a=1, b=2)
        @test meta(pm) === (a=1, b=2)
        @test metafree(pm) === p
        @test propertynames(pm) == (:position, :a, :b)
    end

    @testset "MultiPoint with metadata" begin
        p = collect(Point{2, Float64}(x, x+1) for x in 1:5)
        @test p isa AbstractVector
        mpm = MultiPointMeta(p, a=1, b=2)
        @test coordinates(mpm) == mpm
        @test meta(mpm) === (a=1, b=2)
        @test metafree(mpm) == p
        @test propertynames(mpm) == (:points, :a, :b)
    end

    @testset "LineString with metadata" begin
        linestring = LineStringMeta(Point{2, Int}[(10, 10), (20, 20), (10, 40)], a = 1, b = 2)
        @test linestring isa AbstractVector
        @test meta(linestring) === (a = 1, b = 2)
        @test metafree(linestring) == linestring
        @test propertynames(linestring) == (:lines, :a, :b)
    end

    @testset "MultiLineString with metadata" begin
        linestring1 = LineString(Point{2,Int}[(10, 10), (20, 20), (10, 40)])
        linestring2 = LineString(Point{2,Int}[(40, 40), (30, 30), (40, 20), (30, 10)])
        multilinestring = MultiLineString([linestring1, linestring2])
        multilinestringmeta = MultiLineStringMeta([linestring1, linestring2]; boundingbox = Rectangle(Point(1.0,1.0), Vec(2.0,2.0)))
        @test multilinestringmeta isa AbstractVector
        @test meta(multilinestringmeta) === (boundingbox = Rectangle(Point(1.0,1.0), Vec(2.0,2.0)),)
        @test metafree(multilinestringmeta) == multilinestring
        @test propertynames(multilinestringmeta) == (:linestrings, :boundingbox)
    end

    @testset "Mesh with metadata" begin
       m = triangle_mesh(HyperSphere(Point3f(0,0,0), 1.0f0))
       m_meta = MeshMeta(m; boundingbox = Rectangle(Point(1.0,1.0), Vec(2.0,2.0)))
       @test meta(m_meta) === (boundingbox = Rectangle(Point(1.0,1.0), Vec(2.0,2.0)),)
       @test metafree(m_meta) === m
       @test propertynames(m_meta) == (:mesh, :boundingbox)
   end
end

@testset "embedding MetaT" begin
    @testset "MetaT{Polygon}" begin
        polys = [Polygon(rand(Point2f, 20)) for i in 1:10]
        multipol = MultiPolygon(polys)
        pnames = [randstring(4) for i in 1:10]
        numbers = LinRange(0.0, 1.0, 10)
        bin = rand(Bool, 10)
        # create a polygon
        poly = MetaT(polys[1], name = pnames[1], value = numbers[1], category = bin[1])
        # create a MultiPolygon with the right type & meta information!
        multipoly = MetaT(multipol, name = pnames, value = numbers, category = bin)
        @test multipoly isa MetaT
        @test poly isa MetaT

        @test Meshes.getcolumn(poly, :name) == pnames[1]
        @test Meshes.getcolumn(multipoly, :name) == pnames

        meta_p = MetaT(polys[1], boundingbox = Rectangle(Point(0,0), Vec(2,2)))
        @test meta_p.boundingbox === Rectangle(Point(0,0), Vec(2,2))
        @test Meshes.metafree(meta_p) == polys[1]
        @test Meshes.metafree(poly) == polys[1]
        @test Meshes.metafree(multipoly) == multipol
        @test Meshes.meta(meta_p) == (boundingbox = Rectangle(Point(0,0), Vec(2, 2)),)
        @test Meshes.meta(poly) == (name = pnames[1], value = 0.0, category = bin[1])
        @test Meshes.meta(multipoly) == (name = pnames, value = numbers, category = bin)
    end

    @testset "MetaT{Point}" begin
        p = Point(1.1, 2.2)
        pm = MetaT(Point(1.1, 2.2); a=1, b=2)
        @test pm.meta === (a=1, b=2)
        @test pm.main === p
        @test propertynames(pm) == (:main, :a, :b)
        @test Meshes.metafree(pm) == p
        @test Meshes.meta(pm) == (a = 1, b = 2)
    end

    @testset "MetaT{MultiPoint}" begin
        ps = [Point2(x, x+1) for x in 1:5]
        @test ps isa AbstractVector
        mpm = MetaT(MultiPoint(ps); a=1, b=2)
        @test mpm.main == ps
        @test mpm.meta === (a=1, b=2)
        @test propertynames(mpm) == (:main, :a, :b)
        @test Meshes.metafree(mpm) == ps
        @test Meshes.meta(mpm) == (a = 1, b = 2)
    end

    @testset "MetaT{LineString}" begin
        linestring = MetaT(LineString(Point{2, Int}[(10, 10), (20, 20), (10, 40)]), a = 1, b = 2)
        @test linestring isa MetaT
        @test linestring.meta === (a = 1, b = 2)
        @test propertynames(linestring) == (:main, :a, :b)
        @test Meshes.metafree(linestring) == LineString(Point{2, Int}[(10, 10), (20, 20), (10, 40)])
        @test Meshes.meta(linestring) == (a = 1, b = 2)
    end

    @testset "MetaT{MultiLineString}" begin
        linestring1 = LineString(Point{2, Int}[(10, 10), (20, 20), (10, 40)])
        linestring2 = LineString(Point{2, Int}[(40, 40), (30, 30), (40, 20), (30, 10)])
        multilinestring = MultiLineString([linestring1, linestring2])
        multilinestringmeta = MetaT(MultiLineString([linestring1, linestring2]); boundingbox = Rectangle(Point2(1,1), Vec2(2,2)))
        @test multilinestringmeta isa MetaT
        @test multilinestringmeta.meta === (boundingbox = Rectangle(Point(1.0,1.0), Vec(2.0,2.0)),)
        @test multilinestringmeta.main == multilinestring
        @test propertynames(multilinestringmeta) == (:main, :boundingbox)
        @test Meshes.metafree(multilinestringmeta) == multilinestring
        @test Meshes.meta(multilinestringmeta) == (boundingbox = Rectangle(Point(1.0,1.0), Vec(2.0,2.0)),)
    end

    #=
    So mesh works differently for MetaT
    since `MetaT{Point}` not subtyped to `AbstractPoint`
    =#

   @testset "MetaT{Mesh}" begin
        @testset "per vertex attributes" begin
            points = rand(Point3, 8)
            tfaces = TetrahedronFace{Int}[(1, 2, 3, 4), (5, 6, 7, 8)]
            normals = rand(Vec3, 8)
            stress = LinRange(0, 1, 8)
            mesh_nometa = Mesh(points, tfaces)
            mesh = MetaT(mesh_nometa, normals = normals, stress = stress)

            @test hasproperty(mesh, :stress)
            @test hasproperty(mesh, :normals)
            @test mesh.stress == stress
            @test mesh.normals == normals
            @test Meshes.faces(mesh.main) == tfaces
            @test propertynames(mesh) == (:main, :normals, :stress)
        end
    end
end

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
        x = connect([1, 2, 3, 4, 5, 6, 7, 8], Point2)
        tetrahedra = connect(x, NSimplex{4})
        @test tetrahedra == [Tetrahedron(x[1], x[2], x[3], x[4])]

        @testset "matrix non-copy point views" begin
            # point in row
            points = [1 2; 1 4; 66 77]
            comparison = [Point(1, 2), Point(1, 4), Point(66, 77)]
            @test connect(points, Point{2}) == comparison
            # point in column
            points = [1 1 66; 2 4 77]
            # huh, reinterpret array doesn't seem to like `==`
            @test all(((a,b),)-> a==b, zip(connect(points, Point{2}), comparison))
        end
    end

    @testset "face views" begin
        numbers = [1, 2, 3, 4, 5, 6]
        points = connect(numbers, Point{2})
        faces = connect([1, 2, 3], TriangleFace)
        triangles = connect(points, faces)
        @test triangles == [Triangle(Point(1, 2), Point(3, 4), Point(5, 6))]
        x = Point3(1,1,1)
        triangles = connect([x], [TriangleFace(1, 1, 1)])
        @test triangles == [Triangle(x, x, x)]
        points = connect([1, 2, 3, 4, 5, 6, 7, 8], Point2)
        faces = connect([1, 2, 3, 4], TetrahedronFace{Int})
        triangles = connect(points, faces)
        @test triangles == [Tetrahedron(points...)]
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

        faces = [LineFace(1, 2)
        , LineFace(3, 4)]
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

        faces = [LineFace(1, 2), LineFace(3, 4)]
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
        mesh = Mesh([x], [TriangleFace(1, 1, 1)])
        @test mesh == [Triangle(x, x, x)]

        points = connect([1, 2, 3, 4, 5, 6, 7, 8], Point2)
        sfaces = connect([1, 2, 3, 4], TetrahedronFace{Int})
        mesh = Mesh(points, sfaces)
        @test mesh == [Tetrahedron(points...)]

        points = rand(Point3f, 8)
        tfaces = [TriangleFace(1, 2, 3), TriangleFace(5, 6, 7)]
        mesh = Mesh(points, tfaces)
        @test mesh isa PlainMesh

        t = Tesselation(Rectangle(Point2f(0,0), Vec2f(2,2)), (30, 30))
        m = Meshes.mesh(t, pointtype=Point2f, facetype=QuadFace)
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
    primitive = HyperSphere(Point3f(0,0,0), 1.0f0)
    @test ndims(primitive) === 3
    mesh = triangle_mesh(primitive)
    @test decompose(Point3f, mesh) isa Vector{Point3f}
    @test decompose(Point3f, primitive) isa Vector{Point3f}

    primitive = Rectangle(Point2(0, 0), Vec2(1, 1))
    mesh = triangle_mesh(primitive)

    @test decompose(Point2f, mesh) isa Vector{Point2f}
    @test decompose(Point{2,Int}, primitive) isa Vector{Point{2,Int}}

    primitive = Rectangle(Point3(0,0,0), Vec3(1,1,1))
    triangle_mesh(primitive)

    points = decompose(Point2f, HyperSphere(Point2f(0, 0), 1.0f0))
    m = Meshes.mesh(points)
    @test coordinates(m) == points

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

@testset "modifying meta" begin
    xx = rand(10)
    points = rand(Point3f, 10)
    m = Meshes.Mesh(meta(points, xx=xx), TriangleFace[(1,2,3), (3,4,5)])
    color = rand(10)
    m = pointmeta(m; color=color)

    @test hasproperty(m, :xx)
    @test hasproperty(m, :color)
    @test_throws ErrorException Meshes.MetaType(Simplex)
    @test_throws ErrorException Meshes.MetaFree(Simplex)

    @test m.xx === xx
    @test m.color === color

    m, colpopt = Meshes.pop_pointmeta(m, :color)
    m, xxpopt = Meshes.pop_pointmeta(m, :xx)

    @test propertynames(m) == (:position,)
    @test colpopt === color
    @test xxpopt === xx

    @testset "creating meta" begin
        x = Point3f[(1,3,4)]
        # no meta gets added, so should stay the same
        @test meta(x) === x
        @test meta(x, value=[1]).position === x
    end
    pos = Point2f[(10, 2)]
    m = Mesh(meta(pos, uv=[Vec2f(1, 1)]), [TriangleFace(1, 1, 1)])
    @test m.position === pos
end

@testset "mesh conversion" begin
    s = HyperSphere(Point3(0,0,0), 1.0)
    m = Meshes.mesh(s, pointtype=Point3)
    @test m isa Mesh{3,Float64}
    @test coordinates(m) isa Vector{Point3}
    @test Meshes.faces(m) isa Vector{TriangleFace}
    points1 = coordinates(m)
    points2 = decompose(Point3, m)
    @test coordinates.(points1) ≈ coordinates.(points2)

    m = Meshes.mesh(s, pointtype=Point3f)
    tmesh = triangle_mesh(m)
    points1 = coordinates(tmesh)
    points2 = decompose(Point3f, tmesh)
    @test tmesh isa PlainMesh
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

@testset "Offsetintegers" begin
    x = 1
    @test Meshes.raw(x) isa Int64
    @test Meshes.value(x) == x

    x = ZeroIndex(1)
    @test eltype(x) == Int64

    x = OffsetInteger{0}(1)
    @test typeof(x) == OffsetInteger{0,Int64}

    x1 = OffsetInteger{0}(2)
    @test Meshes.pure_max(x, x1) == x1
    @test promote_rule(typeof(x), typeof(x1)) == OffsetInteger{0,Int64}
    x2 = 1
    @test promote_rule(typeof(x2), typeof(x1)) == Int64
    @test Base.to_index(x1) == 2
    @test -(x1) == OffsetInteger{0,Int64}(-2)
    @test abs(x1) == OffsetInteger{0,Int64}(2)
    @test +(x, x1) == OffsetInteger{0,Int64}(3)
    @test *(x, x1) == OffsetInteger{0,Int64}(2)
    @test -(x, x1) == OffsetInteger{0,Int64}(-1)
    #test for /
    @test div(x, x1) == OffsetInteger{0,Int64}(0)
    @test !==(x, x1)
    @test !>=(x, x1)
    @test <=(x, x1)
    @test !>(x, x1)
    @test <(x, x1)
end

@testset "MetaT and heterogeneous data" begin
    ls = [LineString([Point(i*1.0, (i+1)^2/6), Point(i*0.86,i+5.0), Point(i/3, i/7)]) for i in 1:10]
    mls = MultiLineString([LineString([Point(i+1.0, (i)^2/6), Point(i*0.75,i+8.0), Point(i/2.5, i/6.79)]) for i in 5:10])
    poly = Polygon(Point{2, Int}[(40, 40), (20, 45), (45, 30), (40, 40)])
    geom = [ls..., mls, poly]
    prop = Any[(country_states = "India$(i)", rainfall = (i*9)/2) for i in 1:11]
    push!(prop, (country_states = 12, rainfall = 1000))   # a pinch of heterogeneity

    feat = [MetaT(i, j) for (i,j) = zip(geom, prop)]
    sa = meta_table(feat)

    @test nameof(eltype(feat)) == :MetaT
    @test eltype(sa) === MetaT{Any,(:country_states, :rainfall),Tuple{Any,Float64}}
    @test propertynames(sa) === (:main, :country_states, :rainfall)
    @test getproperty(sa, :country_states) isa Array{Any}
    @test getproperty(sa, :main) == geom

    @test Meshes.getnamestypes(typeof(feat[1])) ==
    (LineString{2,Float64,Point{2,Float64},Base.ReinterpretArray{Meshes.Ngon{2,Float64,2,Point{2,Float64}},1,Tuple{Point{2,Float64},Point{2,Float64}},TupleView{Tuple{Point{2,Float64},Point{2,Float64}},2,1,Array{Point{2,Float64},1}}}},
    (:country_states, :rainfall), Tuple{String,Float64})

    @test StructArrays.createinstance(typeof(feat[1]), LineString([Point(1.0, (2)^2/6), Point(1*0.86,6.0), Point(1/3, 1/7)]), "Mumbai", 100) isa typeof(feat[1])

    @test Base.getindex(feat[1], 1) isa Line
    @test Base.size(feat[1]) == (2,)
end

@testset "Tests from GeometryTypes" begin
    include("geometrytypes.jl")
end

end
