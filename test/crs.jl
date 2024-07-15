@testset "CRS" begin
  @testset "Projected" begin
    g = merc(1, 1)
    @test crs(g) <: Mercator{WGS84Latest}
    g = Ray(merc(0, 0), vector(1, 1))
    @test crs(g) <: Mercator{WGS84Latest}
    g = Line(merc(0, 0), merc(1, 1))
    @test crs(g) <: Mercator{WGS84Latest}
    g = BezierCurve(merc(0, 0), merc(1, 1), merc(2, 0))
    @test crs(g) <: Mercator{WGS84Latest}
    g = Box(merc(0, 0), merc(1, 1))
    @test crs(g) <: Mercator{WGS84Latest}
    g = Ball(merc(0, 0), T(1))
    @test crs(g) <: Mercator{WGS84Latest}
    g = Sphere(merc(0, 0), T(1))
    @test crs(g) <: Mercator{WGS84Latest}
    g = Segment(merc(0, 0), merc(1, 1))
    @test crs(g) <: Mercator{WGS84Latest}
    g = Rope(merc(0, 0), merc(1, 0), merc(0, 1))
    @test crs(g) <: Mercator{WGS84Latest}
    g = Ring(merc(0, 0), merc(1, 0), merc(0, 1))
    @test crs(g) <: Mercator{WGS84Latest}
    g = Triangle(merc(0, 0), merc(1, 0), merc(0, 1))
    @test crs(g) <: Mercator{WGS84Latest}
    g = Quadrangle(merc(0, 0), merc(1, 0), merc(1, 1), merc(0, 1))
    @test crs(g) <: Mercator{WGS84Latest}
    g = PolyArea(merc(0, 0), merc(1, 0), merc(0, 1))
    @test crs(g) <: Mercator{WGS84Latest}
    g = Multi([merc(0, 0), merc(1, 1)])
    @test crs(g) <: Mercator{WGS84Latest}
    t1 = Triangle(merc(0, 0), merc(1, 0), merc(0, 1))
    t2 = Triangle(merc(1, 1), merc(2, 1), merc(1, 2))
    d = GeometrySet([t1, t2])
    @test crs(d) <: Mercator{WGS84Latest}
    d = PointSet([merc(0, 0), merc(1, 1)])
    @test crs(d) <: Mercator{WGS84Latest}
    d = CartesianGrid((10, 10), merc(0, 0), (T(1), T(1)))
    @test crs(d) <: Mercator{WGS84Latest}
    p = merc.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    @test crs(d) <: Mercator{WGS84Latest}
  end
  
  @testset "Geographic" begin
    g = latlon(1, 1)
    @test crs(g) <: LatLon{WGS84Latest}
    g = Ray(latlon(0, 0), vector(1, 1))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Line(latlon(0, 0), latlon(1, 1))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Plane(latlon(90, 0), vector(1, 0, 0), vector(0, 1, 0))
    @test crs(g) <: LatLon{WGS84Latest}
    g = BezierCurve(latlon(0, 0), latlon(1, 1), latlon(0, 2))
    @test crs(g) <: Mercator{WGS84Latest}
  end
end
