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
    g = Ray(latlon(0, 0), vector(1, 1, 1))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Line(latlon(0, 0), latlon(1, 1))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Plane(latlon(0, 0), vector(1, 0, 0), vector(0, 1, 0))
    @test crs(g) <: LatLon{WGS84Latest}
    g = BezierCurve(latlon(0, 0), latlon(1, 1), latlon(0, 2))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Box(latlon(0, 180), latlon(45, 90))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Ball(latlon(0, 0), T(1))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Sphere(latlon(0, 0), T(1))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Ellipsoid((T(3), T(2), T(1)), latlon(0, 0))
    @test crs(g) <: LatLon{WGS84Latest}
    p = Plane(latlon(0, 0), vector(1, 0, 0), vector(0, 1, 0))
    g = Disk(p, T(2))
    @test crs(g) <: LatLon{WGS84Latest}
    p = Plane(latlon(0, 0), vector(1, 0, 0), vector(0, 1, 0))
    g = Circle(p, T(2))
    @test crs(g) <: LatLon{WGS84Latest}
    b = Plane(latlon(90, 0), vector(1, 0, 0), vector(0, 1, 0))
    t = Plane(latlon(-90, 0), vector(1, 0, 0), vector(0, 1, 0))
    g = Cylinder(b, t, T(5))
    @test crs(g) <: LatLon{WGS84Latest}
    b = Plane(latlon(-90, 0), vector(1, 0, 0), vector(0, 1, 0))
    t = Plane(latlon(90, 0), vector(1, 0, 0), vector(0, 1, 0))
    g = CylinderSurface(b, t, T(5))
    @test crs(g) <: LatLon{WGS84Latest}
    g = ParaboloidSurface(latlon(0, 0), T(1), T(2))
    @test crs(g) <: LatLon{WGS84Latest}
    p = Plane(latlon(-90, 0), vector(1, 0, 0), vector(0, 1, 0))
    d = Disk(p, T(2))
    a = latlon(90, 0)
    g = Cone(d, a)
    @test crs(g) <: LatLon{WGS84Latest}
    p = Plane(latlon(-90, 0), vector(1, 0, 0), vector(0, 1, 0))
    d = Disk(p, T(2))
    a = latlon(90, 0)
    g = ConeSurface(d, a)
    @test crs(g) <: LatLon{WGS84Latest}
    pb = Plane(latlon(-90, 0), vector(1, 0, 0), vector(0, 1, 0))
    db = Disk(pb, T(1))
    pt = Plane(latlon(90, 0), vector(1, 0, 0), vector(0, 1, 0))
    dt = Disk(pt, T(2))
    g = Frustum(db, dt)
    @test crs(g) <: LatLon{WGS84Latest}
    pb = Plane(latlon(-90, 0), vector(1, 0, 0), vector(0, 1, 0))
    db = Disk(pb, T(1))
    pt = Plane(latlon(90, 0), vector(1, 0, 0), vector(0, 1, 0))
    dt = Disk(pt, T(2))
    g = FrustumSurface(db, dt)
    @test crs(g) <: LatLon{WGS84Latest}
    g = Torus(latlon(0, 0), vector(1, 0, 0), T(2), T(1))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Segment(latlon(0, 0), latlon(1, 1))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Rope(latlon(0, 0), latlon(0, 1), latlon(1, 0))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Ring(latlon(0, 0), latlon(0, 1), latlon(1, 0))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Triangle(latlon(0, 0), latlon(0, 1), latlon(1, 0))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Quadrangle(latlon(0, 0), latlon(0, 1), latlon(1, 1), latlon(1, 0))
    @test crs(g) <: LatLon{WGS84Latest}
    g = PolyArea(latlon(0, 0), latlon(0, 1), latlon(1, 0))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Tetrahedron(latlon(0, 0), latlon(0, 90), latlon(0, -90), latlon(90, 0))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Hexahedron(
      latlon(0, 45),
      latlon(0, 135),
      latlon(0, -135),
      latlon(0, -45),
      latlon(1, 45),
      latlon(1, 135),
      latlon(1, -135),
      latlon(1, -45)
    )
    @test crs(g) <: LatLon{WGS84Latest}
    g = Pyramid(latlon(0, 45), latlon(0, 135), latlon(0, -135), latlon(0, -45), latlon(90, 0))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Wedge(latlon(0, 0), latlon(0, 90), latlon(0, -90), latlon(1, 0), latlon(1, 90), latlon(1, -90))
    @test crs(g) <: LatLon{WGS84Latest}
    g = Multi([latlon(0, 0), latlon(1, 1)])
    @test crs(g) <: LatLon{WGS84Latest}
    t1 = Triangle(latlon(0, 0), latlon(0, 1), latlon(1, 0))
    t2 = Triangle(latlon(1, 1), latlon(1, 2), latlon(2, 1))
    d = GeometrySet([t1, t2])
    @test crs(d) <: LatLon{WGS84Latest}
    d = PointSet([latlon(0, 0), latlon(1, 1)])
    @test crs(d) <: LatLon{WGS84Latest}
    d = CartesianGrid((10, 10, 10), latlon(0, 0), (T(1), T(1), T(1)))
    @test crs(d) <: LatLon{WGS84Latest}
    p = latlon.([(0, 0), (0, 1), (1, 0), (1, 1), (0.5, 0.5)])
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    @test crs(d) <: LatLon{WGS84Latest}
    d = CylindricalTrajectory([latlon(0, 0), latlon(1, 1), latlon(0, 2)])
    @test crs(d) <: LatLon{WGS84Latest}
  end
end
