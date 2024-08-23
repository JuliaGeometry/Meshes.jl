@testitem "rand" begin
  p = rand(Point)
  @test p isa Point
  @test crs(p) <: Cartesian3D
  @test Meshes.lentype(p) === Meshes.Met{Float64}
  p = rand(Point, crs=Cartesian2D)
  @test p isa Point
  @test crs(p) <: Cartesian2D
  @test Meshes.lentype(p) === Meshes.Met{Float64}

  r = rand(Ray)
  @test r isa Ray
  @test crs(r) <: Cartesian3D
  @test Meshes.lentype(r) === Meshes.Met{Float64}
  r = rand(Ray, crs=Cartesian2D)
  @test r isa Ray
  @test crs(r) <: Cartesian2D
  @test Meshes.lentype(r) === Meshes.Met{Float64}

  l = rand(Line)
  @test l isa Line
  @test crs(l) <: Cartesian3D
  @test Meshes.lentype(l) === Meshes.Met{Float64}
  l = rand(Line, crs=Cartesian2D)
  @test l isa Line
  @test crs(l) <: Cartesian2D
  @test Meshes.lentype(l) === Meshes.Met{Float64}

  p = rand(Plane)
  @test p isa Plane
  @test crs(p) <: Cartesian3D
  @test Meshes.lentype(p) === Meshes.Met{Float64}
  p = rand(Plane, crs=LatLon)
  @test p isa Plane
  @test crs(p) <: LatLon
  @test Meshes.lentype(p) === Meshes.Met{Float64}

  b = rand(BezierCurve)
  @test b isa BezierCurve
  @test crs(b) <: Cartesian3D
  @test Meshes.lentype(b) === Meshes.Met{Float64}
  b = rand(BezierCurve, crs=Cartesian2D)
  @test b isa BezierCurve
  @test crs(b) <: Cartesian2D
  @test Meshes.lentype(b) === Meshes.Met{Float64}

  b = rand(Box)
  @test b isa Box
  @test crs(b) <: Cartesian3D
  @test Meshes.lentype(b) === Meshes.Met{Float64}
  b = rand(Box, crs=Cartesian2D)
  @test b isa Box
  @test crs(b) <: Cartesian2D
  @test Meshes.lentype(b) === Meshes.Met{Float64}

  b = rand(Ball)
  @test b isa Ball
  @test crs(b) <: Cartesian3D
  @test Meshes.lentype(b) === Meshes.Met{Float64}
  b = rand(Ball, crs=Cartesian2D)
  @test b isa Ball
  @test crs(b) <: Cartesian2D
  @test Meshes.lentype(b) === Meshes.Met{Float64}

  s = rand(Sphere)
  @test s isa Sphere
  @test crs(s) <: Cartesian3D
  @test Meshes.lentype(s) === Meshes.Met{Float64}
  s = rand(Sphere, crs=Cartesian2D)
  @test s isa Sphere
  @test crs(s) <: Cartesian2D
  @test Meshes.lentype(s) === Meshes.Met{Float64}

  e = rand(Ellipsoid)
  @test e isa Ellipsoid
  @test crs(e) <: Cartesian3D
  @test Meshes.lentype(e) === Meshes.Met{Float64}
  e = rand(Ellipsoid, crs=LatLon)
  @test e isa Ellipsoid
  @test crs(e) <: LatLon
  @test Meshes.lentype(e) === Meshes.Met{Float64}

  d = rand(Disk)
  @test d isa Disk
  @test crs(d) <: Cartesian3D
  @test Meshes.lentype(d) === Meshes.Met{Float64}
  d = rand(Disk, crs=LatLon)
  @test d isa Disk
  @test crs(d) <: LatLon
  @test Meshes.lentype(d) === Meshes.Met{Float64}

  c = rand(Circle)
  @test c isa Circle
  @test crs(c) <: Cartesian3D
  @test Meshes.lentype(c) === Meshes.Met{Float64}
  c = rand(Circle, crs=LatLon)
  @test c isa Circle
  @test crs(c) <: LatLon
  @test Meshes.lentype(c) === Meshes.Met{Float64}

  c = rand(Cylinder)
  @test c isa Cylinder
  @test crs(c) <: Cartesian3D
  @test Meshes.lentype(c) === Meshes.Met{Float64}
  c = rand(Cylinder, crs=LatLon)
  @test c isa Cylinder
  @test crs(c) <: LatLon
  @test Meshes.lentype(c) === Meshes.Met{Float64}

  c = rand(CylinderSurface)
  @test c isa CylinderSurface
  @test crs(c) <: Cartesian3D
  @test Meshes.lentype(c) === Meshes.Met{Float64}
  c = rand(CylinderSurface, crs=LatLon)
  @test c isa CylinderSurface
  @test crs(c) <: LatLon
  @test Meshes.lentype(c) === Meshes.Met{Float64}

  p = rand(ParaboloidSurface)
  @test p isa ParaboloidSurface
  @test crs(p) <: Cartesian3D
  @test Meshes.lentype(p) === Meshes.Met{Float64}
  p = rand(ParaboloidSurface, crs=LatLon)
  @test p isa ParaboloidSurface
  @test crs(p) <: LatLon
  @test Meshes.lentype(p) === Meshes.Met{Float64}

  c = rand(Cone)
  @test c isa Cone
  @test crs(c) <: Cartesian3D
  @test Meshes.lentype(c) === Meshes.Met{Float64}
  c = rand(Cone, crs=LatLon)
  @test c isa Cone
  @test crs(c) <: LatLon
  @test Meshes.lentype(c) === Meshes.Met{Float64}

  c = rand(ConeSurface)
  @test c isa ConeSurface
  @test crs(c) <: Cartesian3D
  @test Meshes.lentype(c) === Meshes.Met{Float64}
  c = rand(ConeSurface, crs=LatLon)
  @test c isa ConeSurface
  @test crs(c) <: LatLon
  @test Meshes.lentype(c) === Meshes.Met{Float64}

  f = rand(Frustum)
  @test f isa Frustum
  @test crs(f) <: Cartesian3D
  @test Meshes.lentype(f) === Meshes.Met{Float64}
  f = rand(Frustum, crs=LatLon)
  @test f isa Frustum
  @test crs(f) <: LatLon
  @test Meshes.lentype(f) === Meshes.Met{Float64}

  f = rand(FrustumSurface)
  @test f isa FrustumSurface
  @test crs(f) <: Cartesian3D
  @test Meshes.lentype(f) === Meshes.Met{Float64}
  f = rand(FrustumSurface, crs=LatLon)
  @test f isa FrustumSurface
  @test crs(f) <: LatLon
  @test Meshes.lentype(f) === Meshes.Met{Float64}

  t = rand(Torus)
  @test t isa Torus
  @test crs(t) <: Cartesian3D
  @test Meshes.lentype(t) === Meshes.Met{Float64}
  t = rand(Torus, crs=LatLon)
  @test t isa Torus
  @test crs(t) <: LatLon
  @test Meshes.lentype(t) === Meshes.Met{Float64}

  s = rand(Segment)
  @test s isa Segment
  @test crs(s) <: Cartesian3D
  @test Meshes.lentype(s) === Meshes.Met{Float64}
  s = rand(Segment, crs=Cartesian2D)
  @test s isa Segment
  @test crs(s) <: Cartesian2D
  @test Meshes.lentype(s) === Meshes.Met{Float64}

  r = rand(Rope)
  @test r isa Rope
  @test crs(r) <: Cartesian3D
  @test Meshes.lentype(r) === Meshes.Met{Float64}
  r = rand(Rope, crs=Cartesian2D)
  @test r isa Rope
  @test crs(r) <: Cartesian2D
  @test Meshes.lentype(r) === Meshes.Met{Float64}

  r = rand(Ring)
  @test r isa Ring
  @test crs(r) <: Cartesian3D
  @test Meshes.lentype(r) === Meshes.Met{Float64}
  r = rand(Ring, crs=Cartesian2D)
  @test r isa Ring
  @test crs(r) <: Cartesian2D
  @test Meshes.lentype(r) === Meshes.Met{Float64}

  NGONS = [Triangle, Quadrangle, Pentagon, Hexagon, Heptagon, Octagon, Nonagon, Decagon]
  for NGON in NGONS
    n = rand(NGON)
    @test n isa NGON
    @test crs(n) <: Cartesian3D
    @test Meshes.lentype(n) === Meshes.Met{Float64}
    n = rand(NGON, crs=Cartesian2D)
    @test n isa NGON
    @test crs(n) <: Cartesian2D
    @test Meshes.lentype(n) === Meshes.Met{Float64}
  end

  p = rand(PolyArea)
  @test p isa PolyArea
  @test crs(p) <: Cartesian3D
  @test Meshes.lentype(p) === Meshes.Met{Float64}
  p = rand(PolyArea, crs=Cartesian2D)
  @test p isa PolyArea
  @test crs(p) <: Cartesian2D
  @test Meshes.lentype(p) === Meshes.Met{Float64}

  t = rand(Tetrahedron)
  @test t isa Tetrahedron
  @test crs(t) <: Cartesian3D
  @test Meshes.lentype(t) === Meshes.Met{Float64}
  t = rand(Tetrahedron, crs=LatLon)
  @test t isa Tetrahedron
  @test crs(t) <: LatLon
  @test Meshes.lentype(t) === Meshes.Met{Float64}

  h = rand(Hexahedron)
  @test h isa Hexahedron
  @test crs(h) <: Cartesian3D
  @test Meshes.lentype(h) === Meshes.Met{Float64}
  h = rand(Hexahedron, crs=LatLon)
  @test h isa Hexahedron
  @test crs(h) <: LatLon
  @test Meshes.lentype(h) === Meshes.Met{Float64}

  p = rand(Pyramid)
  @test p isa Pyramid
  @test crs(p) <: Cartesian3D
  @test Meshes.lentype(p) === Meshes.Met{Float64}
  p = rand(Pyramid, crs=LatLon)
  @test p isa Pyramid
  @test crs(p) <: LatLon
  @test Meshes.lentype(p) === Meshes.Met{Float64}

  w = rand(Wedge)
  @test w isa Wedge
  @test crs(w) <: Cartesian3D
  @test Meshes.lentype(w) === Meshes.Met{Float64}
  w = rand(Wedge, crs=LatLon)
  @test w isa Wedge
  @test crs(w) <: LatLon
  @test Meshes.lentype(w) === Meshes.Met{Float64}

  # vector of random geometries
  ps = rand(Point, 10)
  @test eltype(ps) <: Point
end
