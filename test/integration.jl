@testitem "integral" setup = [Setup] begin
  # Ray
  ray = Ray(cart(0, 0, 0), vector(1, 1, 1))
  @test integral(ray) do p
    r = ustrip(norm(to(p)))
    exp(-r^2) * u"A"
  end ≈ √(T(π)) / 2 * u"A*m" rtol = 1e-3

  # Line
  line = Line(cart(0, 0, 0), cart(1, 1, 1))
  @test integral(line) do p
    r = ustrip(norm(to(p)))
    exp(-r^2) * u"A"
  end ≈ √(T(π)) * u"A*m" rtol = 1e-3

  # BezierCurve
  bezier = BezierCurve([cart(t, sin(t), 0) for t in range(-π, π, length=361)])
  @test integral(bezier) do p
    x = ustrip(coords(p).x)
    (1 / √(1 + cos(x)^2)) * u"Ω"
  end ≈ T(2π) * u"Ω*m" rtol = 1e-2

  # Plane
  plane = Plane(cart(0, 0, 0), vector(0, 0, 1))
  @test integral(plane) do p
    r = ustrip(norm(to(p)))
    exp(-r^2) * u"A"
  end ≈ T(π) * u"A*m^2" rtol = 1e-1

  # Box 1D
  a = T(π)
  box = Box(cart(0), cart(a))
  @test integral(box) do p
    x = ustrip(coords(p).x)
    √(a^2 - x^2) * u"A"
  end ≈ T(π) * a^2 / 4 * u"A*m" rtol = 1e-3

  # Box 2D
  box = Box(cart(0, 0), cart(a, a))
  @test integral(box) do p
    x, y = ustrip.(to(p))
    (√(a^2 - x^2) + √(a^2 - y^2)) * u"A"
  end ≈ 2a * (T(π) * a^2 / 4) * u"A*m^2" rtol = 1e-2

  # Box 3D
  box = Box(cart(0, 0, 0), cart(a, a, a))
  @test integral(box) do p
    x, y, z = ustrip.(to(p))
    (√(a^2 - x^2) + √(a^2 - y^2) + √(a^2 - z^2)) * u"A"
  end ≈ 3a^2 * (T(π) * a^2 / 4) * u"A*m^3" rtol = 1e-3

  # integral that is exactly zero doesn't hang
  box = Box(cart(0, 0, 0), cart(a, a, a))
  @test integral(box) do p
    x, y, z = ustrip.(to(p))
    (cos(x) + cos(y) + cos(z)) * u"A"
  end ≈ zero(T) * u"A*m^3" atol = 1e-3 * u"A*m^3"

  # Ball 2D
  r = T(2.8)
  ball = Ball(cart(0, 0), r)
  @test integral(ball) do p
    r = ustrip(norm(to(p)))
    exp(-r^2) * u"A"
  end ≈ (T(π) - T(π) * exp(-r^2)) * u"A*m^2" rtol = 1e-3

  # Ellipsoid
  R = r₁ = r₂ = r₃ = T(4.1)
  ellip = Ellipsoid((r₁, r₂, r₃), cart(0, 0, 0))
  @test integral(ellip) do p
    z = ustrip(coords(p).z)
    z^2 * u"A"
  end ≈ (T(4π) * R^4 / 3) * u"A*m^2" rtol = 1e-2

  # Disk
  r = T(2.5)
  o = cart(1, 2, 3)
  n = vector(1 / 2, 1 / 2, √(2) / 2)
  disk = Disk(Plane(o, n), r)
  @test integral(disk) do p
    r = ustrip(norm(p - o))
    exp(-r^2) * u"A"
  end ≈ (T(π) - T(π) * exp(-r^2)) * u"A*m^2" rtol = 1e-3

  # Circle
  r = T(4.4)
  o = cart(1, 2, 3)
  n = vector(1 / 2, 1 / 2, √(2) / 2)
  circle = Circle(Plane(o, n), r)
  @test integral(circle) do p
    r = ustrip(norm(p - o))
    exp(-r^2) * u"A"
  end ≈ T(2π) * r * exp(-r^2) * u"A*m" rtol = 1e-3

  # Cylinder
  h = T(8.5)u"m"
  ρ = T(1.3)u"m"
  a = cart(0, 0, 0)
  b = cart(0u"m", 0u"m", h)
  cyl = Cylinder(a, b, ρ)
  @test integral(cyl) do p
    c = convert(Cylindrical, coords(p))
    ρ = c.ρ
    ϕ = c.ϕ
    z = c.z
    ρ^(-1) * (ρ + ϕ * u"m" + z) * u"A"
  end ≈ ((T(π) * h * ρ^2) + (T(π) * h^2 * ρ) + (T(2π) * T(π) * u"m" * h * ρ)) * u"A" rtol = 1e-3

  # CylinderSurface
  h = T(8.5)u"m"
  ρ = T(1.3)u"m"
  a = cart(0, 0, 0)
  b = cart(0u"m", 0u"m", h)
  cylsurf = CylinderSurface(a, b, ρ)
  @test integral(cylsurf) do p
    c = convert(Cylindrical, coords(p))
    ρ = c.ρ
    ϕ = c.ϕ
    z = c.z
    ρ^(-1) * (ρ + ϕ * u"m" + z) * u"A"
  end ≈ let
    A1 = (T(2π) * h * ρ) + (T(π) * ρ^2) + (T(π) * u"m" * ρ * T(2π))
    A2 = (T(π) * ρ^2) + (T(π) * u"m" * ρ * T(2π))
    A3 = (T(2π) * h * ρ) + (2T(π)^2 * u"m" * h) + (T(π) * h^2)
    (A1 + A2 + A3) * u"A"
  end rtol = 1e-3

  # Cone
  r = T(2.5)u"m"
  h = T(3.5)u"m"
  base = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), r)
  apex = cart(0u"m", 0u"m", h)
  cone = Cone(base, apex)
  @test integral(cone) do p
    T(1) * u"A"
  end ≈ (T(π) * r^2 * h / 3) * u"A" rtol = 1e-3

  # ConeSurface
  r = T(2.5)u"m"
  h = T(3.5)u"m"
  base = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), r)
  apex = cart(0u"m", 0u"m", h)
  conesurf = ConeSurface(base, apex)
  @test integral(conesurf) do p
    T(1) * u"A"
  end ≈ ((T(π) * r^2) + (T(π) * r * hypot(h, r))) * u"A" rtol = 1e-3

  # Frustum
  r = T(2.5)u"m"
  h = T(3.5)u"m"
  disk1 = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), r)
  disk2 = Disk(Plane(cart(0.0u"m", 0.0u"m", h / 2), vector(0, 0, 1)), r / 2)
  frustum = Frustum(disk1, disk2)
  @test integral(frustum) do p
    T(1) * u"A"
  end ≈ (T(7) / T(8)) * (T(π) * r^2 * h / T(3)) * u"A" rtol = 1e-3

  # FrustumSurface
  rbot = T(2.5)u"m"
  rtop = T(1.25)u"m"
  height = T(2π) * u"m"
  origin = cart(0, 0, 0)
  normal = vector(0, 0, 1)
  planebot = Plane(origin, normal)
  diskbot = Disk(planebot, rbot)
  centertop = cart(0.0u"m", 0.0u"m", height / 2)
  planetop = Plane(centertop, normal)
  disktop = Disk(planetop, rtop)
  frustumsurf = FrustumSurface(diskbot, disktop)
  @test integral(frustumsurf) do p
    T(1) * u"A"
  end ≈ let
    A1 = T(π) * rbot * hypot(height, rbot)
    A2 = T(π) * rtop * hypot(height / 2, rtop)
    A3 = T(π) * rtop^2
    A4 = T(π) * rbot^2
    (A1 - A2 + A3 + A4) * u"A"
  end rtol = 1e-3

  # Segment
  ϕ = 7T(pi) / 6
  θ = T(pi) / 3
  a = cart(0, 0, 0)
  b = cart(sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))
  ka = T(7.1)
  kb = T(4.6)
  seg = Segment(a, b)
  @test integral(seg) do p
    r = ustrip(norm(to(p)))
    exp(r * log(ka) + (1 - r) * log(kb)) * u"A"
  end ≈ ((ka - kb) / (log(ka) - log(kb))) * u"A*m" rtol = 1e-3

  # Rope
  a = cart(0, 0, 0)
  b = cart(1, 0, 0)
  c = cart(1, 1, 0)
  d = cart(1, 1, 1)
  rope = Rope(a, b, c, d)
  @test integral(rope) do p
    x, y, z = ustrip.(to(p))
    (x + 2y + 3z) * u"A"
  end ≈ T(7.0)u"A*m" rtol = 1e-3

  # Ring
  a = cart(0, 0, 0)
  b = cart(1, 0, 0)
  c = cart(1, 1, 0)
  d = cart(1, 1, 1)
  ring = Ring(a, b, c, d, c, b)
  @test integral(ring) do p
    x, y, z = ustrip.(to(p))
    (x + 2y + 3z) * u"A"
  end ≈ T(14.0)u"A*m" rtol = 1e-3

  # Pgon
  a, b, c, z = T(0.4), T(0.6), T(1.0), T(0.0)
  outer = [(z, z), (c, z), (c, c), (z, c)]
  hole = [(a, a), (a, b), (b, b), (b, a)]
  poly = Pgon([outer, hole])
  @test integral(poly) do p
    x = ustrip(coords(p).x)
    2x * u"A"
  end ≈ (c^2 - (b - a) * (b^2 - a^2)) * u"A*m^2" rtol = 1e-3

  # Triangle
  a = cart(0, 0, 0)
  b = cart(1, 0, 0)
  c = cart(0, 1, 0)
  tri = Triangle(a, b, c)
  @test integral(tri) do p
    x, y, z = ustrip.(to(p))
    (x + 2y + 3z) * u"A"
  end ≈ T(0.5) * u"A*m^2" rtol = 1e-3

  # Quadrangle
  quad = Quadrangle(cart(-1.0, 0.0), cart(-1.0, 1.0), cart(1.0, 1.0), cart(1.0, 0.0))
  @test integral(quad) do p
    r = ustrip(norm(to(p)))
    exp(-r^2) * u"A"
  end ≈ T(π) * T(0.8427007929497149)^2 / 2 * u"A*m^2" rtol = 1e-3 # erf(1) = 0.8427007929497149

  # Tetrahedron
  a = cart(0, 0, 0)
  b = cart(1, 0, 0)
  c = cart(0, 1, 0)
  d = cart(0, 0, 1)
  tetra = Tetrahedron(a, b, c, d)
  @test integral(tetra) do p
    x, y, z = ustrip.(to(p))
    (x + 2y + 3z) * u"A"
  end ≈ T(0.25) * u"A*m^3" rtol = 1e-3

  # Hexahedron
  a = T(π)
  box = Box(cart(0, 0, 0), cart(a, a, a))
  hexa = convert(Hexahedron, box)
  @test integral(hexa) do p
    x, y, z = ustrip.(to(p))
    (√(a^2 - x^2) + √(a^2 - y^2) + √(a^2 - z^2)) * u"A"
  end ≈ 3a^2 * (T(π) * a^2 / 4) * u"A*m^3" rtol = 1e-3

  # Multi
  box = Box(cart(0, 0), cart(1, 1))
  ball = Ball(cart(5, 5), T(1))
  multi = Multi([box, ball])
  fmulti(p) = sum(to(p))
  @test integral(fmulti, multi) ≈ integral(fmulti, box) + integral(fmulti, ball)

  # GeometrySet
  box = Box(cart(0, 0), cart(1, 1))
  ball = Ball(cart(5, 5), T(1))
  gset = GeometrySet([box, ball])
  fgset(p) = sum(to(p))
  @test integral(fgset, gset) ≈ integral(fgset, box) + integral(fgset, ball)

  # SimpleMesh
  points = [cart(0, 0), cart(1, 0), cart(0, 1), cart(1, 1), cart(0.25, 0.5), cart(0.75, 0.5)]
  tris = connect.([(1, 5, 3), (4, 6, 2)], Triangle)
  quads = connect.([(1, 2, 6, 5), (4, 3, 5, 6)], Quadrangle)
  mesh = SimpleMesh(points, [tris; quads])
  fmesh(p) = T(1) * u"A"
  @test integral(fmesh, mesh) ≈ sum(integral(fmesh, elem) for elem in mesh)
  @test integral(fmesh, mesh) ≈ T(1) * u"A * m^2" rtol = 1e-3

  # Grid
  grid = cartgrid(10, 10)
  fgrid(p) = T(1) * u"A"
  @test integral(fgrid, grid) ≈ 100 * integral(fgrid, first(grid))
end
