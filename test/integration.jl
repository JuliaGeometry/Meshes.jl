@testitem "integral" setup = [Setup] begin
  # Bezier Curve
  bezier = BezierCurve([cart(t, sin(t), 0) for t in range(-π, π, length=361)])
  function funbezier(p)
    ux = ustrip(coords(p).x)
    (1 / sqrt(1 + cos(ux)^2)) * u"Ω"
  end
  solution = T(2π) * u"Ω*m"
  @test integral(funbezier, bezier, n=10) ≈ solution rtol = 1e-2

  # Box 1D
  a = T(π)
  box = Box(cart(0), cart(a))
  function funbox1(p)
    x₁ = only(ustrip.(to(p)))
    √(a^2 - x₁^2) * u"A"
  end
  solution = T(π) * a^2 / 4 * u"A*m"
  @test integral(funbox1, box, n=10) ≈ solution rtol = 1e-3

  # Box 2D
  a = T(π)
  box = Box(cart(0, 0), cart(a, a))
  function funbox2(p)
    x₁, x₂ = ustrip.(to(p))
    (√(a^2 - x₁^2) + √(a^2 - x₂^2)) * u"A"
  end
  solution = 2a * (T(π) * a^2 / 4) * u"A*m^2"
  @test integral(funbox2, box, n=10) ≈ solution rtol = 1e-3

  # Box 3D
  a = T(π)
  box = Box(cart(0, 0, 0), cart(a, a, a))
  function funbox3(p)
    x₁, x₂, x₃ = ustrip.(to(p))
    (√(a^2 - x₁^2) + √(a^2 - x₂^2) + √(a^2 - x₃^2)) * u"A"
  end
  solution = 3a^2 * (T(π) * a^2 / 4) * u"A*m^3"
  @test integral(funbox3, box, n=10) ≈ solution rtol = 1e-3

  # Ball 2D
  origin = cart(0, 0)
  radius = T(2.8)
  ball = Ball(origin, radius)
  function funball2(p)
    r = ustrip(u"m", norm(to(p)))
    exp(-r^2) * u"A"
  end
  solution = (T(π) - T(π) * exp(-radius^2)) * u"A*m^2"
  @test integral(funball2, ball, n=10) ≈ solution rtol = 1e-3

  # Ellipsoid
  origin = cart(0, 0, 0)
  R = r₁ = r₂ = r₃ = T(4.1)
  ellipsoid = Ellipsoid((r₁, r₂, r₃), origin)
  function funellips(p)
    x, y, z = ustrip.(u"m", to(p))
    (z^2) * u"A"
  end
  solution = (T(4π) * R^4 / 3) * u"A*m^2"
  @test integral(funellips, ellipsoid, n=10) ≈ solution rtol = 1e-3

  # Disk
  center = cart(1, 2, 3)
  normal = vector(1 / 2, 1 / 2, sqrt(2) / 2)
  plane = Plane(center, normal)
  radius = T(2.5)
  disk = Disk(plane, radius)
  function fundisk(p)
    offset = p - center
    r = ustrip(u"m", norm(offset))
    exp(-r^2) * u"A"
  end
  solution = (T(π) - T(π) * exp(-radius^2)) * u"A*m^2"
  @test integral(fundisk, disk, n=10) ≈ solution rtol = 1e-3

  # Circle
  center = cart(1, 2, 3)
  normal = vector(1 / 2, 1 / 2, sqrt(2) / 2)
  plane = Plane(center, normal)
  radius = T(4.4)
  circle = Circle(plane, radius)
  function funcircle(p)
    offset = p - center
    r = ustrip(u"m", norm(offset))
    exp(-r^2) * u"A"
  end
  solution = T(2π) * radius * exp(-radius^2) * u"A*m"
  @test integral(funcircle, circle, n=10) ≈ solution rtol = 1e-3

  # Cylinder
  h = T(8.5)u"m"
  ρ = T(1.3)u"m"
  a = cart(0, 0, 0)
  b = cart(0u"m", 0u"m", h)
  cyl = Cylinder(a, b, ρ)
  function funcylinder(p)
    c = convert(Cylindrical, coords(p))
    ρ = c.ρ
    φ = c.ϕ
    z = c.z
    ρ^(-1) * (ρ + φ * u"m" + z) * u"A"
  end
  solution = ((T(π) * h * ρ^2) + (T(π) * h^2 * ρ) + (T(2π) * T(π) * u"m" * h * ρ)) * u"A"
  @test integral(funcylinder, cyl, n=10) ≈ solution

  # CylinderSurface
  h = T(8.5)u"m"
  ρ = T(1.3)u"m"
  a = cart(0, 0, 0)
  b = cart(0u"m", 0u"m", h)
  cylsurf = CylinderSurface(a, b, ρ)
  function funcylsurf(p)
    c = convert(Cylindrical, coords(p))
    ρ = c.ρ
    φ = c.ϕ
    z = c.z
    ρ^(-1) * (ρ + φ * u"m" + z) * u"A"
  end
  solution = let
    A1 = (T(2π) * h * ρ) + (T(π) * ρ^2) + (T(π) * u"m" * ρ * T(2π))
    A2 = (T(π) * ρ^2) + (T(π) * u"m" * ρ * T(2π))
    A3 = (T(2π) * h * ρ) + (T(2π)^2 * u"m" * h) + (T(π) * h^2)
    (A1 + A2 + A3) * u"A"
  end
  @test_broken integral(funcylsurf, cylsurf) ≈ solution

  # Cone
  r = T(2.5)u"m"
  h = T(3.5)u"m"
  origin = cart(0, 0, 0)
  plane = Plane(origin, vector(0, 0, 1))
  base = Disk(plane, r)
  apex = cart(0u"m", 0u"m", h)
  cone = Cone(base, apex)
  funcone(p) = T(1.0)u"A"
  solution = (T(π) * r^2 * h / 3) * u"A"
  @test integral(funcone, cone) ≈ solution

  # ConeSurface
  r = T(2.5)u"m"
  h = T(3.5)u"m"
  origin = cart(0, 0, 0)
  plane = Plane(origin, vector(0, 0, 1))
  base = Disk(plane, r)
  apex = cart(0u"m", 0u"m", h)
  conesurf = ConeSurface(base, apex)
  funconesurf(p) = T(1.0)u"A"
  solution = ((T(π) * r^2) + (T(π) * r * hypot(h, r))) * u"A"
  @test integral(funconesurf, conesurf) ≈ solution

  # Frustum
  r = T(2.5)u"m"
  h = T(3.5)u"m"
  origin = cart(0, 0, 0)
  normal = vector(0, 0, 1)
  midpoint = cart(0.0u"m", 0.0u"m", h / 2)
  base = Disk(Plane(origin, normal), r)
  disk = Disk(Plane(midpoint, normal), r / 2)
  frustum = Frustum(base, disk)
  funfrustum(p) = T(1.0)u"A"
  solution = (T(7) / T(8)) * (T(π) * r^2 * h / T(3)) * u"A"
  @test integral(funfrustum, frustum) ≈ solution

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
  funfrustumsurf(p) = T(1.0)u"A"
  solution = let
    A1 = T(π) * rbot * hypot(height, rbot)
    A2 = T(π) * rtop * hypot(height / 2, rtop)
    A3 = T(π) * rtop^2
    A4 = T(π) * rbot^2
    (A1 - A2 + A3 + A4) * u"A"
  end
  @test integral(funfrustumsurf, frustumsurf) ≈ solution
end
