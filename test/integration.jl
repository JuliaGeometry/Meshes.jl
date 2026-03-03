@testitem "integral" setup = [Setup] begin
  # Bezier Curve
  curve = BezierCurve([cart(t, sin(t), 0) for t in range(-π, π, length=361)])
  function integrand(p)
    ux = ustrip(coords(p).x)
    (1 / sqrt(1 + cos(ux)^2)) * u"Ω"
  end
  solution = T(2π) * u"Ω*m"
  @test integral(integrand, curve, n=10) ≈ solution rtol = 1e-2

  # Box 1D
  a = T(π)
  box = Box(cart(0), cart(a))
  function integrand(p)
    x₁ = only(ustrip.(to(p)))
    √(a^2 - x₁^2) * u"A"
  end
  solution = T(π) * a^2 / 4 * u"A*m"
  @test integral(integrand, box, n=10) ≈ solution rtol = 1e-3

  # Box 2D
  a = T(π)
  box = Box(cart(0, 0), cart(a, a))
  function integrand(p)
    x₁, x₂ = ustrip.(to(p))
    (√(a^2 - x₁^2) + √(a^2 - x₂^2)) * u"A"
  end
  solution = 2a * (T(π) * a^2 / 4) * u"A*m^2"
  @test integral(integrand, box, n=10) ≈ solution rtol = 1e-3

  # Box 3D
  a = T(π)
  box = Box(cart(0, 0, 0), cart(a, a, a))
  function integrand(p)
    x₁, x₂, x₃ = ustrip.(to(p))
    (√(a^2 - x₁^2) + √(a^2 - x₂^2) + √(a^2 - x₃^2)) * u"A"
  end
  solution = 3a^2 * (T(π) * a^2 / 4) * u"A*m^3"
  @test integral(integrand, box, n=10) ≈ solution rtol = 1e-3

  # Ball 2D
  origin = cart(0, 0)
  radius = T(2.8)
  ball = Ball(origin, radius)
  function integrand(p)
    r = ustrip(u"m", norm(to(p)))
    exp(-r^2) * u"A"
  end
  solution = (T(π) - T(π) * exp(-radius^2)) * u"A*m^2"
  @test integral(integrand, ball, n=10) ≈ solution rtol = 1e-3

  # Ellipsoid
  origin = cart(0, 0, 0)
  R = r₁ = r₂ = r₃ = T(4.1)
  ellipsoid = Ellipsoid((r₁, r₂, r₃), origin)
  function integrand(p)
    x, y, z = ustrip.(u"m", to(p))
    (z^2) * u"A"
  end
  solution = (T(4π) * R^4 / 3) * u"A*m^2"
  @test integral(integrand, ellipsoid, n=10) ≈ solution rtol = 1e-3

  # Disk
  center = cart(1, 2, 3)
  normal = vector(1 / 2, 1 / 2, sqrt(2) / 2)
  plane = Plane(center, normal)
  radius = T(2.5)
  disk = Disk(plane, radius)
  function integrand(p)
    offset = p - center
    r = ustrip(u"m", norm(offset))
    exp(-r^2) * u"A"
  end
  solution = (T(π) - T(π) * exp(-radius^2)) * u"A*m^2"
  @test integral(integrand, disk, n=10) ≈ solution rtol = 1e-3

  # Circle
  center = cart(1, 2, 3)
  normal = vector(1 / 2, 1 / 2, sqrt(2) / 2)
  plane = Plane(center, normal)
  radius = T(4.4)
  circle = Circle(plane, radius)
  function integrand(p)
    offset = p - center
    r = ustrip(u"m", norm(offset))
    exp(-r^2) * u"A"
  end
  solution = T(2π) * radius * exp(-radius^2) * u"A*m"
  @test integral(integrand, circle, n=10) ≈ solution rtol = 1e-3
end

@testitem "localintegral" setup = [Setup] begin
  # TODO
end
