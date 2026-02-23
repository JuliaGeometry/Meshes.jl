@testitem "Jacobian" setup = [Setup] begin
  method = FiniteDifference(T(1e-6))

  # 1D geometry (segment): constant tangent vector
  a = cart(T(1), T(2))
  b = cart(T(3), T(6))
  seg = Segment(a, b)
  Jmid = jacobian(seg, (T(0.5),), method)
  Jleft = jacobian(seg, (T(0.0),), method)
  Jright = jacobian(seg, (T(1.0),), method)
  @test length(Jmid) == 1
  @test Jmid[1] ≈ b - a
  @test Jleft[1] ≈ b - a
  @test Jright[1] ≈ b - a

  # invalid number of parametric coordinates
  @test_throws ArgumentError jacobian(seg, (T(0.1), T(0.2)), method)

  # 2D geometry (triangle): barycentric map derivatives
  v₁ = cart(T(0), T(0), T(0))
  v₂ = cart(T(2), T(0), T(0))
  v₃ = cart(T(0), T(3), T(1))
  tri = Triangle(v₁, v₂, v₃)
  Jtri = jacobian(tri, (T(0.2), T(0.3)), method)
  @test length(Jtri) == 2
  @test Jtri[1] ≈ v₂ - v₁
  @test Jtri[2] ≈ v₃ - v₁

  # 3D geometry (box): axis-aligned derivatives
  pmin = cart(T(0), T(0), T(0))
  pmax = cart(T(2), T(3), T(4))
  box = Box(pmin, pmax)
  Jbox = jacobian(box, (T(0.4), T(0.6), T(0.7)), method)
  sides = pmax - pmin
  @test length(Jbox) == 3
  @test Jbox[1] ≈ vector(sides[1], T(0), T(0))
  @test Jbox[2] ≈ vector(T(0), sides[2], T(0))
  @test Jbox[3] ≈ vector(T(0), T(0), sides[3])

  # non-constant Jacobian
  # q(u,v) = (1-u)(1-v)c₀₀ + u(1-v)c₀₁ + (1-u)v c₁₀ + uv c₁₁
  # ∂q/∂u = -(1-v)c₀₀ + (1-v)c₀₁ - v c₁₀ + v c₁₁ = (1-v)(c₀₁-c₀₀) + v(c₁₁-c₁₀)
  # ∂q/∂v = -(1-u)c₀₀ - u c₀₁ + (1-u)c₁₀ + u c₁₁ = (1-u)(c₁₀-c₀₀) + u(c₁₁-c₀₁)
  c₀₀ = cart(T(0), T(0))
  c₀₁ = cart(T(2), T(0))
  c₁₁ = cart(T(3), T(3))
  c₁₀ = cart(T(0), T(2))
  quad = Quadrangle(c₀₀, c₀₁, c₁₁, c₁₀)

  # corner (0,0): ∂u = c₀₁ - c₀₀, ∂v = c₁₀ - c₀₀
  J00 = jacobian(quad, (T(0), T(0)), method)
  @test J00[1] ≈ c₀₁ - c₀₀
  @test J00[2] ≈ c₁₀ - c₀₀

  # corner (1,1): ∂u = c₁₁ - c₁₀, ∂v = c₁₁ - c₀₁
  J11 = jacobian(quad, (T(1), T(1)), method)
  @test J11[1] ≈ c₁₁ - c₁₀
  @test J11[2] ≈ c₁₁ - c₀₁

  # center (0.5, 0.5): derivatives are averages
  J_center = jacobian(quad, (T(0.5), T(0.5)), method)
  @test J_center[1] ≈ ((c₀₁ - c₀₀) + (c₁₁ - c₁₀)) / 2
  @test J_center[2] ≈ ((c₁₀ - c₀₀) + (c₁₁ - c₀₁)) / 2
end

@testitem "Differential" setup = [Setup] begin
  method = FiniteDifference(T(1e-6))

  # line element
  a = cart(T(1), T(2))
  b = cart(T(3), T(6))
  seg = Segment(a, b)
  @test differential(seg, (T(0.25),), method) ≈ measure(seg)

  # surface element
  v₁ = cart(T(0), T(0), T(0))
  v₂ = cart(T(2), T(0), T(0))
  v₃ = cart(T(0), T(3), T(1))
  tri = Triangle(v₁, v₂, v₃)
  expectedarea = norm((v₂ - v₁) × (v₃ - v₁))
  @test differential(tri, (T(0.3), T(0.2)), method) ≈ expectedarea

  # volume element
  pmin = cart(T(0), T(0), T(0))
  pmax = cart(T(2), T(3), T(4))
  box = Box(pmin, pmax)
  @test differential(box, (T(0.1), T(0.2), T(0.3)), method) ≈ measure(box)

  # non-constant differential element (bilinear quadrangle)
  c₀₀ = cart(T(0), T(0))
  c₀₁ = cart(T(2), T(0))
  c₁₁ = cart(T(3), T(3))
  c₁₀ = cart(T(0), T(2))
  quad = Quadrangle(c₀₀, c₀₁, c₁₁, c₁₀)

  # differential varies with position due to non-constant Jacobian
  d00 = differential(quad, (T(0), T(0)), method)
  d11 = differential(quad, (T(1), T(1)), method)
  dcenter = differential(quad, (T(0.5), T(0.5)), method)

  # compute expected values from cross product magnitudes
  expected_d00 = norm((c₀₁ - c₀₀) × (c₁₀ - c₀₀))
  expected_d11 = norm((c₁₁ - c₁₀) × (c₁₁ - c₀₁))
  expected_dcenter = norm(((c₀₁ - c₀₀) + (c₁₁ - c₁₀)) / 2 × ((c₁₀ - c₀₀) + (c₁₁ - c₀₁)) / 2)

  @test d00 ≈ expected_d00
  @test d11 ≈ expected_d11
  @test dcenter ≈ expected_dcenter
end
