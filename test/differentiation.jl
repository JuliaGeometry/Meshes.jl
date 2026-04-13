@testitem "Jacobian" setup = [Setup] begin
  import DifferentiationInterface as DI
  import Enzyme

  # 1D geometry (segment): constant tangent vector
  a = cart(T(1), T(2))
  b = cart(T(3), T(6))
  seg = Segment(a, b)
  Jl = jacobian(seg, (T(0.0),))
  Jr = jacobian(seg, (T(1.0),))
  Jm = jacobian(seg, (T(0.5),))
  @test Jl[1] ≈ b - a
  @test Jr[1] ≈ b - a
  @test Jm[1] ≈ b - a

  Jm = jacobian(seg, (T(0.5),); dbackend=DI.AutoEnzyme(mode=Enzyme.Forward))
  @test Jm[1] ≈ b - a
  Jm = jacobian(seg, (T(0.5),); dbackend=DI.AutoEnzyme(mode=Enzyme.Reverse))
  @test Jm[1] ≈ b - a

  # invalid number of parametric coordinates
  @test_throws ArgumentError jacobian(seg, (T(0.1), T(0.2)))

  # 2D geometry (triangle): barycentric map derivatives
  p₁ = cart(T(0), T(0), T(0))
  p₂ = cart(T(2), T(0), T(0))
  p₃ = cart(T(0), T(3), T(1))
  tri = Triangle(p₁, p₂, p₃)
  J = jacobian(tri, (T(0.2), T(0.3)))
  @test J[1] ≈ p₂ - p₁
  @test J[2] ≈ p₃ - p₁
  # 3D geometry (box): axis-aligned derivatives
  pmin = cart(T(0), T(0), T(0))
  pmax = cart(T(2), T(3), T(4))
  box = Box(pmin, pmax)
  J = jacobian(box, (T(0.4), T(0.6), T(0.7)))
  sides = pmax - pmin
  @test J[1] ≈ vector(sides[1], T(0), T(0))
  @test J[2] ≈ vector(T(0), sides[2], T(0))
  @test J[3] ≈ vector(T(0), T(0), sides[3])

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
  J₀₀ = jacobian(quad, (T(0), T(0)))
  @test J₀₀[1] ≈ c₀₁ - c₀₀
  @test J₀₀[2] ≈ c₁₀ - c₀₀

  # corner (1,1): ∂u = c₁₁ - c₁₀, ∂v = c₁₁ - c₀₁
  J₁₁ = jacobian(quad, (T(1), T(1)))
  @test J₁₁[1] ≈ c₁₁ - c₁₀
  @test J₁₁[2] ≈ c₁₁ - c₀₁

  # center (0.5, 0.5): derivatives are averages
  Jc = jacobian(quad, (T(0.5), T(0.5)))
  @test Jc[1] ≈ ((c₀₁ - c₀₀) + (c₁₁ - c₁₀)) / 2
  @test Jc[2] ≈ ((c₁₀ - c₀₀) + (c₁₁ - c₀₁)) / 2
end

@testitem "Differential" setup = [Setup] begin
  # line element
  a = cart(T(1), T(2))
  b = cart(T(3), T(6))
  seg = Segment(a, b)
  @test differential(seg, (T(0.25),)) ≈ length(seg)

  # surface element
  p₁ = cart(T(0), T(0), T(0))
  p₂ = cart(T(2), T(0), T(0))
  p₃ = cart(T(0), T(3), T(1))
  tri = Triangle(p₁, p₂, p₃)
  @test differential(tri, (T(0.3), T(0.2))) ≈ 2area(tri)

  # volume element
  pmin = cart(T(0), T(0), T(0))
  pmax = cart(T(2), T(3), T(4))
  box = Box(pmin, pmax)
  @test differential(box, (T(0.1), T(0.2), T(0.3))) ≈ volume(box)

  # non-constant differential element (bilinear quadrangle)
  c₀₀ = cart(T(0), T(0))
  c₀₁ = cart(T(2), T(0))
  c₁₁ = cart(T(3), T(3))
  c₁₀ = cart(T(0), T(2))
  quad = Quadrangle(c₀₀, c₀₁, c₁₁, c₁₀)

  # differential varies with position due to non-constant Jacobian
  d₀₀ = differential(quad, (T(0), T(0)))
  d₁₁ = differential(quad, (T(1), T(1)))
  dc = differential(quad, (T(0.5), T(0.5)))

  # compute expected values from cross product magnitudes
  @test d₀₀ ≈ norm((c₀₁ - c₀₀) × (c₁₀ - c₀₀))
  @test d₁₁ ≈ norm((c₁₁ - c₁₀) × (c₁₁ - c₀₁))
  @test dc ≈ norm(((c₀₁ - c₀₀) + (c₁₁ - c₁₀)) / 2 × ((c₁₀ - c₀₀) + (c₁₁ - c₀₁)) / 2)

  # https://github.com/JuliaGeometry/Meshes.jl/issues/1336
  b = BezierCurve([cart(t, sin(t), 0) for t in range(-π, π, length=361)])
  t = T(0.7878788f0)
  @test !isnan(differential(b, (t,)))
end
