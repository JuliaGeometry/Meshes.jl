@testitem "boundarypoints" setup = [Setup] begin
  p = cart(0, 0)
  @test boundarypoints(p) == [cart(0, 0)]

  sphere = Sphere(cart(0, 0), T(1))
  points = boundarypoints(sphere)
  @test all(∈(sphere), points)

  ball = Ball(cart(0, 0), T(1))
  points = boundarypoints(ball)
  @test all(∈(boundary(ball)), points)

  verts = [cart(0, 0), cart(1, 1)]
  segment = Segment(verts...)
  points = boundarypoints(segment)
  @test points == verts

  verts = [cart(0, 0), cart(1, 0), cart(1, 1)]
  triangle = Triangle(verts...)
  points = boundarypoints(triangle)
  @test points == verts

  verts = [cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1)]
  quadrangle = Quadrangle(verts...)
  points = boundarypoints(quadrangle)
  @test points == verts

  tri = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
  quad = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  multi = Multi([tri, quad])
  points = boundarypoints(multi)
  @test points == [boundarypoints(tri); boundarypoints(quad)]

  box = Box(cart(0, 0), cart(1, 1))
  trans = Translate(T(1), T(2))
  tbox = TransformedGeometry(box, trans)
  points = boundarypoints(tbox)
  @test points == trans.(boundarypoints(box))

  box = Box(latlon(0, 0), latlon(45, 45))
  trans = Proj(Mercator)
  tbox = TransformedGeometry(box, trans)
  points = boundarypoints(tbox)
  @test points == vertices(discretize(boundary(tbox)))
end
