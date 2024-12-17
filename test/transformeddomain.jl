@testitem "TransformedDomain" setup = [Setup] begin
  pset = PointSet(randpoint2(10))
  trans1 = Scale(T(2), T(2))
  trans2 = Translate(T(4), T(4))
  tpset1 = TransformedDomain(pset, trans1)
  @test parent(tpset1) == pset
  @test Meshes.transform(tpset1) == trans1
  tpset2 = TransformedDomain(tpset1, trans2)
  @test parent(tpset2) == pset
  @test Meshes.transform(tpset2) == (trans1 → trans2)

  pset = PointSet(randpoint2(10))
  trans = Identity()
  tpset = TransformedDomain(pset, trans)
  @test nelements(tpset) == nelements(pset)
  @test element(tpset, 1) == element(pset, 1)
  @test tpset == pset
  trans = Translate(T(2), T(2))
  tpset = TransformedDomain(pset, trans)
  @test nelements(tpset) == nelements(pset)
  @test element(tpset, 1) == trans(element(pset, 1))
  @test tpset == trans(pset)

  t1 = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
  t2 = Triangle(cart(1, 1), cart(2, 1), cart(1, 2))
  gset = GeometrySet([t1, t2])
  trans = Identity()
  tgset = TransformedDomain(gset, trans)
  @test nelements(tgset) == nelements(gset)
  @test element(tgset, 1) == element(gset, 1)
  @test tgset == gset
  trans = Translate(T(2), T(2))
  tgset = TransformedDomain(gset, trans)
  @test nelements(tgset) == nelements(gset)
  @test element(tgset, 1) == trans(element(gset, 1))
  @test tgset == trans(gset)

  grid = cartgrid(10, 10)
  trans = Identity()
  tgrid = TransformedDomain(grid, trans)
  @test nelements(tgrid) == nelements(grid)
  @test element(tgrid, 1) == element(grid, 1)
  @test tgrid == grid
  trans = Translate(T(2), T(2))
  tgrid = TransformedDomain(grid, trans)
  @test nelements(tgrid) == nelements(grid)
  @test element(tgrid, 1) == trans(element(grid, 1))
  @test tgrid == trans(grid)

  grid = cartgrid(10, 10)
  subgrid = view(grid, 1:10)
  trans = Identity()
  tsubgrid = TransformedDomain(subgrid, trans)
  @test nelements(tsubgrid) == nelements(subgrid)
  @test element(tsubgrid, 1) == element(subgrid, 1)
  @test tsubgrid == subgrid
  trans = Translate(T(2), T(2))
  tsubgrid = TransformedDomain(subgrid, trans)
  @test nelements(tsubgrid) == nelements(subgrid)
  @test element(tsubgrid, 1) == trans(element(subgrid, 1))
  @test tsubgrid == trans(subgrid)

  # transforms that change the CRS
  pset = PointSet(latlon(0, 0), latlon(0, 45), latlon(45, 0))
  trans = Proj(Mercator)
  tpset = TransformedDomain(pset, trans)
  @test manifold(tpset) === 𝔼{2}
  @test crs(tpset) <: Mercator
  @test nelements(tpset) == nelements(pset)
  @test element(tpset, 1) == trans(element(pset, 1))
  @test tpset == trans(pset)
end
