@testitem "TranformedDomain" setup = [Setup] begin
  pset = PointSet(randpoint2(10))
  trans = Identity()
  tpset = TranformedDomain(pset, trans)
  @test nelements(tpset) == nelements(pset)
  @test element(tpset, 1) == element(pset, 1)
  @test tpset == pset
  pts = randpoint2(10)
  pset = PointSet(pts)
  trans = Translate(T(2), T(2))
  tpset = TranformedDomain(pset, trans)
  @test nelements(tpset) == nelements(pset)
  @test element(tpset, 1) == trans(element(pset, 1))
  @test tpset == PointSet(trans.(pts))
end
