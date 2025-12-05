@testitem "TransformedMesh" setup = [Setup] begin
  grid = cartgrid(10, 10)
  trans = Identity()
  tgrid = TransformedGrid(grid, trans)
  @test nelements(tgrid) == nelements(grid)
  @test element(tgrid, 1) == element(grid, 1)
  @test tgrid == grid
  trans = Translate(T(2), T(2))
  tgrid = TransformedGrid(grid, trans)
  @test nelements(tgrid) == nelements(grid)
  @test element(tgrid, 1) == trans(element(grid, 1))
  @test tgrid == trans(grid)

  grid = cartgrid(10, 10)
  subgrid = view(grid, 1:10)
  trans = Identity()
  tsubgrid = TransformedSubGrid(subgrid, trans)
  @test nelements(tsubgrid) == nelements(subgrid)
  @test element(tsubgrid, 1) == element(subgrid, 1)
  @test tsubgrid == subgrid
  trans = Translate(T(2), T(2))
  tsubgrid = TransformedSubGrid(subgrid, trans)
  @test nelements(tsubgrid) == nelements(subgrid)
  @test element(tsubgrid, 1) == trans(element(subgrid, 1))
  @test tsubgrid == trans(subgrid)
end
