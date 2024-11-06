@testitem "RegularCoarsening" setup = [Setup] begin
  # 2D grids
  grid = CartesianGrid(cart(0.0, 0.0), cart(10.0, 10.0), dims=(20, 20))
  tgrid = CartesianGrid(cart(0.0, 0.0), cart(10.0, 10.0), dims=(10, 10))
  @test coarsen(grid, RegularCoarsening(2)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  trgrid = convert(RectilinearGrid, tgrid)
  @test coarsen(rgrid, RegularCoarsening(2)) == trgrid
  sgrid = convert(StructuredGrid, grid)
  tsgrid = convert(StructuredGrid, tgrid)
  @test coarsen(sgrid, RegularCoarsening(2)) == tsgrid
  tfgrid = TransformedGrid(grid, Identity())
  @test coarsen(tfgrid, RegularCoarsening(2)) == coarsen(grid, RegularCoarsening(2))

  grid = CartesianGrid(cart(0.0, 0.0), cart(10.0, 10.0), dims=(20, 20))
  tgrid = CartesianGrid(cart(0.0, 0.0), cart(10.0, 10.0), dims=(10, 5))
  @test coarsen(grid, RegularCoarsening(2, 4)) == tgrid

  # non-multiple dimensions
  grid = CartesianGrid(cart(0, 0), cart(11, 17), dims=(11, 17))
  tgrid = CartesianGrid(cart(0, 0), cart(11, 17), dims=(3, 6))
  @test coarsen(grid, RegularCoarsening(5, 3)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  @test size(coarsen(rgrid, RegularCoarsening(5, 3)) ) == (3, 6)
  sgrid = convert(StructuredGrid, grid)
  @test size(coarsen(sgrid, RegularCoarsening(5, 3))) == (3, 6)
  tfgrid = TransformedGrid(grid, Identity())
  @test size(coarsen(tfgrid, RegularCoarsening(5, 3))) == (3, 6)

  # 3D grids
  grid = cartgrid(100, 100, 100)
  tgrid = CartesianGrid(minimum(grid), maximum(grid), dims=(50, 25, 20))
  @test coarsen(grid, RegularCoarsening(2, 4, 5)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  trgrid = convert(RectilinearGrid, tgrid)
  @test coarsen(rgrid, RegularCoarsening(2, 4, 5)) == trgrid
  sgrid = convert(StructuredGrid, grid)
  tsgrid = convert(StructuredGrid, tgrid)
  @test coarsen(sgrid, RegularCoarsening(2, 4, 5)) == tsgrid
  tfgrid = TransformedGrid(grid, Identity())
  @test coarsen(tfgrid, RegularCoarsening(2, 4, 5)) == coarsen(grid, RegularCoarsening(2, 4, 5))

  # non-multiple dimensions
  grid = CartesianGrid(cart(0, 0, 0), cart(11, 17, 23), dims=(11, 17, 23))
  tgrid = CartesianGrid(cart(0, 0, 0), cart(11, 17, 23), dims=(2, 4, 8))
  @test coarsen(grid, RegularCoarsening(7, 5, 3)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  @test size(coarsen(rgrid, RegularCoarsening(7, 5, 3)) ) == (2, 4, 8)
  sgrid = convert(StructuredGrid, grid)
  @test size(coarsen(sgrid, RegularCoarsening(7, 5, 3))) == (2, 4, 8)
  tfgrid = TransformedGrid(grid, Identity())
  @test size(coarsen(tfgrid, RegularCoarsening(7, 5, 3))) == (2, 4, 8)
end
