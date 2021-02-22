suite["points"] = BenchmarkGroup(["microbenchmark", "points"]) # add tags for eventual filtering

add_benchmark!(suite["points"], [
  :(Point(0., 1.)),
  :(Point([0., 1.])),
  :(Point([0., 1])),
  :(Point(0., convert(Float64, 1))),
])
