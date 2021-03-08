SUITE["intersections"] = BenchmarkGroup(["microbenchmark", "intersections"]) # add tags for eventual filtering

add_benchmark!(SUITE["intersections"], [
  quote
    s1 = Segment((0.0,0.0), (1.0,0.0))
    s2 = Segment((0.5,0.0), (2.0,0.0))
    s1 ∩ s2
  end
])
