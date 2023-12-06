# Input/Output

The [GeoIO.jl](https://github.com/JuliaEarth/GeoIO.jl) package can be
used to load/save mesh data in various different formats, including
VTK, GIS, PLY and many other formats. The package provides two functions
`GeoIO.load` and `GeoIO.save`, which are self-explanatory:

```julia
geotable = GeoIO.load("data.vtr")
```

```julia
GeoIO.save("data.vtu", geotable)
```

Please check the [Geospatial Data Science with Julia](https://juliaearth.github.io/geospatial-data-science-with-julia)
book for more information.
