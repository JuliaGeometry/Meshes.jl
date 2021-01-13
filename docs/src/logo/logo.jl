using Luxor, Random

Random.seed!(2021)

function meshes(s)
  Δ = s/600 # scale factor, original design was 600units
  scale(Δ)
  # background
  squircle(O, 285, 285, :clip, rt=0.2)
  sethue("slateblue4")
  paint()

  # build point array
  pts = polysample(box(O, 600, 600, vertices=true), 16)
  for pt in ngon(O, 100, 3, π/6, vertices=true)
    for θ in 0:π/8:6π
      for i in 160:40:265
        push!(pts, pt  + polar(i, θ))
      end
    end
  end

  # make triangulation
  tris = polytriangulate(pts)

  # draw triangles
  setline(.5)
  for tri in tris
    sethue("grey")
    poly(tri, :stroke, close=true)
  end

  # clipping for julia circles
  @layer begin
    for pt in ngon(O + (0, 35), 150, 3, π/6, vertices=true)
      circlepath(pt, 110, :path)
    end
    sethue("white")
    setline(3)
    strokepreserve()
    clip()
    setline(1)
    for tri in tris
      sethue([Luxor.julia_red, Luxor.julia_green, Luxor.julia_purple][rand(1:end)])
      poly(tri, :fillpreserve, close=true)
      sethue("white")
      strokepath()
    end
    sethue("white")
    circle.(pts, 3, :fill)
    clipreset()
  end

  # draw nodes
  sethue("black")
  circle.(pts, 3, :fill)

  #outline
  setline(4)
  sethue("black")
  squircle(O, 285, 285, :stroke, rt=0.2)
end

function logo(s, fname)
  Drawing(s, s, fname)
  origin()
  meshes(s)
  finish()
  preview()
end

function logotext(w, h, fname)
  Drawing(w, h, fname)
  origin()
  table = Table([h], [h, w - h])
  @layer begin
    translate(table[1])
    meshes(h)
  end
  @layer begin
    translate(table[2])
    background("white")
    sethue("black")
    # find all fonts available on Linux with `fc-list | -f 2 -d ":"`
    fontface("Julius Sans One")
    fontsize(h/2.5)
    text("Meshes.jl", halign=:center, valign=:middle)
  end
  @layer begin
    translate(table[1])
    meshes(h)
  end
  finish()
  preview()
end

for ext in [".svg", ".png"]
  logo(240, "../assets/logo"*ext)
  logotext(700, 200, "../assets/logo-text"*ext)
end
