# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# helper function to print the actual name of
# the geometry inside a deep type hierarchy
prettyname(geom) = prettyname(typeof(geom))
function prettyname(G::Type)
  name = string(G)
  name = replace(name, r"{.*" => "")
  replace(name, r".+\." => "")
end

# helper function to print a large iterator
# in multiple lines with a given tabulation
function io_lines(itr, tab="")
  vec = collect(itr)
  N = length(vec)
  I, J = N > 10 ? (5, N - 4) : (N - 1, N)
  lines = [
    ["$(tab)├─ $(vec[i])" for i in 1:I]
    (N > 10 ? ["$(tab)⋮"] : [])
    ["$(tab)├─ $(vec[i])" for i in J:(N - 1)]
    "$(tab)└─ $(vec[N])"
  ]
  join(lines, "\n")
end

function printverts(io::IO, verts)
  ioctx = IOContext(io, :compact => true)
  if length(verts) > 3
    join(ioctx, (first(verts), "...", last(verts)), ", ")
  else
    join(ioctx, verts, ", ")
  end
end
