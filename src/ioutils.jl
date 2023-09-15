# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# helper function to print the actual name of
# the object type inside a deep type hierarchy
prettyname(obj) = prettyname(typeof(obj))
function prettyname(T::Type)
  name = string(T)
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

function printinds(io::IO, inds)
  print(io, "[")
  if length(inds) > 8
    join(io, @view(inds[1:4]), ", ")
    print(io, ", ..., ")
    join(io, @view(inds[(end - 3):end]), ", ")
  else
    join(io, inds, ", ")
  end
  print(io, "]")
end

printinds(io::IO, inds::AbstractRange) = print(io, inds)
