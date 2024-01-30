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

# helper function to print a large indexable collection
# in multiple lines with a given tabulation
function printelms(io::IO, elms, tab="")
  N = length(elms)
  I, J = N > 10 ? (5, N - 4) : (N - 1, N)
  for i in 1:I
    println(io, "$(tab)├─ $(elms[i])")
  end
  if N > 10
    println(io, "$(tab)⋮")
  end
  for i in J:(N - 1)
    println(io, "$(tab)├─ $(elms[i])")
  end
  print(io, "$(tab)└─ $(elms[N])")
end

# helper function to print a large iterable
# calling the printelms function
printitr(io::IO, itr, tab="") = printelms(io, collect(itr), tab)

# helper function to print the polygons vertices
function printverts(io::IO, verts)
  ioctx = IOContext(io, :compact => true)
  if length(verts) > 3
    join(ioctx, (first(verts), "...", last(verts)), ", ")
  else
    join(ioctx, verts, ", ")
  end
end

# helper function to print the view indices
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

function printfields(io, obj; compact=false)
  fnames = fieldnames(typeof(obj))
  if compact
    ioctx = IOContext(io, :compact => true)
    vals = map(fnames) do field
      val = getfield(obj, field)
      str = repr(val, context=ioctx)
      "$field: $str"
    end
    join(io, vals, ", ")
  else
    len = length(fnames)
    for (i, field) in enumerate(fnames)
      div = i == len ? "\n└─ " : "\n├─ "
      val = getfield(obj, field)
      str = repr(val, context=io)
      print(io, "$div$field: $str")
    end
  end
end
