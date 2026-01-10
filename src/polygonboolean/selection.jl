# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# select segments based on fill information and operation
function _selectsegments(fillsegs, operation)
  P = eltype(vertices(first(fillsegs)))

  selected = Vector{Tuple{P,P}}()
  fills = Vector{UInt8}()
  seen = Set{Tuple{P,P}}()

  for i in 1:length(fillsegs)
    a, b = vertices(fillsegs[i])
    a, b = a < b ? (a, b) : (b, a)
    if (a, b) ∈ seen
      continue
    end

    bits = _filltobits(fillsegs, i)

    # check if this segment is included based on operation
    isfilledabove = _filled(operation, bits, true)
    isfilledbelow = _filled(operation, bits, false)

    if isfilledabove ⊻ isfilledbelow
      push!(seen, (a, b))
      push!(selected, (a, b))
      push!(fills, bits)
    end
  end
  selected, fills
end
