# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
# annotate fill segments with fill information using a sweep line algorithm
function _annotatefill!(fillsegs::FillSegments{S}, nsegsfirst::Int) where {S}
  isempty(fillsegs) && return
  # build event queue
  events = _buildevents(fillsegs, nsegsfirst)
  status = BinaryTrees.AVLTree{SegmentIndex{S},Bool}()

  while !BinaryTrees.isempty(events)
    node = BinaryTrees.minnode(events)
    p, evs = BinaryTrees.key(node), BinaryTrees.value(node)

    # sort: segments ending first, then by index
    sort!(evs, by=ev -> (SegmentIndex{S}(fillsegs, ev.ind), !ev.isstart))

    # remove ending segments
    for ev in evs
      !ev.isstart && BinaryTrees.delete!(status, SegmentIndex{S}(fillsegs, ev.ind))
    end

    # process start events
    startevs = filter(ev -> ev.isstart, evs)

    # pass 1: compute subjectfill (which where is inside of subject polygon relative to segments)
    for ev in startevs
      ind, issubject = ev.ind, ev.issubject
      below, _ = BinaryTrees.prevnext(status, SegmentIndex{S}(fillsegs, ind))

      # find closest segment from same polygon below
      belowind = nothing
      curr = below
      while !isnothing(curr)
        if BinaryTrees.value(curr) == issubject
          belowind = BinaryTrees.key(curr).ind
          break
        end
        curr, _ = BinaryTrees.prevnext(status, BinaryTrees.key(curr))
      end

      if isnothing(belowind) # exterior segment, set filled above, below = false
        setfills!(fillsegs, true, false, ind, issubject)
      else
        # otherwise reverse fill info from below segment
        fillabove, fillbelow = getfills(fillsegs, belowind, issubject)
        setfills!(fillsegs, fillbelow, fillabove, ind, issubject)
      end
      BinaryTrees.insert!(status, SegmentIndex{S}(fillsegs, ind), issubject)
    end

    # pass 2: compute fill information relative to the other polygon for remaining segments
    for ev in startevs
      ind, issubject = ev.ind, ev.issubject
      fillsegs.computed[ind] && continue

      below, _ = BinaryTrees.prevnext(status, SegmentIndex{S}(fillsegs, ind))
      if isnothing(below)
        # outside clip polygon
        setfills!(fillsegs, false, false, ind, !issubject)
      else
        # segment below
        bind = BinaryTrees.key(below).ind
        bissubject = BinaryTrees.value(below)
        # inside subject polygon?
        if issubject == bissubject
          # below segment is from subject polygon, use its clip fill
          inside = getfills(fillsegs, bind, !bissubject)[1]
        else
          # below segment is from clip polygon, use its subject fill
          inside = getfills(fillsegs, bind, bissubject)[1]
        end
        # mark segment fill accordingly
        setfills!(fillsegs, inside, inside, ind, !issubject)
      end
      # mark as computed
      _setcomputed!(fillsegs, ind)
    end
    # pass 3: swap fill information for duplicate segments
    for i in 1:length(startevs)
      evi = startevs[i]
      issubjecti = evi.issubject
      !issubjecti && continue # only check from subject polygon

      indi = evi.ind
      ai, bi = vertices(fillsegs[indi])

      for j in 1:length(startevs)
        evj = startevs[j]
        issubjectj = evj.issubject
        issubjecti == issubjectj && continue # skip subject polygon

        indj = evj.ind
        aj, bj = vertices(fillsegs[indj])

        sharestart = ai == aj
        samedir = sharestart && bi == bj
        revdir = ai == bj && bi == aj

        # if segments are collinear and sharing start point, swap fill info
        if samedir || revdir || (sharestart && sideof(bj, Line(ai, bi)) == ON)
          _swapfills!(fillsegs, indi, indj)
          break
        end
      end
    end

    BinaryTrees.delete!(events, p)
  end
end

# build event queue from segments
function _buildevents(fillsegs::FillSegments{S}, nsegsfirst::Int) where {S}
  P = eltype(vertices(first(fillsegs)))
  events = BinaryTrees.AVLTree{P,Vector{SegmentEvent}}()
  isempty(fillsegs) && return events

  for i in 1:length(fillsegs)
    issubject = i <= nsegsfirst
    s = fillsegs[i]
    a, b = vertices(s)
    a > b && ((a, b) = (b, a))
    _addevent!(events, a, SegmentEvent(i, true, issubject))
    _addevent!(events, b, SegmentEvent(i, false, issubject))
  end
  events
end

function _addevent!(events, p, ev)
  node = BinaryTrees.search(events, p)
  isnothing(node) ? BinaryTrees.insert!(events, p, [ev]) : push!(BinaryTrees.value(node), ev)
end

# swap fill information between two equal segments
function _swapfills!(fillsegs, i, j)
  # segment i is from subject polygon
  # segment j is from clip polygon
  # they overlap, so they should share fill information

  # get what each segment knows about its own polygon
  sai, sbi = _getsubjectfills(fillsegs, i)
  caj, cbj = _getclipfills(fillsegs, j)

  # segment i needs to know about clip polygon (copy from j)
  _setclipfills!(fillsegs, caj, cbj, i)
  # segment j needs to know about subject polygon (copy from i)
  _setsubjectfills!(fillsegs, sai, sbi, j)

  _setcomputed!(fillsegs, i)
  _setcomputed!(fillsegs, j)
end
