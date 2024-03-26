# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    isclosed(chain)

Tells whether or not the `chain` is closed.

A closed `chain` is also known as a ring.
"""
isclosed(c::Chain) = isclosed(typeof(c))

isclosed(::Type{<:Segment}) = false

isclosed(::Type{<:Rope}) = false

isclosed(::Type{<:Ring}) = true
