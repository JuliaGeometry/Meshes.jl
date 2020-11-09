# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
Ordering convention used to construct k-faces of a polytope from vertices (0-faces).
"""
abstract type Ordering end

"""
Connectivity function provided by an `Ordering` subtype used to construct k-faces of a polytope from vertices (0-faces).
"""
function connectivity end