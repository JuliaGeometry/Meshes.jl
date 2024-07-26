# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Manifold

A manifold where geometries and domains are defined.
"""
abstract type Manifold end

"""
    𝔼{Dim}

Euclidean manifold with dimension `Dim`.
"""
abstract type 𝔼{Dim} <: Manifold end

"""
    🌐

Ellipsoid manifold for geodesic geometry.
"""
abstract type 🌐 <: Manifold end
