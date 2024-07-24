# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

abstract type AbstractManifold end

"""
    𝔼{Dim}

Euclidean manifold with dimension `Dim`.
"""
abstract type 𝔼{Dim} <: AbstractManifold end

"""
    🌐

Ellipsoid manifold for geodesic geometry.
"""
abstract type 🌐 <: AbstractManifold end
