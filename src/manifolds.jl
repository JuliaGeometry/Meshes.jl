# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    𝔼{Dim}

Euclidean manifold with dimension `Dim`.
"""
abstract type 𝔼{Dim} <: AbstractManifold{ℝ} end

manifold_dimension(::Type{𝔼{Dim}}) where {Dim} = Dim

"""
    🌐

Ellipsoid manifold for geodesic geometry.
"""
abstract type 🌐 <: AbstractManifold{ℝ} end

manifold_dimension(::Type{🌐}) = 2
