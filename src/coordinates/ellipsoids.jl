# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# WGS84 ellipsoid
const a = 6378137.0
const f⁻¹ = 298.257223563
const f = inv(f⁻¹)
const b = a * (1 - f)
const e² = (2 - f) / f⁻¹
const e = √e²
