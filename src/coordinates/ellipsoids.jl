# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const wgs84 = let 
  a = 6378137.0
  f⁻¹ = 298.257223563
  f = inv(f⁻¹)
  b = a * (1 - f)
  e² = (2 - f) / f⁻¹
  e = √e²
  (; a = a * u"m", b = b * u"m", e, e², f, f⁻¹)
end
