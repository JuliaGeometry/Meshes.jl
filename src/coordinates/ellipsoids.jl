# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const wgs84 = let 
  a = 6378137.0 * u"m"
  f⁻¹ = 298.257223563
  f = inv(f⁻¹)
  b = a * (1 - f)
  e² = (2 - f) / f⁻¹
  e = √e²
  (; a , b , e, e², f, f⁻¹)
end
