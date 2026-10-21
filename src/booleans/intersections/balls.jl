# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

intersection(f, ball₁::Ball{𝔼{2}}, ball₂::Ball{𝔼{2}}) = _polygonintersect(f, ball₁, ball₂)

intersection(f, ball₁::Ball{🌐}, ball₂::Ball{🌐}) = _polygonintersect(f, ball₁, ball₂)

intersection(f, ball::Ball, poly::Polygon) = _polygonintersect(f, ball, poly)
