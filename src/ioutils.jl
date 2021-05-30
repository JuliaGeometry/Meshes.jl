# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# helper function to print a large iterator
# in multiple lines with a given tabulation
function io_lines(itr, tab="  ")
  vec = collect(itr)
  N = length(vec)
  I, J = N > 10 ? (5, N-4) : (N, N+1)
  lines = [["$(tab)└─$(vec[i])" for i in 1:I]
           (N > 10 ? ["$(tab)⋮"] : [])
           ["$(tab)└─$(vec[i])" for i in J:N]]
  join(lines, "\n")
end
