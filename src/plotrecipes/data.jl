# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(data::Data, vars=nothing)
  # valid variables
  validvars = Tables.schema(values(data)).names

  # plot all variables by default
  isnothing(vars) && (vars = validvars)
  @assert vars ⊆ validvars "invalid variable names"

  # underlying domain where variables live
  dom = domain(data)

  # shared plot specs
  layout --> length(vars)

  for (i, var) in enumerate(vars)
    vals = data[var]

    # handle categorical values
    l = vals isa CategoricalArray ? levelcode.(vals) : copy(vals)

    # handle missing values
    replace!(l, missing => NaN)

    @series begin
      subplot := i
      title --> string(var)
      legend --> false
      dom, vals
    end
  end
end
