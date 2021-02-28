# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(partition::Partition)
  marker --> :auto
  for object in partition
    @series begin
      object
    end
  end
end
