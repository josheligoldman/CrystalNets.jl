function (pgt::PeriodicGraphTransformation3D)(crystal::Crystal{Nothing})
    transformed_pge = pgt(crystal.pge)
    transformed_types = crystal.types[pgt.vertex_permutation]
    return Crystal{Nothing}(transformed_pge, transformed_types, crystal.options)
end
