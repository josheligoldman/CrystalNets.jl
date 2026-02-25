function (pgt::PeriodicGraphTransformation3D)(crystal::CrystalNets.Crystal{Nothing})
    transformed_pge = pgt(crystal.pge)
    transformed_types = crystal.types[pgt.vertex_permutation]
    return CrystalNets.Crystal(transformed_pge, transformed_types, crystal.options)
end
