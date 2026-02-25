function apply_transform(pgt::PeriodicGraphTransformation3D, crystal::CrystalNets.Crystal{Nothing})
    transformed_pge = apply_transform(pgt, crystal.pge)
    transformed_types = crystal.types[pgt.vertex_permutation]
    return CrystalNets.Crystal(transformed_pge, transformed_types, crystal.options)
end

(pgt::PeriodicGraphTransformation3D)(crystal::CrystalNets.Crystal{Nothing}) = apply_transform(pgt, crystal)
