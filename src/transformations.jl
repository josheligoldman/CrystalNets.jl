"""
    equiv_mapping(pg1, pg2) -> Union{Nothing, PeriodicGraphTransformation{3,9}}

Return the `PeriodicGraphTransformation` mapping `pg1` to `pg2` via their shared canonical
form, or `nothing` if the two graphs are not topologically equivalent or if CrystalNets
could not compute the canonicalizing transformation. Computed as `inv(T2) * T1` where `T1`
(resp. `T2`) is the transformation CrystalNets applies to canonicalize `pg1` (resp. `pg2`).
"""
function equiv_mapping(pg1::PeriodicGraph{3}, pg2::PeriodicGraph{3})::Union{Nothing, PeriodicGraphTransformation{3,9}}
    nv(pg1) != nv(pg2) && return nothing
    ne(pg1) != ne(pg2) && return nothing
    opts = CrystalNets.Options(
        skip_minimize=true,
        export_input=false,
        export_trimmed=false,
        export_attributions=false,
        export_clusters=false,
        export_net=false,
        export_subnets=false,
    )
    topo1 = last(only(only(topological_genome(pg1, opts))))
    topo2 = last(only(only(topological_genome(pg2, opts))))
    string(topo1.genome) != string(topo2.genome) && return nothing
    (isnothing(topo1.transformation) || isnothing(topo2.transformation)) && return nothing
    return inv(topo2.transformation) * topo1.transformation
end

function isequiv(pg1::PeriodicGraph{3}, pg2::PeriodicGraph{3})::Bool
    return !isnothing(equiv_mapping(pg1, pg2))
end