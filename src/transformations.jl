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
        compute_pgt=true,
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

"""
    isequiv(pge1::PeriodicGraphEmbedding{3}, pge2::PeriodicGraphEmbedding{3}) -> Bool

Return `true` if the two periodic graph embeddings are equivalent.

Two paths:
- **Fast path**: when `pge1.g == pge2.g` byte-for-byte, the only transformation needed
  is the identity. Compare cells (`isapprox` on `cell.mat`) and per-vertex fractional
  positions (`isapprox`) directly.
- **Slow path**: otherwise, run `equiv_mapping(pge1.g, pge2.g)` to find a transformation
  through the shared canonical form, apply it to `pge1`, and compare cells and positions.

The cell comparison uses `rtol=sqrt(eps(Float64))` rather than the default BigFloat
tolerance, because `cell.mat` is typed `BigFloat` but in practice often stores Float64-
quality data (e.g. constructed by widening an `SMatrix{3,3,Float64}`). The default
BigFloat `isapprox` tolerance (`sqrt(eps(BigFloat)) ≈ 1e-38`) would otherwise reject
cells that differ only by Float64 round-off.

Caveats (slow path only):
- No continuous translation is applied — embeddings related by a global fractional shift
  are not considered equivalent.
- The basis change in the returned `PeriodicGraphTransformation` comes from canonicalization,
  not from any user-supplied basis. When the underlying graph has nontrivial automorphisms,
  two embeddings related by a basis change may be rejected because the canonicalizing PGT
  does not select that particular basis. To avoid this, ensure the graphs are byte-for-byte
  equal so the fast path applies.
"""
function isequiv(pge1::PeriodicGraphEmbedding{3}, pge2::PeriodicGraphEmbedding{3})::Bool
    # Fast path: when the graphs are byte-for-byte equal, the only transformation
    # needed is the identity. The canonical-form route (`equiv_mapping`) might
    # otherwise return a non-identity automorphism that re-permutes vertices or
    # rotates the cell — even though the embeddings match — because automorphic
    # graphs admit multiple valid canonical labelings.
    if pge1.g == pge2.g
        isapprox(pge1.cell.mat, pge2.cell.mat; rtol=sqrt(eps(Float64))) || return false
        length(pge1.pos) == length(pge2.pos) || return false
        return all(isapprox(p1, p2) for (p1, p2) in zip(pge1.pos, pge2.pos))
    end
    pgt = equiv_mapping(pge1.g, pge2.g)
    isnothing(pgt) && return false
    pge1_ = pgt(pge1)
    isapprox(pge1_.cell.mat, pge2.cell.mat; rtol=sqrt(eps(Float64))) || return false
    length(pge1_.pos) == length(pge2.pos) || return false
    return all(isapprox(p1, p2) for (p1, p2) in zip(pge1_.pos, pge2.pos))
end