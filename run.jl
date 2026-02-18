using StaticArrays
using CrystalNets

cif_path = "data/zeolites/8000347.cif"
structure = StructureType.Zeolite       # Zeolite, MOF
clustering = Clustering.EachVertex      # SingleNodes, AllNodes, EachVertex

track_mapping = Int[]
basis_mapping = zero(MMatrix{3, 3, Int, 9})
options = CrystalNets.Options(
    structure=structure,
    clusterings=[clustering],
    detect_organiccycles=false,
    detect_pe=false,
    export_input=false,
    export_trimmed=false,
    export_attributions=false,
    export_clusters=false,
    export_net=false,
    export_subnets=false,
    skip_minimize=true,
    track_mapping=track_mapping,
    basis_mapping=basis_mapping,
)

topology_result = determine_topology(cif_path, options)
@show basis_mapping
