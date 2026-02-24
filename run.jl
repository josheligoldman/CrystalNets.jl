using PeriodicGraphs
using CrystalNets

x = PeriodicGraph3D(2, PeriodicEdge3D[
    (1, 2, (-1, 0, -1)), (1, 2, (-1, 0, 0)),
    (1, 2, (0, 0, -1)),  (1, 2, (0, 0, 0)),
    (1, 2, (0, 1, -1)),  (1, 2, (0, 1, 0)),
])

result = topological_genome(
    x;
    skip_minimize=true,
    export_input=false,
    export_trimmed=false,
    export_attributions=false,
    export_clusters=false,
    export_net=false,
    export_subnets=false,
)
result = first(only(result).results)
@show result

pgt = result.transformation
genome = result.genome

@show pgt

@show apply_transform(pgt, x)
@show genome

@show genome == apply_transform(pgt, x)
@show genome == x

@assert apply_transform(pgt, x) == genome
