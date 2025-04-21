using CSV, DataFrames, DelaunayTriangulation, LinearAlgebra, CairoMakie

# 1. Chargement du nuage initial
df = CSV.read("Pts CSV contact 13216.csv", DataFrame; delim=';')
X = [Tuple(row) for row in eachrow(df)]  # (x, y, z)
xy = [(x, y) for (x, y, z) in X]
tri = DelaunayTriangulation.delaunay(xy)

# 2. Calcul des normales
function normal(tri_pts)
    p1, p2, p3 = tri_pts
    v1, v2 = p2 .- p1, p3 .- p1
    n = cross(v1, v2)
    return norm(n) == 0 ? n : n / norm(n)
end
triangle_normals = [normal([X[i] for i in t]) for t in tri]

# 3. Construction du dictionnaire arête → triangles
edge_map = Dict{Tuple{Int, Int}, Vector{Int}}()
for (i, t) in enumerate(tri)
    for e in [(t[1], t[2]), (t[2], t[3]), (t[3], t[1])]
        push!(get!(edge_map, Tuple(sort(e)), Int[]), i)
    end
end

# 4. Création du nuage de points dérivé
edge_points = Vector{NamedTuple{(:x, :y, :z, :angle_mrad), Tuple{Float64, Float64, Float64, Float64}}}()

for ((i, j), tris) in edge_map
    if length(tris) == 2
        n1, n2 = triangle_normals[tris[1]], triangle_normals[tris[2]]
        angle_rad = acos(clamp(dot(n1, n2), -1.0, 1.0))
        angle_mrad = angle_rad * 1000

        p1, p2 = X[i], X[j]
        mid = 0.5 .* (p1 .+ p2)
        push!(edge_points, (x=mid[1], y=mid[2], z=mid[3], angle_mrad=angle_mrad))
    end
end

# 5. Sauvegarde dans un fichier CSV
df_edges = DataFrame(edge_points)
#CSV.write("arêtes_avec_pentes.csv", df_edges)

# 6. Visualisation 3D avec couleurs
x = df_edges.x
y = df_edges.y
z = df_edges.z
a = df_edges.angle_mrad

f = Figure()
ax = Axis3(f[1, 1], xlabel="x", ylabel="y", zlabel="z", title="Ruptures locales (angle entre triangles)")
scatter!(ax, x, y, z; color=a, markersize=10, colormap=:inferno)
Colorbar(f[1, 2], ax.scene.plots[1], label="Angle (mrad)")
f
