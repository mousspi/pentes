using CSV
using DataFrames
using DelaunayTriangulation
using DelaunayTriangulation: triangulate
using GeometryBasics
using LinearAlgebra
using CairoMakie

# === CHARGEMENT DU FICHIER CSV ===
filepath = "Exemple Mesure.csv"
df = CSV.read(filepath, DataFrame; delim=';', header=false)
rename!(df, [:x, :y, :z])

# === CONVERSION EN POINTS 2D ===
points2d = Point2f0.(df.x, df.y)
zvals = df.z

# === TRIANGULATION DELAUNAY ===
tri = triangulate(points2d)

# === FONCTION POUR CALCULER LA PENTE DANS UN TRIANGLE ===
function triangle_gradient(triangle, tri, z)
    # Vérifiez que les indices sont valides
    inds = [triangle[1], triangle[2], triangle[3]]
    if any(ind -> ind < 1 || ind > length(tri.points), inds)
        error("Indice invalide dans le triangle: $inds")
    end

    p1, p2, p3 = tri.points[inds[1]], tri.points[inds[2]], tri.points[inds[3]]
    z1, z2, z3 = z[inds[1]], z[inds[2]], z[inds[3]]

    # Vecteurs dans le plan
    v1 = Point2f0(p2 .- p1)
    v2 = Point2f0(p3 .- p1)
    A = hcat(Vec(v1), Vec(v2))
    dz = [z2 - z1, z3 - z1]

    # Gradient = dérivée partielle par rapport à x et y
    grad = A \ dz
    return norm(grad)
end

# === CALCUL DES PENTES MAXIMALES PAR TRIANGLE ===
pentes = [triangle_gradient(t, tri, zvals) for t in tri.triangles]

# === AFFICHAGE DE LA CARTE ===
f = Figure()
ax = Axis(f[1, 1], title="Carte des pentes locales (triangulation)", xlabel="x", ylabel="y")

# Créer une collection de triangles colorés
poly = [Triangle(tri.points[t[1]], tri.points[t[2]], tri.points[t[3]]) for t in tri.triangles]
polycollection!(ax, poly, color=pentes, colormap=:viridis)

# Barre de couleur
Colorbar(f[1, 2], limits=extrema(pentes), colormap=:viridis, label="Pente maximale")

# Sauvegarde optionnelle
# save("carte_pentes.png", f)

f
