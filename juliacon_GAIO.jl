const a, b = 1.35, 0.3

f( (x,y) ) = ( 1 - a*x^2 + y,  b*x )

orbit = NTuple{2,Float64}[ (1.0, 1.0) ]
for k = 1:10_000
    v = orbit[end]
    push!( orbit, f(v) )
end

using Plots
scatter(
    orbit, 
    marker_z = 1:length(orbit),
    color=:haline,
    alpha=0.4,
    lab=false
)

using GAIO

center, radius = (0,0), (3,3)
domain = Box(center, radius)

P = BoxPartition(domain, (2,2))
F = BoxMap(:montecarlo, f, domain)

function relative_attractor(F::BoxMap, B::BoxSet, steps)
    for k = 1:steps
        B = subdivide(B)
        B = B ∩ F(B)
    end
    return B
end

B = cover(P, :)
B = relative_attractor(F, B, 6)

plot(B)

F♯ = TransferOperator(F, B, B)
λs, μs, _ = eigs(F♯)
μ = log ∘ abs ∘ μs[1]

plot(μ, lab=false)
lens!(
    [0.1, 0.5], [0.15, 0.25], 
    inset=(1, bbox(0.1,0.2,0.4,0.5))
)