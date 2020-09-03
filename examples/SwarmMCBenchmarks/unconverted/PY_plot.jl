
using SwarmMC, PyPlot

figure()
ax1 = axes()
figure()
ax2 = axes()
figure()
ax3 = axes()

colors = Dict(0.0 => "k", 0.2 => "b", 0.3 => "r", 0.4 => "g")

for phi in [0.0, 0.2, 0.3, 0.4]
	restrictstr = Regex("PercusYevick_Phi$phi")
	@show restrictstr
	data = SwarmMC.GetAll(:ETd, :eps, :bulk_W, :bulk_DL, :bulk_DT, :n0, restrict=restrictstr, n0_factors=true, sortby=:ETd)
	data2 = convert(Array{Float64,2},data[2:end,:])
	ax1[:plot](data2[1,:], data2[3,:], color=colors[phi], label="Phi=$phi")
	ax2[:plot](data2[1,:], data2[4,:], color=colors[phi], label="Phi=$phi")
	ax2[:plot](data2[1,:], data2[5,:], color=colors[phi], linestyle="--")
	ax3[:plot](data2[1,:], data2[2,:], color=colors[phi], linestyle="--", label="Phi=$phi")
end

A = readcsv("W.csv", skipstart=1)
ax1[:plot](A[:,1], A[:,2], "o", label="Wade digitized")

ax1[:xscale]("log")
ax2[:xscale]("log")
ax3[:xscale]("log")
ax1[:yscale]("log")
ax2[:yscale]("log")
ax3[:yscale]("log")

ax1[:legend]()
ax2[:legend]()
ax3[:legend]()
