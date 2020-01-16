# npp_plot_TS.jl
using PyCall, PyPlot, NaNMath
cmo = pyimport("cmocean")
plt = PyPlot;
nm = NaNMath;
import Gong

if (@isdefined ctd) == false
    include("run_npp_load_ctd.jl")
end

function load_TS(ctd)
    gind = [];
    ctdtemp = Array{Float64}[];
    ctdsalt = Array{Float64}[];
    ctdsigma0 = Array{Float64}[];
    ctdoxy = Array{Float64}[];
    ctdlon = Array{Float64}[];
    ctdlat = Array{Float64}[];

    for i = 1:length(ctd)
        gind = push!(gind,findall(ctd[i].dataflag .> 0));
        #ctd[i].temp[gind];
    end

    for i = 1:length(ctd)
        ctdtemp = cat(ctdtemp, ctd[i].temp[gind[i]], dims = 1);
        ctdsalt = cat(ctdsalt, ctd[i].salt[gind[i]], dims = 1);
        ctdsigma0 = cat(ctdsigma0, ctd[i].sigma0[gind[i]], dims = 1);
        ctdoxy = cat(ctdoxy, ctd[i].oxy[gind[i]], dims = 1);
        ctdlon = cat(ctdlon, ctd[i].lon[gind[i]], dims = 1);
        ctdlat = cat(ctdlat, ctd[i].lat[gind[i]], dims = 1);
    end

    return Float64.(ctdtemp), Float64.(ctdsalt), Float64.(ctdsigma0), Float64.(ctdoxy), Float64.(ctdlon), Float64.(ctdlat)
end

(temp, salt, sigma0, oxy, lon, lat) = load_TS(ctd);

(tempscl, tempmin, tempmax) = Gong.minmaxscaler(temp,[0 1]);
(saltscl, saltmin, saltmax) = Gong.minmaxscaler(salt,[0 1]);
(oxyscl, oxymin, oxymax) = Gong.minmaxscaler(oxy,[0 1]);

datascl = permutedims([tempscl saltscl oxyscl]);



figTS = plt.figure(1)
clf()
ax = figTS.add_subplot(1, 1, 1)
#axT = fig.add_subplot(2, 2, 1)
hpTS = plt.scatter3D(salt, temp, oxy, s=1.0, c=sigma0, cmap=cmo.cm.thermal) #mpl.cm.coolwarm
#hcT = ax.contour(transect[i].xgrid, transect[i].ygrid, transect[i].sigma0grid, colors="k")
#ax.clabel(hcT, inline=1, fontsize=9);
ax.set_title("NPP T/S Diagram")
ax.set_xlabel("Salinity")
ax.set_ylabel("Temperature")
#ax.set_zlabel("Oxygen")
#ax.set_xlim(xmin,xmax)
#ax.set_ylim(ymin,ymax)
#hpTS.set_clim(-2,2)
figTS.colorbar(hpTS,ax=ax)
#fig.colorbar(hpT1,ax=axT)
#plt.savefig(string(figoutdir, "ctmp_", i, npp.sectname[i], "_full.png"), dpi = 300)
