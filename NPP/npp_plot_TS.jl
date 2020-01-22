# npp_plot_TS.jl
using PyCall, PyPlot, NaNMath, Clustering
plt = PyPlot;
nm = NaNMath;
cmo = pyimport("cmocean")
#cluster = pyimport("sklearn.cluster")
import Gong

calcflag = 0
plotflag = 2

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

if calcflag == 1
    (temp, salt, sigma0, oxy, lon, lat) = load_TS(ctd);

    (tempscl, tempmin, tempmax) = Gong.minmaxscaler(temp,[0 1]);
    (saltscl, saltmin, saltmax) = Gong.minmaxscaler(salt,[0 1]);
    (oxyscl, oxymin, oxymax) = Gong.minmaxscaler(oxy,[0 1]);

    datascl = permutedims([tempscl saltscl oxyscl]);
    clusters = dbscan(datascl[:,1:4:end], 0.05, min_neighbors = 100, min_cluster_size = 4000);
end

if plotflag == 2
    figCluster = plt.figure(2)
    clf()
    axC = figCluster.add_subplot(1, 1, 1)
    #for i = 1:length(clusters)
    for i = 1:1
        println(i)
        corei = clusters[i].core_indices;
        bndyi = clusters[i].boundary_indices;
        hpC = plt.scatter(salt[corei],temp[corei],s=1.0)
        hpB = plt.scatter(salt[bndyi],temp[bndyi],s=1.0)
    end
    axC.set_title("NPP T/S Diagram")
    axC.set_xlabel("Salinity")
    axC.set_ylabel("Temperature")
end


#datascl = [tempscl saltscl oxyscl];
#db = cluster.dbscan(datascl[1:4:end,:], eps=0.1, min_samples=10)

if plotflag == 1
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
end
