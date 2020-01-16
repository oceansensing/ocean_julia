# this scripts plot the hydrographic sections for NPP
using PyCall, PyPlot, NaNMath
cmo = pyimport("cmocean")
plt = PyPlot;
nm = NaNMath;

include("npp_path_setup.jl") # for donglai's computer, please use file specific to your computer

try !isassigned(transect)
catch e
    import npps: transect, npp
end

for i = 1:length(transect)
#for i = 9:9

    # only plot if there is a section to plot, meaning more than one station.
    if length(npp.section[i]) > 1
        # find the xlim range for plotting

        if i != 9
            xrd = round.(ceil.(vec(npp.s[i])) ./ 10.0) .* 10.0; # along track distance x rounded to the 10 km
            xrd = append!(xrd, round.(floor.(vec(npp.s[i])) ./ 10.0) .* 10.0);
            xmin = nm.minimum(xrd);
            xmax = nm.maximum(xrd);
        else
            xmin = floor(minimum(vec(npp.s[i])));
            xmax = ceil(maximum(npp.s[i]));
        end

        # find the ylim range for plotting
        ymin = nm.minimum(unique(transect[i].ygrid))-10.0;
        ymax = nm.maximum(unique(transect[i].ygrid));

        # plotting conservative temperature, full water column
        #fig = plt.figure(1)
        figT = plt.figure(1)
        clf()
        axT = figT.add_subplot(1, 1, 1)
        #axT = fig.add_subplot(2, 2, 1)

        hpT1 = axT.pcolor(transect[i].xgrid, transect[i].ygrid, transect[i].ctmpgrid, cmap=cmo.cm.thermal)
        hpT2 = axT.scatter(transect[i].x, transect[i].y, s=1.0, c=transect[i].ctmp, cmap=cmo.cm.thermal) #mpl.cm.coolwarm
        hcT = axT.contour(transect[i].xgrid, transect[i].ygrid, transect[i].sigma0grid, colors="k")
        axT.clabel(hcT, inline=1, fontsize=9);
        axT.scatter(transect[i].s, repeat([0.5], length(transect[i].s)), c="black", marker="|")
        axT.set_title(string(npp.sectname[i], " (CT)"))
        axT.set_xlabel("Along track distance (km)")
        axT.set_ylabel("Depth (m)")
        axT.set_xlim(xmin,xmax)
        axT.set_ylim(ymin,ymax)
        hpT1.set_clim(-2,2)
        hpT2.set_clim(-2,2)
        figT.colorbar(hpT1,ax=axT)
        #fig.colorbar(hpT1,ax=axT)
        plt.savefig(string(figoutdir, "ctmp_", i, npp.sectname[i], "_full.png"), dpi = 300)

        # plotting conservative temperature, upper 100 m
        figTs = plt.figure(2)
        clf()
        axTs = figTs.add_subplot(1, 1, 1)
        #axTs = figT.add_subplot(2, 2, 2)

        hpTs1 = axTs.pcolor(transect[i].xgrid, transect[i].ygrid, transect[i].ctmpgrid, cmap=cmo.cm.thermal)
        hpTs2 = axTs.scatter(transect[i].x, transect[i].y, s=1.0, c=transect[i].ctmp, cmap=cmo.cm.thermal) #mpl.cm.coolwarm
        hcTs = axTs.contour(transect[i].xgrid, transect[i].ygrid, transect[i].sigma0grid, colors="k")
        axTs.clabel(hcTs, inline=1, fontsize=9);
        axTs.scatter(transect[i].s, repeat([0.5], length(transect[i].s)), c="black", marker="|")
        axTs.set_title(string(npp.sectname[i], " (CT)"))
        axTs.set_xlabel("Along track distance (km)")
        axTs.set_ylabel("Depth (m)")
        axTs.set_xlim(xmin,xmax)
        axTs.set_ylim(-100,0)
        hpTs1.set_clim(-2,2)
        hpTs2.set_clim(-2,2)
        figTs.colorbar(hpTs1,ax=axTs)
        #fig.colorbar(hpTs1,ax=axTs)
        plt.savefig(string(figoutdir, "ctmp_", i, npp.sectname[i], "_100.png"), dpi = 300)

        # plotting absolute salinity, full water column
        figS = plt.figure(3)
        clf()
        axS = figS.add_subplot(1, 1, 1)
        #axS = fig.add_subplot(2, 2, 3)

        hpS1 = axS.pcolor(transect[i].xgrid, transect[i].ygrid, transect[i].saltAgrid, cmap=cmo.cm.haline)
        hpS2 = axS.scatter(transect[i].x, transect[i].y, s=1.0, c=transect[i].saltA, cmap=cmo.cm.haline) #mpl.cm.coolwarm
        hcS = axS.contour(transect[i].xgrid, transect[i].ygrid, transect[i].sigma0grid, colors="k")
        axS.clabel(hcS, inline=1, fontsize=9);
        axS.scatter(transect[i].s, repeat([0.5], length(transect[i].s)), c="black", marker="|")
        axS.set_title(string(npp.sectname[i], " (SA)"))
        axS.set_xlabel("Along track distance (km)")
        axS.set_ylabel("Depth (m)")
        axS.set_xlim(xmin,xmax)
        axS.set_ylim(ymin,ymax)
        hpS1.set_clim(28,35)
        hpS2.set_clim(28,35)
        figS.colorbar(hpS1,ax=axS)
        #figS.colorbar(hpS1,ax=axS)
        plt.savefig(string(figoutdir, "saltA_", i, npp.sectname[i], "_full.png"), dpi = 300)

        # plotting absolute salinity, upper 100m
        figSs = plt.figure(4)
        clf()
        axSs = figSs.add_subplot(1, 1, 1)
        #axSs = fig.add_subplot(2, 2, 4)

        hpSs1 = axSs.pcolor(transect[i].xgrid, transect[i].ygrid, transect[i].saltAgrid, cmap=cmo.cm.haline)
        hpSs2 = axSs.scatter(transect[i].x, transect[i].y, s=1.0, c=transect[i].saltA, cmap=cmo.cm.haline) #mpl.cm.coolwarm
        hcSs = axSs.contour(transect[i].xgrid, transect[i].ygrid, transect[i].sigma0grid, colors="k")
        axSs.clabel(hcSs, inline=1, fontsize=9);
        axSs.scatter(transect[i].s, repeat([0.5], length(transect[i].s)), c="black", marker="|")
        axSs.set_title(string(npp.sectname[i], " (SA)"))
        axSs.set_xlabel("Along track distance (km)")
        axSs.set_ylabel("Depth (m)")
        axSs.set_xlim(xmin,xmax)
        axSs.set_ylim(-100,0)
        hpSs1.set_clim(28,35)
        hpSs2.set_clim(28,35)
        figSs.colorbar(hpSs1,ax=axSs)
        #fig.colorbar(hpSs1,ax=axSs)
        plt.savefig(string(figoutdir, "saltA_", i, npp.sectname[i], "_100.png"), dpi = 300)
    end #if

    close()
end #for
