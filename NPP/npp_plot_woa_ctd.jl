using PyPlot
plt = PyPlot;

# if NPP CTD data is not yet loaded, load it.
if (@isdefined ctd) == false
    import npp_load_ctd: ctd, npp
end

# if WOA data is not yet loaded, load it.
i = 3; #1: winter, 2: spring, 3: summer, 4: fall
if (@isdefined load_woa) == false
    import WOA: npp_load_woa, load_woa, season
end
season[i]
profile = npp_load_woa(i);

# load profile plotting flag
include("npp_flags.jl")
#woaplotflag = 1 # change this to 1 if you want to make NPP CTD & WOA comparison plots, not active, use npp_flags.jl

if woaplotflag == 1
        # make plots of comparison profiles between NPP CTD and closest WOA location
    for j = 1:length(ctd)
        println(string("Plotting Station ", j))
        gind = findall(ctd[j].dataflag .> 0);

        figT = plt.figure(1)
        clf()
        axT = figT.add_subplot(1, 1, 1)

        axT.plot(profile[j].temp, profile[j].z, c="black")
        axT.plot(ctd[j].temp[gind], ctd[j].z[gind], c="red")
        axT.set_title(string("Station ", j, " (Temperature)"))
        axT.set_xlabel("Temperature (C)")
        axT.set_ylabel("Depth (m)")
        axT.legend(["WOA", "NPP"])
        plt.savefig(string(figoutdir, "temp_sta_", string(j/10000)[end-1:end], "_", season[i], ".png"), dpi = 300)

        figS = plt.figure(2)
        clf()
        axS = figS.add_subplot(1, 1, 1)

        axS.plot(profile[j].salt, profile[j].z, c="black")
        axS.plot(ctd[j].salt[gind], ctd[j].z[gind], c="red")
        axS.set_title(string("Station ", j, " (Salinity)"))
        axS.set_xlabel("Salinity")
        axS.set_ylabel("Depth (m)")
        axS.legend(["WOA", "NPP"])
        plt.savefig(string(figoutdir, "salt_sta_", string(j/10000)[end-1:end], "_", season[i], ".png"), dpi = 300)
    end
    close()

end
