module npps
# npp_plot_section.jl
# DG 2019-12-03
#
using PyCall
using NaNMath
using FileIO, JLD2
#using Debugger

if (pwd() in LOAD_PATH) == false
    push!(LOAD_PATH, pwd());
end

import C2PO: findNaNmax, findNaNmin, meshgrid
import GridFit: gridfit
import NPP_types: CTD, Expedition, Transect

np = pyimport("numpy")
interpolate = pyimport("scipy.interpolate")
nm = NaNMath;

try !isassigned(ctd)
catch e
    import npp_load_ctd: ctd, npp
end

export transect, gridCTD, ctd, npp

function gridCTD(ctd::Array{CTD,1}, npp::Expedition, decim::Int64=8)
    transect = Transect[];
    for si = 1:length(npp.section)
    #for si = 8:8
        println(string(si, npp.sectname[si]))
        # loading the CTD cast numbers for each section
        sect = vec(npp.section[si]);
        x = []; y = []; ctmp = []; saltA = []; sigma0 = [];

        # finding boundaries
        xgridlim = Array{Float64,1}(undef,2);
        ygridlim = Array{Float64,1}(undef,2);
        bott = Array{Float64,1}(undef,length(sect));

        if length(sect) > 1
            xgridlim[1] = nm.min(xgridlim[1],ceil(findNaNmin(npp.s[si])[1]));
            xgridlim[2] = nm.max(xgridlim[2],floor(findNaNmax(npp.s[si])[1]));
            for ii = 1:length(sect)
                ygridlim[1] = nm.min(ygridlim[1],floor(findNaNmin(ctd[sect[ii]].z)[1]));
                ygridlim[2] = nm.max(ygridlim[2],min(0.0,ceil(findNaNmax(ctd[sect[ii]].z)[1])));
                #println(size(x))
                #append!(x, repeat([npp.s[si][ii]], [length(ctd[sect[ii]].p)]...)[1:decim:end]);
                gd = findall(ctd[sect[ii]].dataflag .> 0);
                x = vcat(x, repeat([npp.s[si][ii]], [length(ctd[sect[ii]].p[gd])]...)[1:decim:end]);
                y = vcat(y, ctd[sect[ii]].z[gd][1:decim:end]);
                ctmp = vcat(ctmp, ctd[sect[ii]].ctmp[gd][1:decim:end]);
                saltA = vcat(saltA, ctd[sect[ii]].saltA[gd][1:decim:end]);
                sigma0 = vcat(sigma0, ctd[sect[ii]].sigma0[gd][1:decim:end]);

                bott[ii] = nm.minimum(npps.ctd[sect[ii]].z);
                println(length(x))
            end
            xnodes = collect(xgridlim[1]:1.0:xgridlim[2]);
            ynodes = collect(ygridlim[1]:1.0:ygridlim[2]);
            (xgrid,ygrid) = meshgrid(xnodes,ynodes);

            # zgrid = gridfit(x,y,z,xgrid,ygrid); using gridfit to calculate the gridded version of the data. this is where the magic happens!!!
            ctmpgrid = gridfit(Float64.(x), Float64.(y), Float64.(ctmp), Float64.(xnodes), Float64.(ynodes), smoothness=[1.0 1.0], regularizer="gradient");
            saltAgrid = gridfit(Float64.(x), Float64.(y), Float64.(saltA), Float64.(xnodes), Float64.(ynodes), smoothness=[1.0 1.0], regularizer="gradient");
            sigma0grid = gridfit(Float64.(x), Float64.(y), Float64.(sigma0), Float64.(xnodes), Float64.(ynodes), smoothness=[1.0 1.0], regularizer="gradient");

            #println(bott)

            #find bottom depth for each grid location so that the portion below the bottom can be masked.
            if length(sect) > 2
                bottf = interpolate.interp1d(npp.s[si], bott, kind="cubic");
            elseif length(sect) == 2
                bottf = interpolate.interp1d(npp.s[si], bott, kind="linear");
            end
            botti = bottf(xnodes);
            #println(botti)

            # set portion of the gridded below the bottom depth to NaN
            for ii = 1:length(xnodes)
                subbi = findall(ygrid[:,ii] .< botti[ii]); # indices for sub-bottom
                if !isempty(subbi)
                    ctmpgrid[subbi,ii] .= NaN;
                    saltAgrid[subbi,ii] .= NaN;
                    sigma0grid[subbi,ii] .= NaN;
                end
            end

            transect = push!(transect, Transect(si, npp.sectname[si], Float64.(x), Float64.(y), Float64.(ctmp), Float64.(saltA), Float64.(sigma0), xgrid, ygrid, ctmpgrid, saltAgrid, sigma0grid, vec(npp.s[si]), bott))
        else
            transect = push!(transect, Transect(si, npp.sectname[si], [NaN], [NaN], [NaN], [NaN], [NaN], [NaN NaN], [NaN NaN], [NaN NaN], [NaN NaN], [NaN NaN], [NaN], [NaN]));
        end #if
    end #for

    #return (x,y,ctmp,saltA,xnodes,ynodes,xgrid,ygrid,ctmpgrid,saltAgrid)
    return transect
end

reloadflag = 1
loadsavedflag = 0
saveflag = 0

# this is the call to the gridCTD function that does the analysis
if reloadflag == 1
    transect = gridCTD(ctd, npp);
end

# not working for some reason
if loadsavedflag == 1 & reloadflag != 1
    workdir = "/Users/gong/Research/NPP/data/CTD/";
    filepath = string(workdir, "NPP_ctd_data_gridded", ".jld2")

    @time project = load(filepath)
    transect = project["transect"];
    npp = project["npp"]
end

# not working right
if saveflag == 1
    workdir = "/Users/gong/Research/NPP/data/CTD/";
    filepath = string(workdir, "NPP_ctd_data_gridded", ".jld2")
    @time @save filepath transect npp
end

end # module
