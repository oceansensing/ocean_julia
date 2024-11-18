using NCDatasets, Dates, GLMakie, Colors, ColorSchemes

include("/Users/gong/GitHub/ocean_julia/C2PO.jl")
import .C2PO: datetime2unix, unix2datetime, datenum2datetime, unix2yearday, yearday2datetime, missing2nan

# reload flags
hycom_rf = false
bathy_rf = false
marietharp_rf = false

# plot flags
plotflag = true
bathy_pf = true
hycom_pf = true
marietharp_pf = true

figoutdir1 = "/Users/gong/oceansensing Dropbox/C2PO/NESMA/marietharp_2024/marietharp_nesma/TSG/";
figoutdir2 = "/Users/gong/oceansensing Dropbox/Donglai Gong/Projects/NESMA-PASSENGERS/ONR_annual_report_TFO/2024/";

t0 = DateTime(2024, 7, 21, 0, 0, 0);
z0 = 10.0;

#latmin, latmax = 36, 42;
#lonmin, lonmax = -67, -55;   
latmin, latmax = 37.5, 40.5;
lonmin, lonmax = -66.5, -59.8;
lat0 = 39;

if hycom_rf == true

    lonmin_hycom, lonmax_hycom = (lonmin < 0 ? lonmin + 360 : lonmin), (lonmax < 0 ? lonmax + 360 : lonmin)
    latmin_hycom, latmax_hycom = latmin, latmax;

    lonmin_hycom, lonmax_hycom = lonmin_hycom - 0.1, lonmax_hycom + 0.1;
    latmin_hycom, latmax_hycom = latmin_hycom - 0.1, latmax_hycom + 0.1;

    hycomurl = "https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0";
    ds = NCDataset(hycomurl);

    t_hycom = ds["time"][:];
    lon_hycom = ds["lon"][:];
    lat_hycom = ds["lat"][:];
    z_hycom = ds["depth"][:];

    latind_hycom = findall(latmin_hycom .<= lat_hycom .<= latmax_hycom);
    lonind_hycom = findall(lonmin_hycom .<= lon_hycom .<= lonmax_hycom);

    #tind = findall(DateTime(2024, 8, 24) .>= t_hycom .>= DateTime(2024, 7, 10))
    tind_hycom = findall(t_hycom .== t0);
    zind_hycom = findall(z_hycom .== z0);
    unixt_hycom = datetime2unix.(t_hycom[tind_hycom]);

    temp_hycom = missing2nan(ds["water_temp"][lonind_hycom, latind_hycom, zind, tind]);
    salt_hycom = missing2nan(ds["salinity"][lonind_hycom, latind_hycom, zind, tind]);
    u_hycom = missing2nan(ds["water_u"][lonind_hycom, latind_hycom, zind, tind]);
    v_hycom = missing2nan(ds["water_v"][lonind_hycom, latind_hycom, zind, tind]);
    el_hycom = missing2nan(ds["surf_el"][lonind_hycom, latind_hycom, tind]);
end

if bathy_rf == true
    ## plot bathymetry 
    bathypath = "/Users/gong/oceansensing Dropbox/C2PO/Data/bathy/ETOPO1/ETOPO_2022_v1_30s_N90W180_surface.nc";
    bathyds = Dataset(bathypath,"r");

    lon_bathy = bathyds["lon"][:];
    lat_bathy = bathyds["lat"][:];

    # extract indices from the bathymetric data file
    latind = findall(latmin .<= lat_bathy .<= latmax);
    lonind = findall(lonmin .<= lon_bathy .<= lonmax);

    z_bathy = Float64.(bathyds["z"][lonind, latind]); # recasting as Float64 to fix a StackOverFlow error seen in GLMakie 0.6.0.
    x_bathy = lon_bathy[lonind];
    y_bathy = lat_bathy[latind];

    pzind = findall(z_bathy .> 1);
    nzind = findall(z_bathy .< -1);
    zzind = findall(-1 .<= z_bathy .<= 1);
end

if marietharp_rf == true
    figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/NESMA/marietharp_2024/marietharp_nesma/TSG/";

    nesma_leg1_path = "/Users/gong/oceansensing Dropbox/C2PO/NESMA/marietharp_2024/marietharp_nesma/TSG/20240720T150306_processed_MicroTSG - NESMA_MT_Leg1_2024.csv"
    nesma_leg2_path = "/Users/gong/oceansensing Dropbox/C2PO/NESMA/marietharp_2024/marietharp_nesma/TSG/20240908T120357_processed_MicroTSG_NESMA_MT_Leg2_2024.csv"
    
    nesma_leg1 = CSV.read(nesma_leg1_path, DataFrame);
    nesma_leg2 = CSV.read(nesma_leg2_path, DataFrame);
    
    t1 = DateTime.(nesma_leg1.Timestamp, "yyyymmddTHHMMSS");
    lon1 = nesma_leg1.Longitude_DD;
    lat1 = nesma_leg1.Latitude_DD;
    temp1 = nesma_leg1.Temperature;
    salt1 = nesma_leg1.Salinity;
    cond1 = nesma_leg1.Conductivity;
    
    t2 = DateTime.(nesma_leg2.Timestamp, "yyyymmddTHHMMSS");
    lon2 = nesma_leg2.Longitude_DD;
    lat2 = nesma_leg2.Latitude_DD;
    temp2 = nesma_leg2.Temperature;
    salt2 = nesma_leg2.Salinity;
    cond2 = nesma_leg2.Conductivity;
end

if plotflag == true
    plottitle = "HYCOM Current Speed (2024-07-21) & SRV MarieTharp TSG Salinity";
    plotname = "NESMA2024_HYCOM_MT.png";

    # approximate x-axis scaling to make it look "normal"
    dlat = latmax - latmin;
    dlon = lonmax - lonmin;
    xfac = sind(90-lat0);
    yres = 1000;
    pres = (abs(ceil(yres * (xfac/(dlat/dlon)))), abs(yres));

    
    fig = Figure(size = pres, fontsize = 32)
    ax = Axis(
        fig[1, 1];
        title = plottitle,
        xlabel = "Longitude",
        ylabel = "Latitude",
    )

    if bathy_pf == true
        # plot bathymetry
        cmin_bathy, cmax_bathy = -6000, -1500;
        Makie.contourf!(x_bathy, y_bathy, z_bathy, colormap = ColorSchemes.gray1, levels = range(cmin_bathy, cmax_bathy, length = 128))
        xlims!(lonmin, lonmax);
        ylims!(latmin, latmax);
    end

    if hycom_pf == true
        # plot hycom data
        x_hycom = repeat(lon_hycom[lonind_hycom], 1, length(latind_hycom));
        y_hycom = repeat(lat_hycom[latind_hycom]', length(lonind_hycom), 1);
        x_hycom = ifelse.(x_hycom .> 180.0, x_hycom .- 360.0, x_hycom);
        #c_hycom = temp_hycom[:,:,1,1];
        c_hycom = sqrt.(u_hycom.^2 + v_hycom.^2);
        cmin_hycom, cmax_hycom = minimum(c_hycom[:]), maximum(c_hycom[:]);
        
        # Create a gradient with 256 colors based on the `:jet` colormap
        cmap = cgrad(ColorSchemes.tempo, 256);  # Create a gradient with 256 colors
        
        # Set desired alpha level (0.0 to 1.0)
        alpha_value = 0.5  # Set desired alpha level (0.0 to 1.0)
        
        # Create a new colormap with the added transparency
        cmap_transparent = [RGBA(c.r, c.g, c.b, alpha_value) for c in cmap.colors]

        Makie.contourf!(ax, x_hycom[:], y_hycom[:], c_hycom[:], levels = range(cmin_hycom, cmax_hycom, length = 16), colormap = cmap_transparent)  # Set alpha to make the fill semi-transparent
        cb_hycom = Colorbar(fig[1, 2], limits = (cmin_hycom, cmax_hycom), colormap = cmap, flipaxis = true, label="Speed (m/s)")
        cb_hycom.labelrotation = 3*π/2
    end

    if marietharp_pf == true
        # plot TSG data
        cmin_mt, cmax_mt = 32, 37;
        GLMakie.scatter!(lon1, lat1, color=salt1, colorrange=(cmin_mt,cmax_mt), colormap=:jet, markersize=10);
        GLMakie.scatter!(lon2, lat2, color=salt2, colorrange=(cmin_mt,cmax_mt), colormap=:jet, markersize=10);
        cb_mt = Colorbar(fig[1, 3], limits = (cmin_mt, cmax_mt), colormap = :jet, flipaxis = true, label="Salinity")
        cb_mt.labelrotation = 3*π/2 # Rotate by 90 degrees (π/2 radians)
    end

    fig
    save(figoutdir1 * plotname, fig)
    save(figoutdir2 * plotname, fig)
#    GLMakie.closeall()
end