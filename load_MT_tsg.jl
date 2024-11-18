using Dates, CSV, DataFrames, NaNMath, NCDatasets, GLMakie

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

bathypath = "/Users/gong/oceansensing Dropbox/C2PO/Data/bathy/ETOPO1/ETOPO_2022_v1_30s_N90W180_surface.nc";
bathyds = Dataset(bathypath,"r");

lon = bathyds["lon"][:];
lat = bathyds["lat"][:];

latmin, latmax = 37, 40.5;
lonmin, lonmax = -66.0, -59.8;

# approximate x-axis scaling to make it look "normal"
dlat = latmax - latmin;
dlon = lonmax - lonmin;
lat0 = NaNMath.mean(lat1);
xfac = sind(90-lat0);
yres = 2000;
pres = (abs(ceil(yres * (xfac/(dlat/dlon)))), abs(yres));

# extract indices from the bathymetric data file
latind = findall(latmin-0.1 .<= lat .<= latmax+0.1);
lonind = findall(lonmin-0.1 .<= lon .<= lonmax+0.1);

z = Float64.(bathyds["z"][lonind, latind]); # recasting as Float64 to fix a StackOverFlow error seen in GLMakie 0.6.0.
x = lon[lonind];
y = lat[latind];

pzind = findall(z .> 1);
nzind = findall(z .< -1);
zzind = findall(-1 .<= z .<= 1);

logzflag = 0
if logzflag == 1
    log10z = deepcopy(z);
    log10z[pzind] .= log10.(z[pzind]);
    log10z[nzind] .= -log10.(-z[nzind]);
    log10z[zzind] .= 0;
    zrange = (-4, 4);
    zp = log10z;
else
    zp = z;
    zrange = (-6000, 6000);
end

plottitle = "SRV Marie Tharp TSG - NESMA 2024";
plotname = "NESMA2024_MarieTharp_TSG.png";

fig = Figure(size = pres, fontsize = 64)
ax = Axis(
    fig[1, 1];
    title = plottitle,
    xlabel = "Longitude",
    ylabel = "Latitude",
)

Makie.contourf!(x, y, zp, colormap = :bukavu, levels = range(zrange[1], zrange[2], length = 128))
xlims!(lonmin, lonmax);
ylims!(latmin, latmax);

#lon, lat, c = lon1, lat1, salt1;
cmin, cmax = 32, 37;
GLMakie.scatter!(lon1, lat1, color=salt1, colorrange=(cmin,cmax), colormap=:jet, markersize=20);
GLMakie.scatter!(lon2, lat2, color=salt2, colorrange=(cmin,cmax), colormap=:jet, markersize=20);
Colorbar(fig[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = true, label="Sigma0")


fig
save(figoutdir * plotname, fig)
GLMakie.closeall()
