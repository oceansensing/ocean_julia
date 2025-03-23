module C2PO

using Distributed
@everywhere begin
    using Dates, Missings
end

export missing2nan, datetime2yearday, yearday2datetime, unix2yearday, yearday2unix, pyrow2jlcol, gc_distance, rad2deg, deg2rad, histc, meshgrid, nan, findNaNmin, findNaNmax, nanfy, oneDize, datetimemissing2unixtimenan, cdnlp, stresslp, intersectalajulia2, intersectalajulia4, runningavg!

function missing2nan(varin::AbstractFloat)
    varin = collect(varin);
    if (typeof(varin) == Vector{Union{Missing, Int64}}) | (typeof(varin) == Matrix{Union{Missing, Int64}}) | (typeof(varin) == Vector{Union{Missing, Int32}}) | (typeof(varin) == Matrix{Union{Missing, Int32}}) | (typeof(varin) == Vector{Union{Missing, Int16}}) | (typeof(varin) == Matrix{Union{Missing, Int16}})
        varout = Array{Float64}(undef,size(collect(varin)));
        varintypes = typeof.(varin);
        notmissind = findall(varintypes .!= Missing);
        missind = findall(varintypes .== Missing); 
        if isempty(notmissind) != true  
            varout[notmissind] .= Float64.(varin[notmissind]);
        end
        if isempty(missind) != true
            varout[missind] .= NaN;
        end
    elseif (typeof(varin) == Vector{Union{Missing, Float64}}) | (typeof(varin) == Vector{Union{Missing, Float32}}) | (typeof(varin) == Matrix{Union{Missing, Float64}}) | (typeof(varin) == Matrix{Union{Missing, Float32}}) | (typeof(varin) == Array{Union{Missing, Float32}, 4}) | (typeof(varin) == Array{Union{Missing, Float64}, 4}) | (typeof(varin) == Array{Union{Missing, Float32}, 3}) | (typeof(varin) == Array{Union{Missing, Float32}, 3})
        varout = Float64.(collect(Missings.replace(varin, NaN)));
    elseif (typeof(varin) == Vector{Missing}) | (typeof(varin) == Matrix{Missing}) | (typeof(varin) == Array{Missing, 3}) | (typeof(varin) == Array{Missing, 4})
        varout = Array{Float64}(undef,size(collect(varin)));
        varout .= NaN; 
    else
        varout = varin;
    end

    return varout
end

# DG 2024-09-06, adopted from seaexplorer
function datetime2yearday(xdt::DateTime)
    year = Dates.year(xdt);
    yday = Dates.dayofyear(xdt);
    seconds_in_day = 86400;
    ydayfrac = yday + (Dates.hour(xdt) * 3600 .+ Dates.minute(xdt) * 60 .+ Dates.second(xdt)) ./ seconds_in_day;
    return ydayfrac;
end

# DG 2024-09-06, with help from ChatGPT 4o
function yearday2datetime(yyyy::Int, yearday::Float64)
    # Separate integer part (day) and fractional part (time of day)
    int_day = Int(floor(yearday))  # Get the integer part of the day
    fractional_day = yearday - int_day  # Get the fractional part of the day

    # Convert fractional day into hours, minutes, and seconds
    total_seconds = Int(round(fractional_day * 86400))  # 86400 seconds in a day
    hours = div(total_seconds, 3600)  # Number of full hours
    minutes = div(total_seconds % 3600, 60)  # Number of full minutes
    seconds = total_seconds % 60  # Remaining seconds

    # Start from January 1st of the given year and add the number of days
    dateint = Date(yyyy, 1, 1) + Day(int_day - 1)  # Subtract 1 because days start from 1

    # Create the final DateTime with the computed hours, minutes, and seconds
    return DateTime(year(dateint), month(dateint), day(dateint), hours, minutes, seconds)
    #return DateTime(date) + Time(hours, minutes, seconds)
end

function unix2yearday(unixt::Float64)
    return datetime2yearday(unix2datetime(unixt));
end

function yearday2unix(yyyy::Int, yearday::Float64)
    return unix2datetime(yearday2datetime(yyyy, yearday));
end

# DG 2024-11-13, code provided by Grok2
function datenum2datetime(matlab_datenum::AbstractFloat)
    julia_datetime = DateTime(0) + Dates.Millisecond(round(Int, (matlab_datenum - 1) * 24 * 60 * 60 * 1000));
    return julia_datetime;
end


function pyrow2jlcol(invar::Matrix{Float64})
    return reverse(rotr90(invar), dims = 2);
end

function gc_distance(lat1deg::Float64,lon1deg::Float64,lat2deg::Float64,lon2deg::Float64)
    # This code implements Vincenty 1975: https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
    # In the future, need to implement: Karney 2013 https://arxiv.org/abs/1109.4448

    a = 6378.1370; # Earth semi major axis in km
    b = 6356.7523142; # Earth semi-minor axis in km
    R1 = 1/3 * (2*a + b);

    f = (a - b) / a; # flattening parameter for an elipsoid
    #f = 1/298.257223563; # defined by WGS84

    lat1 = deg2rad(lat1deg);
    lon1 = deg2rad(lon1deg);
    lat2 = deg2rad(lat2deg);
    lon2 = deg2rad(lon2deg);

    beta1 = atan((1-f)*tan(lat1)); # reduced latitude 1
    beta2 = atan((1-f)*tan(lat2)); # reduced latitude 2
    P = (beta1 + beta2)/2;
    Q = (beta2 - beta1)/2;

    lambda = abs(lon1 - lon2);
    sigma = atan(sqrt((cos(beta2)*sin(lambda))^2 + (cos(beta1)*sin(beta2) - sin(beta1)*cos(beta2)*cos(lambda))^2) / (sin(beta1)*sin(beta2) + cos(beta1)*cos(beta2)*cos(lambda)));

    alpha = asin(cos(beta1) * cos(beta2) * sin(lambda)/sin(sigma));
    sigmam = acos(cos(sigma) - 2*sin(beta1)*sin(beta2) / cos(alpha)^2)/2;


    X = (sigma - sin(sigma)) * (sin(P)^2 * cos(Q)^2) / (cos(sigma/2)^2)
    Y = (sigma + sin(sigma)) * (cos(P)^2 * sin(Q)^2) / (sin(sigma/2)^2)

    dist = a*(sigma - f/2 * (X+Y));
    return dist
end

function rad2deg(rad)
    return rad * 180.0/pi;
end

function deg2rad(deg)
    return deg * pi/180.0;
end

# define a histogram index function similar to matlab histc
function histc(x, xnodes)
    bin = Array{Number}(undef,length(x));
    for i = 1:length(x)
        a = findlast(x[i] .>= xnodes);
        if !isnothing(a)
            bin[i] = a;
        end
    end
    return bin;
end

# emulate the behavior of Matlab's meshgrid
function meshgrid(xgrid::Array{<:AbstractFloat,1},ygrid::Array{<:AbstractFloat,1})
    nx = length(xgrid);
    ny = length(ygrid);
    minx = minimum(xgrid);
    maxx = maximum(xgrid);
    miny = minimum(ygrid);
    maxy = maximum(ygrid);
    dx = (maxx .- minx) ./ (nx-1);
    dy = (maxy .- miny) ./ (ny-1);
    i = [i for j in miny:dy:maxy, i in minx:dx:maxx];
    j = [j for j in miny:dy:maxy, i in minx:dx:maxx];
    return (i,j)
end

function meshgrid(xgrid::Array{<:Integer,1},ygrid::Array{<:Integer,1})
    minx = minimum(xgrid);
    maxx = maximum(xgrid);
    miny = minimum(ygrid);
    maxy = maximum(ygrid);
    i = [i for j in miny:maxy, i in minx:maxx];
    j = [j for j in miny:maxy, i in minx:maxx];
    return (i,j)
end

function meshgrid(xgrid::UnitRange{<:Integer},ygrid::UnitRange{<:Integer})
    #minx = minimum(xgrid);
    #maxx = maximum(xgrid);
    #miny = minimum(ygrid);
    #maxy = maximum(ygrid);
    i = [i for i in xgrid[1]:xgrid[end], j in ygrid[1]:ygrid[end]];
    j = [j for i in xgrid[1]:xgrid[end], j in ygrid[1]:ygrid[end]];
    return (i,j)
end

# emulate the behavior of Matlab's nan function
function nan(m::Unsigned,n::Unsigned)
    nanarray=Array{Float64,2}(undef,m,n);
    nanarray .= NaN;
end

# adopted from https://github.com/mlubin/NaNMath.jl/issues/31
function findNaNmax(x::Array{<:AbstractFloat,1})
	result = convert(eltype(x), NaN)
    indmax = 0
    @inbounds @simd for i in eachindex(x)
    	v = x[i]
        if !isnan(v)
            if (isnan(result) || v > result)
                result = v
                indmax = i
            end
        end
    end
    return (result,indmax) # note the order of result and ind may be different that what's proposed by floswald as of 2019-12-03
end

# adopted from https://github.com/mlubin/NaNMath.jl/issues/31
function findNaNmin(x::Array{<:AbstractFloat,1})
	result = convert(eltype(x), NaN)
    indmin = 0
    @inbounds @simd for i in eachindex(x)
    	v = x[i]
        if !isnan(v)
            if (isnan(result) || v < result)
                result = v
                indmin = i
            end
        end
    end
    return (result,indmin) # note the order of result and ind may be different that what's proposed by floswald as of 2019-12-03
end

function dg_pol2cart(magnitude::AbstractFloat, compassdir::AbstractFloat)
    theta = mod.(360.0 .- compassdir .+ 90.0, 360) .* pi/180.0;
    u = magnitude .* cos.(theta);
    v = magnitude .* sin.(theta);
    return u + v*im
end

# this does not work, left in for reference. can you 'repeat' function with splat or Tuple() to reproduce Matlab's repmat behavior
#function repmat(x, shape)
#    xsize = size(x);
#    xout = reshape(repeat(x, outer=prod(shape)),[collect(size(x));collect(shape)]...);
#    return xout;
#end

function nanfy(datawithmissing)
    collect(Missings.replace(datawithmissing, NaN));
end

# oneDize turns multidimensional array into 1D array
function oneDize(nddata)
    return collect(reshape(nddata, prod(size(nddata))));
end

# datetimemissing2unixtimenan converts arrays with type of Union{Missing,DateTime} to Float64 with NaN for fillvalue
function datetimemissing2unixtimenan(datetimemissing::Array{Union{Missing, DateTime}})
    unixtimenan = Array{Float64}(undef,size(datetimemissing));
    gind = findall(ismissing.(datetimemissing) .== false);
    bind = findall(ismissing.(datetimemissing) .== true);
    unixtimenan[gind] .= Dates.datetime2unix.(disallowmissing(datetimemissing[gind]));
    unixtimenan[bind] .= NaN;
    return unixtimenan
end

function datetimemissing2unixtimenan(datetimemissing::Array{DateTime})
    unixtimenan = Array{Float64}(undef,size(datetimemissing));
    unixtimenan .= Dates.datetime2unix.(datetimemissing);
    return unixtimenan
end

function datetimemissing2unixtimenan(datetimemissing::Array{Missing})
    unixtimenan = Array{Float64}(undef,size(datetimemissing));
    unixtimenan .= NaN;
    return unixtimenan
end

function runningavg!(avgfield, nowfield, i, n)
    avgfield = ((i-1)/n * avgfield + nowfield/n) / (i/n);
    return avgfield;
end

function stresslp(sp, z, rhoa=1.22)
    # stresslp: computes neutral wind stress following Large and Pond (1981). tau = stresslp(sp,z,rhoa) computes the neutral wind stress given the wind speed at height z following Large and Pond (1981), J. Phys. Oceanog., 11, 324-336. Air density is an optional input, otherwise assumed to be constant (1.22 kg/m^3). 

    #INPUT:   sp    - wind speed             [m/s]
    #         z     - measurement height     [m]
    #         rhoa  - air_density (optional) [kg/m^3]

    #OUTPUT:  tau   - wind stress            [N/m^2]

    # 3/8/97: version 1.0
    # 8/26/98: version 1.1 (revised by RP)
    # 4/2/99: version 1.2 (optional air density added by AA)
    # 8/5/99: version 2.0
    # 3/5/23: ported to Julia by Donglai Gong

    include("as_consts.jl");
    if (@isdefined rho_air) == false
        rhoa = rho_air; 
    end
    (cd, u10) = cdnlp(sp, z);
    tau=rhoa .* (cd .* u10 .^ 2);
    return tau;
end

function cdnlp(sp, z)
    # cdnlp: computes neutral drag coefficient following Large&Pond (1981).
    # [cd,u10]=cdnlp(sp,z) computes the neutral drag coefficient and wind speed at 10m given the wind speed at height z following Large and Pond (1981), J. Phys. Oceanog., 11, 324-336. 
    
    # INPUT:   sp - wind speed  [m/s]
    #          z - measurement height [m]
    
    # OUTPUT:  Cd - neutral drag coefficient at 10m
    #          u10 - wind speed at 10m  [m/s]
    
    
    # 3/8/97: version 1.0
    # 8/26/98: version 1.1 (vectorized by RP)
    # 8/5/99: version 2.0
    # 3/5/2023: ported to Julia by Donglai Gong
    
    include("as_consts.jl"); # define physical constants
    
    if typeof(sp) == Float64
        sp = [sp];
    end

    a = log.(z ./ 10.0) ./ kappa;  # log-layer correction factor
    tol = 0.001;            # tolerance for iteration [m/s]
    u10o = zeros(size(sp));
    Cd = 1.15e-3 .* ones(size(sp));
    u10 = sp ./ (1 .+ a .* sqrt.(Cd));
    ii = abs.(u10 .- u10o) .> tol; 
    while any(ii)
        u10o = u10;
        Cd=4.9e-4 .+ 6.5e-5 .* u10o;    # compute cd(u10)
        smallwindi = findall(u10o .< 10.15385);
        Cd[smallwindi] .= 1.15e-3;
        u10 = sp ./ (1 .+ a .* sqrt.(Cd));   # next iteration
        ii = abs.(u10 .- u10o) .> tol;      # keep going until iteration converges
    end
    return (Cd, u10);
end    

# https://discourse.julialang.org/t/indices-of-intersection-of-two-arrays/23043/20
function intersectalajulia2(a,b)
    ia = findall(in(b), a)
    ib = findall(in(view(a,ia)), b)
    return unique(view(a,ia)), ia, ib[indexin(view(a,ia), view(b,ib))]
end

function intersectalajulia4(a,b)
    ab=intersect(a,b)
    ia = [findall(==(e), a) for e in ab]
    ib = [findall(==(e), b) for e in ab]
    return hcat(ab, ia,ib)
end

# linearly scale between xmin and xmax using the range provided by xscl
function minmaxscaler(x, xscl; zeromin::Bool=false)
    if zeromin == false
        xmin = minimum(skipmissing(x));
    else
        xmin = 0.0;
    end

    xmax = maximum(skipmissing(x));
    if xmin != xmax
        x_std = (x .- xmin) ./ (xmax - xmin);
        x_scl = x_std .* (maximum(xscl) - minimum(xscl)) .+ minimum(xscl);
    else
        x_scl = x;
    end
    return x_scl, xmin, xmax;
end

# linearly scale between xmin and xmax using the range provided in xscl
function minmaxscaler(x, xscl, xminmax)
    xmin = minimum(skipmissing(xminmax));
    xmax = maximum(skipmissing(xminmax));
    if minimum(skipmissing(x)) != maximum(skipmissing(x))
        x_std = (x .- xmin) ./ (xmax - xmin);
        x_scl = x_std .* (maximum(xscl) - minimum(xscl)) .+ minimum(xscl);
    else
        x_scl = x;
    end
    return x_scl, xmin, xmax;
end

# find standard deviation of a variable, then scale to +/- n std. dev.
function varscaler(varin,nstd=2.0)
    varstd = std(skipmissing(varin));
    varmean = mean(skipmissing(varin));
    varscl = (varin .- varmean) ./ (nstd * varstd);
    return varscl, varmean, varstd, nstd;
end

# manually scale the input variable using specified mean and standard deviation
function varscaler_manual(varin, varmean, varstd, nstd=2.0)
    return (varin .- varmean) ./ (nstd * varstd);
end

# use specified mean and standard deviation to 'unscale' the input variable
function varunscaler(varscl, varmean, varstd, nstd)
    return varscl .* (nstd * varstd) .+ varmean;
end

# this function returns the indices for tt within a specified time period
function findtind(tt, trange)
    return findall(minimum(trange) .<= tt .<= maximum(trange));
end

end