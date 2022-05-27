module C2PO

using Dates, Missings
export gc_distance, rad2deg, deg2rad, histc, meshgrid, nan, findNaNmin, findNaNmax, nanfy, oneDize, datetimemissing2unixtimenan

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

end
