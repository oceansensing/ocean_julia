module NPP_types

using Dates
export CTD, Expedition, Transect

# defining the CTD data type for storing hydrographic data
mutable struct CTD
    yday::Array{Float64}
    lon::Array{Float64}
    lat::Array{Float64}
    p::Array{Float64}
    z::Array{Float64}
    temp::Array{Float64}
    ptmp::Array{Float64}
    ctmp::Array{Float64}
    salt::Array{Float64}
    saltA::Array{Float64}
    rho::Array{Float64}
    sigma0::Array{Float64}
    oxy::Array{Float64}
    oxysat::Array{Float64}
    fluor::Array{Float64}
    casttype::Array{Int64}
    pind::Array{Int64}
    dataflag::Array{Int64}
end

# define the Expedition data type for storing cruise metadata
struct Expedition
    name::AbstractString
    chiefscientist::AbstractString
    section::AbstractArray
    sectname::Array{AbstractString}
    begindate::Date
    enddate::Date
    s::AbstractArray
end

mutable struct Transect
    number::Integer
    name::AbstractString
    x::Array{<:AbstractFloat,1}
    y::Array{<:AbstractFloat,1}
    ctmp::Array{<:AbstractFloat,1}
    saltA::Array{<:AbstractFloat,1}
    sigma0::Array{<:AbstractFloat,1}
    xgrid::Array{<:AbstractFloat,2}
    ygrid::Array{<:AbstractFloat,2}
    ctmpgrid::Array{<:AbstractFloat,2}
    saltAgrid::Array{<:AbstractFloat,2}
    sigma0grid::Array{<:AbstractFloat,2}
    s::Array{<:AbstractFloat,1}
    ymin::Array{<:AbstractFloat,1}
end

end
