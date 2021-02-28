using PyCall
cdsapi = pyimport("cdsapi")



c = cdsapi.Client()

c.retrieve("reanalysis-era5-pressure-levels",
    Dict("variable" => "temperature", 
        "pressure_level" => "1000", 
        "product_type" => "reanalysis", 
        "year" => "2008",
        "month" => "01",
        "day" => "01",
        "time" => "12:00",
        "format" => "netcdf"),
    "era5output.netcdf"
)
