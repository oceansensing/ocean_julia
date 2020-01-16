# include("run_npp_load_ctd.jl")
#
# Revision: 2019-11-19: DG

if (pwd() in LOAD_PATH) == false
    push!(LOAD_PATH, pwd());
end

import npp_load_ctd: ctd, npp
import npps: transect, npp
import WOA: npp_load_woa, load_woa, season
profile = npp_load_woa(3);
