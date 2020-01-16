module GridFit

using SparseArrays
using Statistics
using StatsBase
using LinearAlgebra

import C2PO: nan, meshgrid, histc
export gridfit

# define parameter type
mutable struct Parameters
    smoothness::Array{Float64}
    interp::String
    regularizer::String
    solver::String
    maxiter::Int64
    tilesize::Float64
    overlap::Float64
    mask::AbstractArray
    autoscale::String
    xscale::Float64
    yscale::Float64
    maskflag::Bool
    extend::String
end

function gridfit(x::Array{Float64}, y::Array{Float64}, z::Array{Float64}, xnodes::Array{Float64}, ynodes::Array{Float64}; smoothness::Array{Float64}=[1.0 1.0], interp::String="triangle", regularizer::String="gradient", solver::String="normal", maxiter::Int64=10000, extend::String="warning", tilesize::Float64=Inf, overlap::Float64=0.2, autoscale::String="on")

    # gridfit: estimates a surface on a 2d grid, based on scattered data
    #          Replicates are allowed. All methods extrapolate to the grid
    #          boundaries. Gridfit uses a modified ridge estimator to
    #          generate the surface, where the bias is toward smoothness.
    #
    #          Gridfit is not an interpolant. Its goal is a smooth surface
    #          that approximates your data, but allows you to control the
    #          amount of smoothing.
    #
    # usage #1: zgrid = gridfit(x,y,z,xnodes,ynodes);
    # usage #2: [zgrid,xgrid,ygrid] = gridfit(x,y,z,xnodes,ynodes);
    # usage #3: zgrid = gridfit(x,y,z,xnodes,ynodes,prop,val,prop,val,...);
    #
    # Arguments: (input)
    #  x,y,z - vectors of equal lengths, containing arbitrary scattered data
    #          The only constraint on x and y is they cannot ALL fall on a
    #          single line in the x-y plane. Replicate points will be treated
    #          in a least squares sense.
    #
    #          ANY points containing a NaN are ignored in the estimation
    #
    #  xnodes - vector defining the nodes in the grid in the independent
    #          variable (x). xnodes need not be equally spaced. xnodes
    #          must completely span the data. If they do not, then the
    #          'extend' property is applied, adjusting the first and last
    #          nodes to be extended as necessary. See below for a complete
    #          description of the 'extend' property.
    #
    #          If xnodes is a scalar integer, then it specifies the number
    #          of equally spaced nodes between the min and max of the data.
    #
    #  ynodes - vector defining the nodes in the grid in the independent
    #          variable (y). ynodes need not be equally spaced.
    #
    #          If ynodes is a scalar integer, then it specifies the number
    #          of equally spaced nodes between the min and max of the data.
    #
    #          Also see the extend property.
    #
    #  Additional arguments follow in the form of property/value pairs.
    #  Valid properties are:
    #    'smoothness', 'interp', 'regularizer', 'solver', 'maxiter'
    #    'extend', 'tilesize', 'overlap'
    #
    #  Any UNAMBIGUOUS shortening (even down to a single letter) is
    #  valid for property names. All properties have default values,
    #  chosen (I hope) to give a reasonable result out of the box.
    #
    #   'smoothness' - scalar or vector of length 2 - determines the
    #          eventual smoothness of the estimated surface. A larger
    #          value here means the surface will be smoother. Smoothness
    #          must be a non-negative real number.
    #
    #          If this parameter is a vector of length 2, then it defines
    #          the relative smoothing to be associated with the x and y
    #          variables. This allows the user to apply a different amount
    #          of smoothing in the x dimension compared to the y dimension.
    #
    #          Note: the problem is normalized in advance so that a
    #          smoothness of 1 MAY generate reasonable results. If you
    #          find the result is too smooth, then use a smaller value
    #          for this parameter. Likewise, bumpy surfaces suggest use
    #          of a larger value. (Sometimes, use of an iterative solver
    #          with too small a limit on the maximum number of iterations
    #          will result in non-convergence.)
    #
    #          DEFAULT: 1
    #
    #
    #   'interp' - character, denotes the interpolation scheme used
    #          to interpolate the data.
    #
    #          DEFAULT: 'triangle'
    #
    #          'bilinear' - use bilinear interpolation within the grid
    #                     (also known as tensor product linear interpolation)
    #
    #          'triangle' - split each cell in the grid into a triangle,
    #                     then linear interpolation inside each triangle
    #
    #          'nearest' - nearest neighbor interpolation. This will
    #                     rarely be a good choice, but I included it
    #                     as an option for completeness.
    #
    #
    #   'regularizer' - character flag, denotes the regularization
    #          paradignm to be used. There are currently three options.
    #
    #          DEFAULT: 'gradient'
    #
    #          'diffusion' or 'laplacian' - uses a finite difference
    #              approximation to the Laplacian operator (i.e, del^2).
    #
    #              We can think of the surface as a plate, wherein the
    #              bending rigidity of the plate is specified by the user
    #              as a number relative to the importance of fidelity to
    #              the data. A stiffer plate will result in a smoother
    #              surface overall, but fit the data less well. I've
    #              modeled a simple plate using the Laplacian, del^2. (A
    #              projected enhancement is to do a better job with the
    #              plate equations.)
    #
    #              We can also view the regularizer as a diffusion problem,
    #              where the relative thermal conductivity is supplied.
    #              Here interpolation is seen as a problem of finding the
    #              steady temperature profile in an object, given a set of
    #              points held at a fixed temperature. Extrapolation will
    #              be linear. Both paradigms are appropriate for a Laplacian
    #              regularizer.
    #
    #          'gradient' - attempts to ensure the gradient is as smooth
    #              as possible everywhere. Its subtly different from the
    #              'diffusion' option, in that here the directional
    #              derivatives are biased to be smooth across cell
    #              boundaries in the grid.
    #
    #              The gradient option uncouples the terms in the Laplacian.
    #              Think of it as two coupled PDEs instead of one PDE. Why
    #              are they different at all? The terms in the Laplacian
    #              can balance each other.
    #
    #          'springs' - uses a spring model connecting nodes to each
    #              other, as well as connecting data points to the nodes
    #              in the grid. This choice will cause any extrapolation
    #              to be as constant as possible.
    #
    #              Here the smoothing parameter is the relative stiffness
    #              of the springs connecting the nodes to each other compared
    #              to the stiffness of a spting connecting the lattice to
    #              each data point. Since all springs have a rest length
    #              (length at which the spring has zero potential energy)
    #              of zero, any extrapolation will be minimized.
    #
    #          Note: The 'springs' regularizer tends to drag the surface
    #          towards the mean of all the data, so too large a smoothing
    #          parameter may be a problem.
    #
    #
    #   'solver' - character flag - denotes the solver used for the
    #          resulting linear system. Different solvers will have
    #          different solution times depending upon the specific
    #          problem to be solved. Up to a certain size grid, the
    #          direct \ solver will often be speedy, until memory
    #          swaps causes problems.
    #
    #          What solver should you use? Problems with a significant
    #          amount of extrapolation should avoid lsqr. \ may be
    #          best numerically for small smoothnesss parameters and
    #          high extents of extrapolation.
    #
    #          Large numbers of points will slow down the direct
    #          \, but when applied to the normal equations, \ can be
    #          quite fast. Since the equations generated by these
    #          methods will tend to be well conditioned, the normal
    #          equations are not a bad choice of method to use. Beware
    #          when a small smoothing parameter is used, since this will
    #          make the equations less well conditioned.
    #
    #          DEFAULT: 'normal'
    #
    #          '\' - uses matlab's backslash operator to solve the sparse
    #                     system. 'backslash' is an alternate name.
    #
    #          'symmlq' - uses matlab's iterative symmlq solver
    #
    #          'lsqr' - uses matlab's iterative lsqr solver
    #
    #          'normal' - uses \ to solve the normal equations.
    #
    #
    #   'maxiter' - only applies to iterative solvers - defines the
    #          maximum number of iterations for an iterative solver
    #
    #          DEFAULT: min(10000,length(xnodes)*length(ynodes))
    #
    #
    #   'extend' - character flag - controls whether the first and last
    #          nodes in each dimension are allowed to be adjusted to
    #          bound the data, and whether the user will be warned if
    #          this was deemed necessary to happen.
    #
    #          DEFAULT: 'warning'
    #
    #          'warning' - Adjust the first and/or last node in
    #                     x or y if the nodes do not FULLY contain
    #                     the data. Issue a warning message to this
    #                     effect, telling the amount of adjustment
    #                     applied.
    #
    #          'never'  - Issue an error message when the nodes do
    #                     not absolutely contain the data.
    #
    #          'always' - automatically adjust the first and last
    #                     nodes in each dimension if necessary.
    #                     No warning is given when this option is set.
    #
    #
    #   'tilesize' - grids which are simply too large to solve for
    #          in one single estimation step can be built as a set
    #          of tiles. For example, a 1000x1000 grid will require
    #          the estimation of 1e6 unknowns. This is likely to
    #          require more memory (and time) than you have available.
    #          But if your data is dense enough, then you can model
    #          it locally using smaller tiles of the grid.
    #
    #          My recommendation for a reasonable tilesize is
    #          roughly 100 to 200. Tiles of this size take only
    #          a few seconds to solve normally, so the entire grid
    #          can be modeled in a finite amount of time. The minimum
    #          tilesize can never be less than 3, although even this
    #          size tile is so small as to be ridiculous.
    #
    #          If your data is so sparse than some tiles contain
    #          insufficient data to model, then those tiles will
    #          be left as NaNs.
    #
    #          DEFAULT: inf
    #
    #
    #   'overlap' - Tiles in a grid have some overlap, so they
    #          can minimize any problems along the edge of a tile.
    #          In this overlapped region, the grid is built using a
    #          bi-linear combination of the overlapping tiles.
    #
    #          The overlap is specified as a fraction of the tile
    #          size, so an overlap of 0.20 means there will be a 20%
    #          overlap of successive tiles. I do allow a zero overlap,
    #          but it must be no more than 1/2.
    #
    #          0 <= overlap <= 0.5
    #
    #          Overlap is ignored if the tilesize is greater than the
    #          number of nodes in both directions.
    #
    #          DEFAULT: 0.20
    #
    #
    #   'autoscale' - Some data may have widely different scales on
    #          the respective x and y axes. If this happens, then
    #          the regularization may experience difficulties.
    #
    #          autoscale = 'on' will cause gridfit to scale the x
    #          and y node intervals to a unit length. This should
    #          improve the regularization procedure. The scaling is
    #          purely internal.
    #
    #          autoscale = 'off' will disable automatic scaling
    #
    #          DEFAULT: 'on'
    #
    #
    # Arguments: (output)
    #  zgrid   - (nx,ny) array containing the fitted surface
    #
    #  xgrid, ygrid - as returned by meshgrid(xnodes,ynodes)
    #
    #
    # Speed considerations:
    #  Remember that gridfit must solve a LARGE system of linear
    #  equations. There will be as many unknowns as the total
    #  number of nodes in the final lattice. While these equations
    #  may be sparse, solving a system of 10000 equations may take
    #  a second or so. Very large problems may benefit from the
    #  iterative solvers or from tiling.
    #
    #
    # Example usage:
    #
    #  x = rand(100,1);
    #  y = rand(100,1);
    #  z = exp(x+2*y);
    #  xnodes = 0:.1:1;
    #  ynodes = 0:.1:1;
    #
    #  g = gridfit(x,y,z,xnodes,ynodes);
    #
    # Note: this is equivalent to the following call:
    #
    #  g = gridfit(x,y,z,xnodes,ynodes, ...
    #              'smooth',1, ...
    #              'interp','triangle', ...
    #              'solver','normal', ...
    #              'regularizer','gradient', ...
    #              'extend','warning', ...
    #              'tilesize',inf);
    #
    #
    # Author: John D'Errico
    # e-mail address: woodchips@rochester.rr.com
    # Release: 2.0
    # Release date: 5/23/06

    testflag = 0
    if testflag == 1
        # create test dataset
        #x = [collect(1:0.5:100.0); 200.0];
        #y = deepcopy(x);
        #z = deepcopy(x);
        #x[10] = NaN;
        #y[30] = NaN;
        #z[70] = NaN;
        #xnodes = collect(1:100.0);
        #ynodes = collect(1:150.0);

        n1 = 15;
        n2 = 20;
        theta = rand(n1,1) .* pi/2;
        r = rand(1,n2);

        x = cos.(theta) * r;
        y = sin.(theta) * r;
        x=x[:];
        y=y[:];

        x = [Array{Float64,1}([0;0;1;1]); x; x; 1 .- x; 1 .- x];
        y = [Array{Float64,1}([0;1;0;1]); y; 1 .- y; y; 1 .- y];

        z = sin.(4 .* x .+ 5 .* y) .* cos.(7 .* (x .- y)) .+ exp.(x .+ y);

        xi = collect(0.0 : 1.0/200.0 : 1.0);
        yi = collect(0.0 : 1.0/100.0 : 1.0);
    end

    # ensure all of x,y,z,xnodes,ynodes are column vectors, also drop any NaN data
    x = x[:]; y = y[:]; z = z[:];
    k = findall((.!isnan.(x) .& .!isnan.(y) .& .!isnan.(z)) .== true);
    x = x[k]; y = y[k]; z = z[k];
    xmin = minimum(x);
    xmax = maximum(x);
    ymin = minimum(y);
    ymax = maximum(y);

    # did they supply a scalar for the nodes?
    if length(xnodes)==1
        xnodes = collect(xmin : (xmax-xmin)/xnodes : xmax);
        xnodes[end] = xmax; # make sure it hits the max
    end
    if length(ynodes)==1
        ynodes = collect(ymin : (ymax-ymin)/ynodes : ymax);
        ynodes[end] = ymax; # make sure it hits the max
    end

    xnodes=xnodes[:];
    ynodes=ynodes[:];
    dx = diff(xnodes);
    dy = diff(ynodes);
    nx = length(xnodes);
    ny = length(ynodes);
    ngrid = nx*ny;

    # set default parameters
    #params = Parameters([1.0 1.5], "triangle", "gradient", "backslash", 10000, Inf, 0.2, Array{Bool,2}(undef,nx,ny), "on", 1.0, 1.0, Bool(0), "always");
    #params.mask .= 1;
    params = Parameters([1.0 1.5], "triangle", "gradient", "backslash", 10000, Inf, 0.2, [], "on", 1.0, 1.0, Bool(0), "always");

    # set the scaling if autoscale was on
    if params.autoscale == "on"
        params.xscale = mean(dx);
        params.yscale = mean(dy);
        params.autoscale = "off";
    end

    if params.tilesize < max(nx,ny)
        # split it into smaller tiles. compute zgrid and ygrid
        # at the very end if requested
        #zgrid = tiled_gridfit(x,y,z,xnodes,ynodes,params);
        zgrid = [];
        return zgrid;
    else
        #it's a single tile.

        # mask must be either an empty array, or a boolean array of the same size as the final grid.
        nmask = size(params.mask);
        if !isempty(params.mask)
            if ((nmask[2] != nx) | (nmask[1] != ny))
                if ((nmask[2]==ny) | (nmask[1]==nx))
                    ErrorException("Mask array is probably transposed from proper orientation.")
                else
                    ErrorException("Mask array must be the same size as the final grid.")
                end
            end
        end
        if !isempty(params.mask)
            params.maskflag = 1;
        else
            params.maskflag = 0;
        end

        # default for maxiter?
        if isempty(params.maxiter)
            params.maxiter = min(10000,nx*ny);
        end

        # check lengths of the data
        n = length(x);
        if (length(y) != n) | (length(z) != n)
            ErrorException("Data vectors are incompatible in size.")
        end
        if n < 3
            ErrorException("Insufficient data for surface estimation.")
        end

        # verify the nodes are distinct
        if any(diff(xnodes) .<= 0) .| any(diff(ynodes) .<= 0)
            ErrorException("xnodes and ynodes must be monotone increasing")
        end

        # do we need to tweak the first or last node in x or y?
        if xmin<xnodes[1]
            if params.extend == "always"
                xnodes[1] = xmin;
            elseif param.extend == "warning"
                display("WARNING: GRIDFIT:extend xnodes[1] was decreased by: " * string(xnodes[1]-xmin) * ", new node = " * string(xmin))
                xnodes[1] = xmin;
            elseif param.extend == "never"
                ErrorException("Some x (" * string(xmin) * ") falls below xnodes[1] by: " * string(xnodes[1]-xmin))
            end
        end
        if xmax>xnodes[end]
            if params.extend == "always"
                xnodes[end] = xmax;
            elseif params.extend == "warning"
                display("WARNING: GRIDFIT:extend xnodes(end) was increased by: " * string(xmax-xnodes[end]) * ", new node = " * string(xmax))
                xnodes[end] = xmax;
            elseif params.extend == "never"
                ErrorException("Some x (" * string(xmax) * ") falls above xnodes[end] by: " * string(xmax-xnodes[end]))
            end
        end
        if ymin<ynodes[1]
            if params.extend == "always"
                ynodes[1] = ymin;
            elseif params.extend == "warning"
                display("WARNING: GRIDFIT:extend ynodes[1] was decreased by: " * string(ynodes[1]-ymin) * ", new node = " * string(ymin))
                ynodes[1] = ymin;
            elseif params.extend == "never"
                ErrorException("Some y (" * string(ymin) * ") falls below ynodes[1] by: " * string(ynodes[1]-ymin))
            end
        end
        if ymax>ynodes[end]
            if params.extend == "always"
                ynodes[end] = ymax;
            elseif params.extend == "warning"
                display("WARNING: GRIDFIT:extend ynodes[end] was increased by: " * string(ymax-ynodes[end]) * ", new node = " * string(ymax))
                ynodes[end] = ymax;
            elseif params.extend == "never"
                ErrorException("Some y (" * string(ymax) * ") falls above ynodes[end] by: " * string(ymax-ynodes[end]))
            end
        end

        # determine which cell in the array each point lies in
        #hx = fit(Histogram, x, xnodes);
        #hy = fit(Histogram, y, ynodes);
        indx = histc(x, xnodes);
        indy = histc(y, ynodes);

        # any point falling at the last node is taken to be inside the last cell in x or y.
        k = findall(indx .== nx);
        indx[k]=indx[k] .- 1;
        k = findall(indy .== ny);
        indy[k]=indy[k] .- 1;
        ind = indy .+ ny .* (indx .- 1);

        # Do we have a mask to apply?
        if params.maskflag
            # if we do, then we need to ensure that every
            # cell with at least one data point also has at
            # least all of its corners unmasked.
            params.mask[ind] .= 1;
            params.mask[ind .+ 1] .= 1;
            params.mask[ind .+ ny] .= 1;
            params.mask[ind .+ ny .+ 1] .= 1;
        end

        # interpolation equations for each point
        tx = min.(1.0,max.(0.0,(x .- xnodes[indx])./dx[indx]));
        ty = min.(1.0,max.(0.0,(y .- ynodes[indy])./dy[indy]));

        if params.interp == "triangle"
            # linear interpolation inside each triangle
            k = findall(tx .> ty);
            L = ones(n);
            L[k] .= ny;

            t1 = min(tx,ty);
            t2 = max(tx,ty);
            #A = sparse(repeat(collect(1:n),[1,3]...), [ind ind .+ ny .+ 1 ind .+ L],[1 .- t2 t1 t2 .- t1], n, ngrid);
            A = sparse(repeat(collect(1:n),[3]...), [ind; ind .+ ny .+ 1; ind .+ L],[1 .- t2; t1; t2 .- t1], n, ngrid);

        elseif params.interp == "nearest"
            # nearest neighbor interpolation in a cell
            k = round.(1 .- ty) + round.(1 .- tx) .* ny;
            A = sparse(collect(1:n), ind .+ k, ones(n), n, ngrid);

        elseif params.interp == "bilinear"
            # bilinear interpolation in a cell
            A = sparse(repeat(collect(1:n),[4]...), [ind; ind .+ 1; ind .+ ny .+ 1; ind .+ ny .+ 1], [(1 .- tx) .* (1 .- ty); (1 .- tx) .* ty; tx .* (1 .- ty); tx .* ty], n, ngrid);
            #A = sparse(repmat((1:n)',1,4),[ind,ind+1,ind+ny,ind+ny+1],[(1-tx).*(1-ty), (1-tx).*ty, tx.*(1-ty), tx.*ty], n,ngrid);
        end
        rhs = z;

        # do we have relative smoothing parameters?
        if length(params.smoothness) == 1
            # it was scalar, so treat both dimensions equally
            smoothparam = params.smoothness;
            xyRelativeStiffness = [1.0; 1.0];
        else
            # It was a vector, so anisotropy reigns.
            # I've already checked that the vector was of length 2
            smoothparam = sqrt(prod(params.smoothness));
            xyRelativeStiffness = params.smoothness[:] ./ smoothparam;
        end

        # Build regularizer. Add del^4 regularizer one day.
        if params.regularizer == "springs"
            # zero "rest length" springs
            #i = [i for i in 1:nx, j in 1:(ny-1)];
            #j = [j for i in 1:nx, j in 1:(ny-1)];
            (i,j) = meshgrid(1:nx,1:(ny-1));
            ind = j[:] .+ ny .* (i[:] .- 1);
            m = nx*(ny-1);
            stiffness = 1 ./ (dy ./ params.yscale);
            Areg = sparse(repeat(collect(1:m),[2]...), [ind; ind .+ 1], (xyRelativeStiffness[2] * stiffness[j[:]] .* [-1 1])[:], m, ngrid);
            #Areg = sparse(repmat((1:m)',1,2),[ind,ind+1], xyRelativeStiffness(2)*stiffness(j(:))*[-1 1], m,ngrid); #matlab version

            #i = [i for i in 1:(nx-1), j in 1:ny];
            #j = [j for i in 1:(nx-1), j in 1:ny];
            (i,j) = meshgrid(1:(nx-1),1:ny);
            ind = j[:] .+ ny .* (i[:] .- 1);
            m = (nx-1)*ny;
            stiffness = 1 ./ (dx ./ params.xscale);
            Areg = [Areg; sparse(repeat(collect(1:m),[2]...), [ind; ind .+ ny], (xyRelativeStiffness[1] .* stiffness[i[:]] .* [-1 1])[:], m, ngrid)];
            #Areg = [Areg;sparse(repmat((1:m)',1,2),[ind,ind+ny], xyRelativeStiffness(1)*stiffness(i(:))*[-1 1],m,ngrid)]; #matlab version

            #i = [i for i in 1:(nx-1), j in 1:(ny-1)];
            #j = [j for i in 1:(nx-1), j in 1:(ny-1)];
            (i,j) = meshgrid(1:(nx-1),1:(ny-1));
            ind = j[:] .+ ny .* (i[:] .- 1);
            m = (nx-1)*(ny-1);
            stiffness = 1 ./ sqrt.((dx[i[:]] ./ params.xscale ./ xyRelativeStiffness[1]) .^ 2 + (dy[j[:]] ./ params.yscale ./ xyRelativeStiffness[2]) .^ 2);
            Areg = [Areg; sparse(repeat(collect(1:m),[2]...), [ind; ind .+ ny .+ 1], (stiffness .* [-1 1])[:], m, ngrid)];
            #Areg = [Areg;sparse(repmat((1:m)',1,2),[ind,ind+ny+1], stiffness*[-1 1],m,ngrid)];
            Areg = [Areg; sparse(repeat(collect(1:m),[2]...), [ind .+ 1; ind .+ ny], (stiffness .* [-1 1])[:], m, ngrid)];
            #Areg = [Areg;sparse(repmat((1:m)',1,2),[ind+1,ind+ny], stiffness*[-1 1],m,ngrid)];

        elseif (params.regularizer == "diffusion") | (params.regularizer == "laplacian")
            # thermal diffusion using Laplacian (del^2)
            #i = [i for i in 1:nx, j in 2:(ny-1)];
            #j = [j for i in 1:nx, j in 2:(ny-1)];
            (i,j) = meshgrid(1:nx,2:(ny-1));
            ind = j[:] .+ ny .* (i[:] .- 1);
            dy1 = dy[j[:] .- 1] ./ params.yscale;
            dy2 = dy[j[:]] ./ params.yscale;
            Areg = sparse(repeat(ind,[3]...), [ind .- 1; ind; ind .+ 1], xyRelativeStiffness[2] .* [-2 ./ (dy1 .* (dy1 .+ dy2)); 2 ./ (dy1 .* dy2); -2 ./ (dy2 .* (dy1 .+ dy2))], ngrid, ngrid);
            #Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], xyRelativeStiffness(2)*[-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), -2./(dy2.*(dy1+dy2))],ngrid,ngrid);

            #i = [i for i in 2:(nx-1), j in 1:ny];
            #j = [j for i in 2:(nx-1), j in 1:ny];
            (i,j) = meshgrid(2:(nx-1),1:ny);
            ind = j[:] .+ ny .* (i[:] .- 1);
            dx1 = dx[i[:] .- 1] ./ params.xscale;
            dx2 = dx[i[:]] ./ params.xscale;
            Areg .= Areg .+ sparse(repeat(ind,[3]...), [ind .- ny; ind; ind .+ ny], xyRelativeStiffness[1] .* [-2 ./ (dx1 .* (dx1 .+ dx2)); 2 ./ (dx1 .* dx2); -2 ./ (dx2 .* (dx1 .+ dx2))], ngrid, ngrid);
            #Areg = Areg + sparse(repmat(ind,1,3),[ind-ny,ind,ind+ny], xyRelativeStiffness(1)*[-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), -2./(dx2.*(dx1+dx2))],ngrid,ngrid);

        elseif params.regularizer == "gradient"
            # Subtly different from the Laplacian. A point for future
            # enhancement is to do it better for the triangle interpolation case.
            #i = [i for i in 1:nx, j in 2:(ny-1)];
            #j = [j for i in 1:nx, j in 2:(ny-1)];
            (i,j) = meshgrid(1:nx,2:(ny-1));
            ind = j[:] .+ ny .* (i[:] .- 1);
            dy1 = dy[j[:] .- 1] ./ params.yscale;
            dy2 = dy[j[:]] ./ params.yscale;

            Areg = sparse(repeat(ind,[3]...), [ind .- 1; ind; ind .+ 1], xyRelativeStiffness[2] .* [-2 ./ (dy1 .* (dy1 .+ dy2)); 2 ./ (dy1 .* dy2); -2 ./ (dy2 .* (dy1 .+ dy2))], ngrid, ngrid);
            #Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], xyRelativeStiffness(2)*[-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), -2./(dy2.*(dy1+dy2))],ngrid,ngrid);

            #i = [i for i in 2:(nx-1), j in 1:ny];
            #j = [j for i in 2:(nx-1), j in 1:ny];
            (i,j) = meshgrid(2:(nx-1),1:ny);
            ind = j[:] .+ ny .* (i[:] .- 1);
            dx1 = dx[i[:] .- 1] ./ params.xscale;
            dx2 = dx[i[:]] ./ params.xscale;
            Areg = [Areg; sparse(repeat(ind,[3]...), [ind .- ny; ind; ind .+ ny], xyRelativeStiffness[1] .* [-2 ./ (dx1 .* (dx1 .+ dx2)); 2 ./ (dx1 .* dx2); -2 ./ (dx2 .* (dx1 .+ dx2))], ngrid, ngrid)];
            #Areg = [Areg;sparse(repmat(ind,1,3),[ind-ny,ind,ind+ny], xyRelativeStiffness(1)*[-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), -2./(dx2.*(dx1+dx2))],ngrid,ngrid)];
        end
        nreg = size(Areg,1);

        # Append the regularizer to the interpolation equations,
        # scaling the problem first. Use the 1-norm for speed.
        NA = norm(A,1);
        NR = norm(Areg,1);
        A = [A; Areg .* (smoothparam*NA/NR)];
        rhs = [rhs; zeros(nreg,1)];
        # do we have a mask to apply?
        if params.maskflag == true
            unmasked = findall(params.mask[:]);
            unmaskedCI = findall(params.mask);
        end
        # solve the full system, with regularizer attached
        if (params.solver == "backslash") | (params.solver == "\\")
            if params.maskflag
                # there is a mask to use
                zgrid=Array{Float64,2}(undef,ny,nx);
                zgrid .= NaN;
                zgrid[unmaskedCI] = A[:,unmasked] \ rhs;
            else
                # no mask
                zgrid = reshape(A \ rhs, (ny, nx));
            end
        elseif params.solver == "normal"
            # The normal equations, solved with \. Can be faster
            # for huge numbers of data points, but reasonably
            # sized grids. The regularizer makes A well conditioned
            # so the normal equations are not a terribly bad thing
            # here.
            if params.maskflag
                # there is a mask to use
                Aunmasked = A[:,unmasked];
                zgrid=Array{Float64,2}(undef,ny,nx);
                zgrid .= NaN;
                zgrid[unmaskedCI] = (Aunmasked' * Aunmasked) \ (Aunmasked' * rhs);
            else
                zgrid = reshape((A' * A) \ (A' * rhs), (ny, nx));
            end
        end  # if params.solver
        return zgrid;
    end

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

# emulate the behavior of Matlab's nan function
function nan(m::Unsigned,n::Unsigned)
    nanarray=Array{Float64,2}(undef,m,n);
    nanarray .= NaN;
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
    i = [i for j in ygrid[1]:ygrid[end], i in xgrid[1]:xgrid[end]];
    j = [j for j in ygrid[1]:ygrid[end], i in xgrid[1]:xgrid[end]];
    return (i,j)
end

end # module GridFit
