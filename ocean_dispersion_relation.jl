L = 1:1000; d = 1:1000;
Lmatrix = Array{Float64}(undef, length(L), length(d));
c(L,d) = (10*L/(2*pi)*tanh(2*pi*d/L))^0.5;

Lmatrix = [fill(convert(Float64,i),(1,length(d))) for i in 1:length(d)]
