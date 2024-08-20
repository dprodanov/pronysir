function [c, R2] = rsapfit(X, Y, sircoeffs)
% c(1) - shape
% c(2) - peak
% c(3) - center
% c(4) - scaling

dx=max(X)- min(X)
ss=max(Y)- min(Y)
g0=1/max(X)
c = fminsearchbnd(@(Params) costfun(Params, X, Y), sircoeffs, [g0; eps; -dx; 1], [Inf;  Inf; dx; Inf]);
        
R2 = sqrt(sum((c(4)* rapprox(X -c(3), c(1), c(2)) - Y).^2));
n=length(Y)-    length(sircoeffs);
R2=1- R2/n/var(Y);
end     
        
function f = costfun(Params, X, Y)
R2 =  sum((Params(4)* rapprox(X -Params(3), Params(1), Params(2)) - Y).^2);
f = R2;
end