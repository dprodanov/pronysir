function zz=sirapprox3(x, g, bb)
% double exponential fit
% g - shape
% bb - peak value
 cc=(1-exp(-bb*x))*g/bb;
 zz=bb*exp(cc-g*x);
end