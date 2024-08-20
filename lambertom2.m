%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIR Outbreak Model analysis package
%
% SIROM
% (c) Dimiter Prodanov, 2020
% Lambert Omega (x \pm i pi)
% k=0 - principlal branch, k=-1  - non principal branch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result=lambertom2(k,x)
 
 n=length(x);
 result=zeros(1,n);
 for i=1:n 
    if x(i)==-1 
        result (i)=1;
    elseif x(i)~=-1
        if k==0
            result(i)=principalh(x(i));
        elseif k==-1
            result(i)=nonprincipalh(x(i));
        end
    end
 end
end


function w= principalh(x)
   niter=16;
   tol=5e-16;
   eps=1e-32;
   w=exp(x);
   if abs(w)<eps 
       w=eps;
   else
      i=0;
      err=2;
      while i<niter && abs(err) >tol
          lw=log(abs(w));
         if ~isfinite(lw)
            lw=1/tol;
         end
         dw= x+w- lw;
         q= -dw./(2.*w) +(1.-w).^2./w;
         if q~= 0
             w=w+dw.*(1.-w)./q;
         else
             w=w+w.*dw./(1.-w);
         end
         err=lw-w-x;
         i=i+1;
      end
   end
end


function w= nonprincipalh(x)
   niter=16;
   tol=5e-16;
   w=-x +log(-x);
   i=0;
   err=abs(1-x);
   while i<niter && abs(err) >tol
        lw=log(abs(w));
        if ~isfinite(lw)
            lw=1/tol;
        end
        dw= x+w- lw;
        q= -dw./(2.*w) +(1.-w).^2./w;
        if q~= 0
             w=w+dw.*(1.-w)./q;
        else
             w=w+w.*dw./(1.-w);
        end
        err=lw-w-x;
        i=i+1;
   end
   
end

