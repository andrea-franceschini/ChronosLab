 function [alpha, beta, gamma] = cpt_LSpol(intv, iter) 
%%--------------------- 
 p0 = initpol(intv);
 ppol = p0; 
 ppol = xmul(ppol, intv);
%%--------------------  
 bet1      = sqrt(dotp(ppol,ppol));
 ppol      = ppol ./ bet1;
 gamma(1)  = dotp(ppol,p0);
 bet = 0.0;
 beta(1) = bet1;
 for i=1:iter
     appol = xmul(ppol,intv); 
     alp = dotp(appol, ppol) ;
     alpha(i) = alp;
     appol = polsum(appol,-alp, ppol);
     if (i >1)
        appol = polsum(appol, -bet, qpol);
     end
     bet = sqrt(dotp(appol,appol));
     beta(i+1) = bet;
     if (bet == 0) 
        error(' division by zero') 
     end
     qpol = ppol;
     ppol = appol/bet;
     gamma(i+1)  = dotp(ppol,p0);
  end
