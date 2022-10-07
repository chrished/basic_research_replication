function hcond = hfunc(h)
  global alg p eq

  m = alg.m;

  if (m == p.M)
    hp = h;
  else
    hp = eq.hvec(m+1);
  end

  rhoa = p.rhoa(m);
  rhob = p.rhob(m);

  za = max(0.0,(1.0+p.rhoa(m))*(eq.abeta+h)/eq.w);
  zb = max(0.0,(1.0+p.rhob(m))*(eq.bbeta+h)/eq.w);
  zm = max(0.0,hp-h);

  if (alg.nobasic == 1)
    zb = 0.0;
  end

  xa = (za/(p.asigma*p.agamma))^(1.0/(p.agamma-1.0));
  xb = (zb/(p.bsigma*p.bgamma))^(1.0/(p.bgamma-1.0));

  % limit overall innovation rate. this is for stability
  if (alg.rlim == 1)
    gsum = (1.0+p.rhoa(m))*xa + (1.0+p.rhob(m))*xb;
    glim = eq.tau + alg.glimfact*(eq.r-eq.g);

    if (glim <= 0.0)
      xa = 0.0;
      xb = 0.0;
    elseif (gsum > glim)
      xa = xa*(glim/gsum);
      xb = xb*(glim/gsum);
    end
  end

  % labor cost
  ca = p.asigma*xa^(p.agamma);
  cb = p.bsigma*xb^(p.bgamma);

  % probability of positive basic research
  dstar = xb*(1.0+rhob)*(eq.bbeta+h)/(eq.w*p.bfcoeff) - cb;
  dstart = (log(max(0.0,dstar)) - p.dmu)/p.dsig;
  epb = normcdf(dstart);
  if (epb == 0.0)
    fb = 0.0;
  else
    fb = p.bfcoeff*exp(p.dmu+0.5*p.dsig^2.0)*normcdf(dstart-p.dsig)/epb;
  end

  pval = eq.w*(xa*za-ca + epb*(xb*zb-cb-fb)) + p.kappa*eq.ebeta + eq.freebase*eq.bbeta + eq.mrate*p.eprob(m)*zm;
  hcond = h*eq.denomh - pval;
end
