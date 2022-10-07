function innovation
  global alg p eq

  % gains from actions
  eq.za = max(0.0,(1.0+p.rhoa).*(eq.abeta+eq.hvec)/eq.w);
  eq.zb = max(0.0,(1.0+p.rhob).*(eq.bbeta+eq.hvec)/eq.w);

  if (alg.nobasic == 1)
    eq.zb(:) = 0.0;
  end

  % innovation intensities
  eq.xa = maxint(eq.za,p.asigma,p.agamma);
  eq.cxb = maxint(eq.zb,p.bsigma,p.bgamma);

  % probability of positive basic research
  cb = p.bsigma*eq.cxb.^(p.bgamma);
  eq.dstar = eq.cxb.*(1.0+p.rhob).*(eq.bbeta+eq.hvec)/(eq.w*p.bfcoeff) - cb;
  eq.dstart = (log(max(0.0,eq.dstar))-p.dmu)/p.dsig;
  eq.epb = normcdf(eq.dstart);
  eq.xb = eq.epb.*eq.cxb;
  eq.efb = p.bfcoeff*exp(p.dmu+0.5*p.dsig^2.0)*normcdf(eq.dstart-p.dsig);

  % limit things for stability and existence
  if (alg.rlim == 1)
    gsum = (1.0+p.rhoa).*eq.xa + (1.0+p.rhob).*eq.xb;
    glim = eq.tau + alg.glimfact*(eq.r-eq.g);

    for m=1:p.M
      % we need to limit x so that H exists - for now just renormalize
      eq.lim(m) = 0;

      if (glim <= 0.0)
        eq.lim(m) = 1;
        eq.xa(m) = 0.0;
        eq.xb(m) = 0.0;
      elseif (gsum(m) > glim)
        eq.lim(m) = 1;
        eq.xa(m) = eq.xa(m)*(glim/gsum(m));
        eq.xb(m) = eq.xb(m)*(glim/gsum(m));
      end
    end
  end

  % output logs
  for m=1:p.M
    fprintf(alg.fid,'%2i: a = %8.5f, b = %8.5f, epb = %f\n',m,eq.xa(m),eq.xb(m),eq.epb(m));
  end
  fprintf(alg.fid,'eq.lim = ');
  fprintf(alg.fid,'%i ',eq.lim);
  fprintf(alg.fid,'\n');
end

function xval = maxint(z,sig,gam)
  xval = (z/(sig*gam)).^(1.0/(gam-1.0));
end
