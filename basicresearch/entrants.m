function entrants
  global alg p eq

  hgain = [(eq.hvec(2:end)-eq.hvec(1)) 0.0];
  eq.vmerge = p.merge*sum(p.eprob.*hgain.*eq.psipdf);
  eq.zout = max(0.0,(eq.abeta+eq.hvec(1)+eq.vmerge)/eq.w);
  eq.xout = maxint(eq.zout,p.asigma,p.agamma);

  fprintf(alg.fid,'eq.zout = %8.5f\n',eq.zout);
  fprintf(alg.fid,'eq.xout = %8.5f\n',eq.xout);
end

function xval = maxint(z,sig,gam)
  xval = (z/(sig*gam)).^(1.0/(gam-1.0));
end
