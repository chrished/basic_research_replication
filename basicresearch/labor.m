function labor
  global alg p eq qd

  % R&D costs
  eq.ca = xcost(eq.xa,p.asigma_ns,p.agamma);
  eq.ccb = xcost(eq.cxb,p.bsigma_ns,p.bgamma);
  eq.cb = eq.epb.*eq.ccb + eq.efb/p.bfcoeff;
  eq.cout = xcost(eq.xout,p.asigma_ns,p.agamma);
  eq.cz = p.massac*(p.acfixed+xcost(eq.xu,p.usigma,p.ugamma)); % always pay academic fixed

  % labor aggregates
  eq.cent = p.massout*eq.cout;
  eq.capp = sum(eq.ca.*eq.phivec);
  eq.cbas = sum(eq.cb.*eq.phivec);
  eq.cinc = eq.capp + eq.cbas;

  eq.prodlab = ((p.epsn-1.0)/p.epsn)/eq.w;
  %eq.prodlab = ((p.epsn-1.0)/p.epsn)*(qd.Qhat/eq.w);
  eq.rndlab = eq.cent + eq.cz + eq.cinc;

  %{
  for m=1:p.M
    fprintf(alg.fid,'%2i: ca = %8.5f, cb = %8.5f\n',m,eq.ca(m),eq.cb(m));
  end
  %}
  fprintf(alg.fid,'eq.cinc = %8.5f\n',eq.cinc);
  fprintf(alg.fid,'eq.cent = %8.5f\n',eq.cent);
  fprintf(alg.fid,'eq.cz = %8.5f\n',eq.cz);

  fprintf(alg.fid,'prodlab = %8.5f, rndlab = %8.5f\n',eq.prodlab,eq.rndlab);
end

function c = xcost(x,sig,gam)
  c = sig*x.^gam;
end
