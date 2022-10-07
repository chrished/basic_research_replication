function growthrate
  global alg p eq qd

  gbase = eq.g*eq.F0*p.qminp;
  ga = eq.taua*(qd.step_a_pmf-qd.even_pmf)*qd.binmids_eps/(p.epsn-1.0);
  gb = eq.taub*(qd.step_b_pmf-qd.even_pmf)*qd.binmids_eps/(p.epsn-1.0);
  eq.newg = (ga+gb)/(1.0-eq.F0*p.qminp);

  if (alg.fid ~= -1)
    fprintf(alg.fid,'ga = %10.8f, gb = %10.8f, gbase = %10.8f, g = %10.8f\n',ga,gb,gbase,eq.newg);
  end

end
