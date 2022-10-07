% generate expected production value of new innovations
function findbetas
  global p alg eq qd

  eq.pibar = 1.0/p.epsn;
  eq.denom2 = eq.r+eq.tau+p.kappa-eq.g;
  eq.denom1 = eq.r+eq.tau+p.kappa+eq.g*(p.epsn-2.0);

  eq.meanq = qd.even_pmf*qd.binmids';

  eq.eqhat_a = qd.step_a_pmf*qd.binmids_eps;
  eq.eqhat_b = qd.step_b_pmf*qd.binmids_eps;
  eq.eqhat_e = qd.even_pmf*qd.binmids_eps;

  eq.abeta = eq.pibar*((eq.eqhat_a-p.qminp)/eq.denom1+p.qminp/eq.denom2);
  eq.bbeta = eq.pibar*((eq.eqhat_b-p.qminp)/eq.denom1+p.qminp/eq.denom2);
  eq.ebeta = eq.pibar*((eq.eqhat_e-p.qminp)/eq.denom1+p.qminp/eq.denom2);

  if (alg.fid ~= -1)
    fprintf(alg.fid,'mean_q = %7.5g\n',eq.meanq);
    fprintf(alg.fid,'eqhat_a = %7.5g, eqhat_b = %7.5g, eqhat_e = %7.5g\n',eq.eqhat_a,eq.eqhat_b,eq.eqhat_e);
    fprintf(alg.fid,'abeta = %10.8f, bbeta = %10.8f, ebeta = %10.8f\n',eq.abeta,eq.bbeta,eq.ebeta);
  end
end
