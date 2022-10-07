function qdistbin
  global alg p eq qd

  alpha_a = 1.0/p.alpha;
  alpha_b = 1.0/p.eta;

  qd.nBpow = 8;
  qd.nB = 2^qd.nBpow;
  nbins = qd.nB;

  %qd.qmax = 0.30;
  qd.qmax = 3.0;
  qd.delt = 0.001;

  qd.delq1 = (qd.qmax-p.qmin)/(nbins-1)-1e-12;
  qd.binmins = [0.0 p.qmin:qd.delq1:(qd.qmax+qd.delq1)];
  qd.binmins(nbins+1) = qd.qmax;
  qd.delq = qd.binmins(2:(nbins+1))-qd.binmins(1:nbins);
  qd.binmids = [p.qmin 0.5*(qd.binmins(2:nbins)+qd.binmins(3:(nbins+1)))];

  % calculate the jump sizes and transition probabilities
  mdown = 1.0/(1.0+eq.g*qd.delt);

  % bumat and bdmat are the basis matrices for up and down transitions,
  % respectively. they depend only on mup and mdown. the probabilites
  % are applied to them after construction
  bumat_a = zeros(nbins,nbins);
  bumat_b = zeros(nbins,nbins);
  bdmat = zeros(nbins,nbins);

  % here we assume that the mass in a bin is distributed uniformly over
  % that bin. this way, the transition probabilities change continuously
  % with jump size. bins 1:nbins are active and bin nbins+1 is inactive.

  % DROP-DOWNS
  bdmat(1,1) = 1.0;
  for i=2:nbins
    % calcuate the upper and lower envelope of the tranlated mass in real numbers
    % and bin numbers
    btlow = qd.binmins(i)*mdown;
    bthigh = qd.binmins(i+1)*mdown;
    btlowi = icut(qd.binmins,btlow,nbins);
    bthighi = icut(qd.binmins,bthigh,nbins);

    bdmat(i,btlowi) = bdmat(i,btlowi) + (qd.binmins(btlowi+1)-btlow)/(qd.delq(i)*mdown);
    if (bthighi-btlowi > 1)
      for j=btlowi+1:bthighi-1
        bdmat(i,j) = bdmat(i,j) + qd.delq(j)/(qd.delq(i)*mdown);
      end
    end
    bdmat(i,bthighi) = bdmat(i,bthighi) + (bthigh-qd.binmins(bthighi))/(qd.delq(i)*mdown);
  end

  % MOVE-UPS APP
  bumat_a(1,2:nbins) = exp_cdf(qd.binmins(3:(nbins+1)),p.qmin,alpha_a)-exp_cdf(qd.binmins(2:nbins),p.qmin,alpha_a);
  bumat_a(1,nbins) = bumat_a(1,nbins) + 1.0-exp_cdf(qd.qmax,p.qmin,alpha_a);
  for i=2:nbins-1
    binlow = qd.binmins(i);
    binhigh = qd.binmins(i+1);

    bumat_a(i,i) = 1.0 - inv_exp_int(binhigh,binlow,binhigh,alpha_a)/qd.delq(i);
    bumat_a(i,(i+1):nbins) = (inv_exp_int(qd.binmins((i+1):nbins),binlow,binhigh,alpha_a)-inv_exp_int(qd.binmins((i+2):(nbins+1)),binlow,binhigh,alpha_a))/qd.delq(i);
    bumat_a(i,nbins) = bumat_a(i,nbins) + inv_exp_int(qd.qmax,binlow,binhigh,alpha_a)/qd.delq(i);
  end
  bumat_a(nbins,nbins) = 1.0;

  % MOVE-UPS BAS
  bumat_b(1,2:nbins) = exp_cdf(qd.binmins(3:(nbins+1)),p.qmin,alpha_b)-exp_cdf(qd.binmins(2:nbins),p.qmin,alpha_b);
  bumat_b(1,nbins) = bumat_b(1,nbins) + 1.0-exp_cdf(qd.qmax,p.qmin,alpha_b);
  for i=2:nbins-1
    binlow = qd.binmins(i);
    binhigh = qd.binmins(i+1);

    bumat_b(i,i) = 1.0 - inv_exp_int(binhigh,binlow,binhigh,alpha_b)/qd.delq(i);
    bumat_b(i,(i+1):nbins) = (inv_exp_int(qd.binmins((i+1):nbins),binlow,binhigh,alpha_b)-inv_exp_int(qd.binmins((i+2):(nbins+1)),binlow,binhigh,alpha_b))/qd.delq(i);
    bumat_b(i,nbins) = bumat_b(i,nbins) + inv_exp_int(qd.qmax,binlow,binhigh,alpha_b)/qd.delq(i);
  end
  bumat_b(nbins,nbins) = 1.0;

  %max(abs(sum(bdmat,2)-1.0))
  %max(abs(sum(bumat_a,2)-1.0))
  %max(abs(sum(bumat_b,2)-1.0))

  tmatsize = 2*nbins;
  tmat = zeros(tmatsize,tmatsize);

  pupa = eq.taua*qd.delt;
  pupb = (1-alg.bloss)*eq.taub*qd.delt;
  pspill = eq.taus*qd.delt;
  pcool = p.zeta*qd.delt;
  pdown = 1.0-(eq.taua+eq.taub+eq.taus+p.zeta)*qd.delt;

  % hot -> hot
  tmat(1:nbins,1:nbins) = tmat(1:nbins,1:nbins) + pupa*bumat_b + pupb*bumat_b + (pspill+pdown)*bdmat;
  % hot -> cold
  tmat(1:nbins,(nbins+1):tmatsize) = tmat(1:nbins,(nbins+1):tmatsize) + pcool*bdmat;
  % cold -> cold
  tmat((nbins+1):tmatsize,(nbins+1):tmatsize) = tmat((nbins+1):tmatsize,(nbins+1):tmatsize) + pupa*bumat_a + (pcool+pdown)*bdmat;
  % cold -> hot
  tmat((nbins+1):tmatsize,1:nbins) = tmat((nbins+1):tmatsize,1:nbins) + pupb*bumat_b + pspill*bdmat;

  %max(abs(sum(tmat,2)-1.0))

  qd.bumat_a = bumat_a;
  qd.bumat_b = bumat_b;
  qd.bdmat = bdmat;
  qd.tmat = tmat;

  % normalize all rows, just in case
  for i=1:tmatsize
    tmat(i,:) = tmat(i,:)/sum(tmat(i,:));
  end

  % transpose
  tmat = tmat'-eye(tmatsize);

  % null space of A-I
  tstart = clock;
  bestevec = null(tmat);
  tend = clock;

  bestevec = bestevec/sum(bestevec);

  qd.hot_pmf = bestevec(1:nbins)';
  qd.cold_pmf = bestevec((nbins+1):tmatsize)';
  qd.even_pmf = qd.hot_pmf + qd.cold_pmf;

  qd.step_a_pmf = qd.hot_pmf*bumat_b + qd.cold_pmf*bumat_a;
  qd.step_b_pmf = qd.hot_pmf*bumat_b + qd.cold_pmf*bumat_b;
  qd.pure_a_pmf = qd.hot_pmf*bumat_a + qd.cold_pmf*bumat_a;

  eq.F0 = qd.even_pmf(1);
  eq.hot = sum(qd.hot_pmf);

  qd.even_cmf = cumsum(qd.even_pmf);
  qd.step_a_cmf = cumsum(qd.step_a_pmf);
  qd.step_b_cmf = cumsum(qd.step_b_pmf);

  % aggregates
  qd.binmids_eps = qd.binmids.^(p.epsn-1.0)';
  qd.Qhat_hot = qd.hot_pmf*qd.binmids_eps;
  qd.Qhat_cold = qd.cold_pmf*qd.binmids_eps;
  qd.Qhat = qd.Qhat_hot + qd.Qhat_cold;

  if (alg.check == 1)
    writetable(table({'nbins';'qmax';'delt'},[nbins;qd.qmax;qd.delt],'VariableNames',{'name','value'}),'logs/qdist_algs.csv');
    writetable(table({'alpha';'eta';'qmin';'zeta';'epsn'},[p.alpha;p.eta;p.qmin;p.zeta;p.epsn],'VariableNames',{'name','value'}),'logs/qdist_pars.csv');
    writetable(table({'taua';'taub';'taus';'g'},[eq.taua;eq.taub;eq.taus;eq.g],'VariableNames',{'name','value'}),'logs/qdist_vars.csv');
    writetable(table(qd.binmids',qd.even_pmf',qd.step_a_pmf',qd.step_b_pmf','VariableNames',{'binmids','even_pmf','step_a_pmf','step_b_pmf'}),'logs/qdist.csv');
  end

  if (alg.fid ~= -1)
    if (any(imag(bestevec) ~= 0.0))
      fprintf(alg.fid,'eigenvector imaginary! ');
      fprintf(alg.fid,'imag mean = %f\n',mean(abs(imag(bestevec))));
    end
    fprintf(alg.fid,'eigtime1 = %f\n',tend(6)-tstart(6));
    fprintf(alg.fid,'eq.F0 = %7.5f\n',eq.F0);
    fprintf(alg.fid,'eq.hot = %f\n',eq.hot);
    fprintf(alg.fid,'qhat_hot = %f, qhat_cold = %f, qhat = %f\n',qd.Qhat_hot,qd.Qhat_cold,qd.Qhat);
  end

end

% map from qhat value to bin number - returns largest i1 such that qlist(i1) < qcut
function iout = icut(qlist,qcut,qlen)
  i1 = 1;
  while ((qlist(i1) < qcut) && (i1 <= qlen))
    i1 = i1 + 1;
  end

  iout = i1-1;
end

function rval = inv_exp_int(z,a,b,alpha)
  rval = 1.0/alpha*(exp(-alpha*(z-b))-exp(-alpha*(z-a)));
end

function rval = exp_cdf(z,a,alpha)
  rval = 1.0-exp(-alpha*(z-a));
end

function rval = inv_exp_int2(z,a,b,alpha)
  c = alpha*z;
  rval = exp(alpha)*((c.*expint2(-c/b)+b*exp(-c/b))-(c.*expint2(-c/a)+a*exp(-c/a)));
end

function rval = exp_cdf2(z,a,alpha)
  rval = 1.0-exp(-alpha*(z./a-1.0));
end

function e = expint2(x)
  e = -real(expint(-x));
end
