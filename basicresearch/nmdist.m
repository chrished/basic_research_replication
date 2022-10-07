function nmdist
  global p alg eq

  % we want to set this so pstay for (alg.prodmax,p.M) is ~0.5
  alg.delt = 0.5/(alg.prodmax*((1.0+p.pa)*eq.xa(p.M)+(1.0+p.pb)*eq.xb(p.M)+p.kappa+eq.tau)+eq.mrate+p.kappa);

  % transition probabilities
  pdown = alg.delt*(eq.tau);
  mup0  = alg.delt*(eq.mrate*(1.0-p.eprob));
  mup1  = alg.delt*(eq.mrate*p.eprob);
  pout  = alg.delt*(p.kappa);

  % spillovers
  Kmax = 10;
  nK = Kmax + 1;
  mvec = 1:p.M;
  kvec = (0:Kmax)';

  sma = p.pa*mvec./(p.pa*mvec+p.M);
  smb = p.pb*mvec./(p.pb*mvec+p.M);
  pta = 1 - sma.^(Kmax+1);
  ptb = 1 - smb.^(Kmax+1);
  pka = (1-sma).*sma.^kvec./pta;
  pkb = (1-smb).*smb.^kvec./ptb;
  eq.pupa = pka; % shifted up one
  eq.pupb = pkb; % shifted up one

  pup = alg.delt*(eq.xa.*eq.pupa+eq.xb.*eq.pupb);
  pup(1,:) = pup(1,:) + alg.delt*(p.kappa+eq.freebase); % free, no spillover

  % cuda solver
  eigv = eigsim(pup,mup0,mup1,pdown,pout,alg.prodmax,nK,0);
  eq.eigm = max(0.0,reshape(eigv/sum(eigv),alg.prodmax,p.M));
  eq.nmass = (1:alg.prodmax)'.*eq.eigm;

  % normalize to mass one of products
  eq.phivec = sum(eq.nmass,1);
  prodtot = sum(eq.phivec);
  eq.eigm = eq.eigm/prodtot;
  eq.nmass = eq.nmass/prodtot;
  eq.phivec = eq.phivec/prodtot;

  % construct psivec
  eq.psivec = sum(eq.eigm,1);
  eq.fsize = sum(eq.psivec);
  eq.psipdf = eq.psivec/eq.fsize;

  % get nsum(1)
  eq.nsum1 = sum(eq.eigm(1,:));
  eq.n1vec = eq.eigm(1,:)./eq.phivec;

  % check imputed fsize
  fsize2 = eq.tau*eq.nsum1/(eq.mrate*(1.0-p.merge)/p.merge-p.kappa);
  fprintf(alg.fid,'fsize = %8.5f, fsize2 = %8.5f\n',eq.fsize,fsize2);
  %{
  for m=1:p.M
    fprintf(alg.fid,'psipdf(%i) = %8.5f\n',m,eq.psipdf(m));
  end
  %}

  % calculate firms at prodmax
  eq.firmpdf = eq.eigm/eq.fsize;
  eq.npdf = sum(eq.firmpdf,2);
  eq.t1firms = 1.0 - eq.npdf(end);

  fprintf(alg.fid,'npdf(end) = %15.10f\n',eq.npdf(end));

  if (alg.check == 1)
    writetable(table({'pdown';'pout';'M';'prodmax'},[pdown;pout;p.M;alg.prodmax],'VariableNames',{'name','value'}),'logs/nmdist_scalar.csv');
    writetable(table((1:p.M)',mup0',mup1','VariableNames',{'m','mup0','mup1'}),'logs/nmdist_vector.csv');
    writetable(array2table([(1:alg.prodmax)',eq.eigm],'VariableNames',{'n','m1','m2','m3','m4','m5','m6','m7','m8','m9','m10'}),'logs/nmdist_eigm.csv');
  end

  %{
  if (alg.check == 1)
    % check the flow equations
    nmcond = zeros(alg.prodmax,p.M);

    eq.xbtot = eq.xb;
    eq.x1tot = eq.xbtot.*(1.0-p.rho)+eq.xa;
    eq.x2tot = eq.xbtot.*p.rho;
    eq.xtot = eq.xbtot+eq.xa;

    % m = 1
    nmcond(1,1) = (eq.tau+eq.xtot(1)+p.kappa+eq.xe(1)+p.kappa)*eq.eigm(1,1) - 2*eq.tau*eq.eigm(2,1) - p.massout*eq.xout;
    nmcond(2,1) = (2*(eq.tau+eq.xtot(1)+p.kappa)+eq.xe(1)+p.kappa)*eq.eigm(2,1) - 3*eq.tau*eq.eigm(3,1) - (eq.x1tot(1)+p.kappa)*eq.eigm(1,1);
    for n=3:(alg.prodmax-2)
      nmcond(n,1) = (n*(eq.tau+eq.xtot(1)+p.kappa)+eq.xe(1)+p.kappa)*eq.eigm(n,1) - (n-2)*(eq.x2tot(1))*eq.eigm(n-2,1) - (n-1)*(eq.x1tot(1)+p.kappa)*eq.eigm(n-1,1) - (n+1)*eq.tau*eq.eigm(n+1,1);
    end
    nmcond(alg.prodmax-1,1) = ((alg.prodmax-1)*(eq.tau+eq.x1tot(1)+p.kappa)+eq.xe(1)+p.kappa)*eq.eigm(alg.prodmax-1,1) - (alg.prodmax-3)*(eq.x2tot(1))*eq.eigm(alg.prodmax-3,1) - (alg.prodmax-2)*(eq.x1tot(1)+p.kappa)*eq.eigm(alg.prodmax-2,1) - (alg.prodmax)*eq.tau*eq.eigm(alg.prodmax,1);
    nmcond(alg.prodmax,1) = (alg.prodmax*eq.tau+eq.xe(1)+p.kappa)*eq.eigm(alg.prodmax,1) - (alg.prodmax-1)*(eq.x1tot(1)+p.kappa)*eq.eigm(alg.prodmax-1,1) - (alg.prodmax-2)*(eq.x2tot(1))*eq.eigm(alg.prodmax-2,1);

    % m > 1
    for m=2:p.M
      nmcond(1,m) = (eq.tau+eq.xtot(m)+p.kappa+eq.xe(m)+p.kappa)*eq.eigm(1,m) - 2*eq.tau*eq.eigm(2,m) - eq.xe(m-1)*eq.eigm(1,m-1);
      nmcond(2,m) = (2*(eq.tau+eq.xtot(m)+p.kappa)+eq.xe(m)+p.kappa)*eq.eigm(2,m) - 3*eq.tau*eq.eigm(3,m) - (eq.x1tot(m)+p.kappa)*eq.eigm(1,m) - eq.xe(m-1)*eq.eigm(2,m-1);
      for n=3:(alg.prodmax-2)
        nmcond(n,m) = (n*(eq.tau+eq.xtot(m)+p.kappa)+eq.xe(m)+p.kappa)*eq.eigm(n,m) - (n-2)*(eq.x2tot(m))*eq.eigm(n-2,m) - (n-1)*(eq.x1tot(m)+p.kappa)*eq.eigm(n-1,m) - (n+1)*eq.tau*eq.eigm(n+1,m) - eq.xe(m-1)*eq.eigm(n,m-1);
      end
      nmcond(alg.prodmax-1,m) = ((alg.prodmax-1)*(eq.tau+eq.x1tot(m)+p.kappa)+eq.xe(m)+p.kappa)*eq.eigm(alg.prodmax-1,m) - (alg.prodmax-3)*(eq.x2tot(m))*eq.eigm(alg.prodmax-3,m) - (alg.prodmax-2)*(eq.x1tot(m)+p.kappa)*eq.eigm(alg.prodmax-2,m) - (alg.prodmax)*eq.tau*eq.eigm(alg.prodmax,m) - eq.xe(m-1)*eq.eigm(alg.prodmax-1,m-1);
      nmcond(alg.prodmax,m) = (alg.prodmax*eq.tau+eq.xe(m)+p.kappa)*eq.eigm(alg.prodmax,m) - (alg.prodmax-1)*(eq.x1tot(m)+p.kappa)*eq.eigm(alg.prodmax-1,m) - (alg.prodmax-2)*(eq.x2tot(m))*eq.eigm(alg.prodmax-2,m) - eq.xe(m-1)*eq.eigm(alg.prodmax,m-1);
    end

    disp(['sscond = ' num2str(max(max(abs(nmcond))))]);
    disp(['sum(phivec) = ' num2str(sum(eq.phivec),8)]);
  end
  %}
end
