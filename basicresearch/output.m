function output
  global alg p eq qd

  out_fid = fopen('output/tables.txt','w');

  % table stats
  if (isfile('logs/optvals.txt'))
    optv = load('logs/optvals.txt');
    alg.fbw = optv(1);
    alg.fbg = optv(2);
    crat = eq.w/alg.fbw;
    ceval = (alg.fbw/eq.w)*((p.disc+eq.g*(p.crra-1.0))/(p.disc+alg.fbg*(p.crra-1.0)))^(1.0/(p.crra-1.0));
  else
    crat = nan;
    ceval = nan;
  end

  % policy
  fprintf(out_fid,'Policy:\n');
  fprintf(out_fid,'asubs = %4f\n',p.asubs);
  fprintf(out_fid,'bsubs = %4f\n',p.bsubs);
  fprintf(out_fid,'That = %4f\n',p.That);
  fprintf(out_fid,'\n');

  % after tax wage share
  alab = sum(eq.ca.*eq.phivec) + eq.cent;
  blab = sum(eq.cb.*eq.phivec);
  taxes = p.asubs*alab + p.bsubs*blab;
  wp = eq.w*(1.0-taxes);
  rhobar = sum(p.rhob.*eq.psivec);
  acbar = p.massac*(1.0+p.pb)*eq.xu;

  % labor breakdown
  ca_inc = sum(eq.ca.*eq.phivec);
  cb_inc = sum(eq.cb.*eq.phivec);

  fprintf(out_fid,'Equilibrium variables:\n');
  fprintf(out_fid,'abar = %11.8f\n',eq.abar*100);
  fprintf(out_fid,'atilde = %11.8f\n',eq.a0*100);
  fprintf(out_fid,'dbar = %11.8f\n',acbar*100);
  fprintf(out_fid,'bbar = %11.8f\n',eq.bbar*100);
  fprintf(out_fid,'Lprod = %11.8f\n',eq.prodlab*100);
  fprintf(out_fid,'Lapp_incumbent = %11.8f\n',ca_inc*100);
  fprintf(out_fid,'Lapp_entrant = %11.8f\n',eq.cent*100);
  fprintf(out_fid,'Lbas = %11.8f\n',cb_inc*100);
  fprintf(out_fid,'Lac = %11.8f\n',eq.cz*100);
  fprintf(out_fid,'psi = %11.8f\n',eq.hot*100);
  fprintf(out_fid,'c-ratio = %11.8f\n',crat*100);
  fprintf(out_fid,'g = %11.8f\n',eq.g*100);
  fprintf(out_fid,'alpha = %11.8f\n',ceval*100);
  fprintf(out_fid,'\n');

  % growth decomp
  xout = p.massout*eq.xout;
  xinc = sum(eq.xa.*eq.phivec);
  step_a = (qd.step_a_pmf-qd.even_pmf)*qd.binmids_eps/(p.epsn-1.0);
  step_b = (qd.step_b_pmf-qd.even_pmf)*qd.binmids_eps/(p.epsn-1.0);

  gbase = eq.g*eq.F0*p.qminp;
  gout = xout*step_a;
  gapp = xinc*step_a;
  gbas = eq.taub*step_b;
  g = gbase+gout+gapp+gbas;

  fprintf(out_fid,'Growth Decoposition 1:\n');
  fprintf(out_fid,'gout = %8.6f%% (%8.6f%%)\n',100*gout,100*gout/g);
  fprintf(out_fid,'gapp = %8.6f%% (%8.6f%%)\n',100*gapp,100*gapp/g);
  fprintf(out_fid,'gbas = %8.6f%% (%8.6f%%)\n',100*gbas,100*gbas/g);
  fprintf(out_fid,'gbase = %8.6f%% (%8.6f%%)\n',100*gbase,100*gbase/g);
  fprintf(out_fid,'\n');

  b0 = sum(eq.xb.*eq.phivec);
  bs = sum(p.rhob.*eq.xb.*eq.phivec);
  pure_a = (qd.pure_a_pmf-qd.even_pmf)*qd.binmids_eps/(p.epsn-1.0);
  step_diff = (qd.step_a_pmf-qd.pure_a_pmf)*qd.binmids_eps/(p.epsn-1.0);

  gbase = eq.g*eq.F0*p.qminp;
  ga = eq.taua*pure_a;
  gb = b0*step_b;
  ghot = eq.taua*step_diff;
  gspill = bs*step_b;

  fprintf(out_fid,'Growth Decomposition 2:\n');
  fprintf(out_fid,'ga0 = %8.6f%% (%8.6f%%)\n',100*ga,100*ga/g);
  fprintf(out_fid,'gb0 = %8.6f%% (%8.6f%%)\n',100*gb,100*gb/g);
  fprintf(out_fid,'gbase = %8.6f%% (%8.6f%%)\n',100*gbase,100*gbase/g);
  fprintf(out_fid,'ghot = %8.6f%% (%8.6f%%)\n',100*ghot,100*ghot/g);
  fprintf(out_fid,'gspill = %8.6f%% (%8.6f%%)\n',100*gspill,100*gspill/g);
  fprintf(out_fid,'\n');

  fprintf(out_fid,'Labor Decomposition:\n');
  fprintf(out_fid,'production = %8.6f%%\n',100*eq.prodlab);
  fprintf(out_fid,'applied incumbent = %8.6f%%\n',100*ca_inc);
  fprintf(out_fid,'basic incumbent = %8.6f%%\n',100*cb_inc);
  fprintf(out_fid,'entrants = %8.6f%%\n',100*eq.cent);
  fprintf(out_fid,'academic = %8.6f%%\n',100*eq.cz);
  fprintf(out_fid,'\n');

  fclose(out_fid);

  %%
  %% citations
  %%

  nC = 50;
  nK = 10;

  % grids
  cvec = 0:nC;
  kvec = 1:nK;
  mvec = 1:p.M;

  % this is the "natural" case
  % relah = eq.tau*x*p.alpha/p.zeta;
  % relb = eq.tau*x*p.eta/p.zeta;
  % decah = relah/(1.0+relah);
  % decb = relb/(1.0+relb);
  % citah = (1.0-decah)*decah.^cvec;
  % cita = [(eq.hot*citah(1)+(1.0-eq.hot)) citah(2:end)];
  % citb = (1.0-decb)*decb.^cvec;

  % spil - M x K
  % citk - K x C
  % cite - M x C

  % this is the "proportional" case
  exph = eq.taub;
  expc = eq.taub;
  inoh = eq.taua*p.citx*p.eta;
  inoc = eq.taua*p.citx*p.alpha;
  probh = inoh/(exph+inoh);
  probc = inoc/(expc+inoc);

  % citations (c) conditional on industry presence (m)
  cith = (1-probh)*probh.^cvec;
  citc = (1-probc)*probc.^cvec;

  % aggregate over hot/hold
  citb = cith;
  cita = eq.hot*cith + (1.0-eq.hot)*citc;

  % aggregate over public private
  citpub = citb;
  citpri = (eq.taua*cita+eq.taub*citb)/eq.tau;

  % save results
  writetable(array2table([cvec;citpri;citpub]','VariableNames',{'cites','private','public'}),'output/citations.csv');

  %%
  %% epb/bri for graphs
  %%

  epb = eq.epb;
  bri = eq.cb./eq.ca;

  cbri = zeros(1,8);
  cbri(1:7) = bri(1:7);
  cbri(8) = sum(eq.psivec(8:end).*bri(8:end))/sum(eq.psivec(8:end));

  cepb = zeros(1,8);
  cepb(1:7) = epb(1:7);
  cepb(8) = sum(eq.psivec(8:end).*epb(8:end))/sum(eq.psivec(8:end));

  writetable(array2table([1:8;cepb;cbri]','VariableNames',{'m','epb','bri'}),'output/bristats.csv');
end
