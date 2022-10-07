function [wval,ceval,info] = welfare(asubs,bsubs,polzout,polmac,lname)
  global alg p eq

  % alg init
  alg = {};
  initalg();

  % policies
  p = {};
  if (~isnan(asubs))
    p.asubs = asubs;
  end
  if (~isnan(bsubs))
    p.bsubs = bsubs;
  end
  if (~isnan(polzout))
    p.polzout = polzout;
  end
  if (~isnan(polmac))
    p.polmac = polmac;
  end

  % load in stuff
  if (isfile('logs/weqvars.txt'))
    eqv_file = 'logs/weqvars.txt';
  else
    eqv_file = alg.eqv_file;
  end
  [par0,names] = parse_params(alg.par_file);
  eqv0 = load(eqv_file);
  loadparams(par0);

  % solver for w
  options = optimset('Display','iter','MaxFunEvals',250,'TypicalX',eqv0,'DiffMinChange',1e-6);
  [eqfin,eqdiff,eflag] = fsolve(@eqfunc,eqv0,options);
  eqout = eqfunc(eqfin);

  if (eflag <= 0)
    disp('Failed.');
    wval = nan;
    ceval = nan;
    info = 'FAILED';
    return
  end

  err = mean(abs(eqout));
  if (err > 1e-6)
    disp(['Error: err = ' num2str(err)]);
    wval = nan;
    ceval = nan;
    info = 'ERROR';
    return
  end

  optv = load('logs/optvals.txt');
  save('logs/weqvars.txt','eqfin','-ascii','-double');

  alg.fbw = optv(1);
  alg.fbg = optv(2);

  wval = log(p.disc+eq.g*(p.crra-1.0)) - (p.crra-1.0)*log(eq.w);
  ceval = (alg.fbw/eq.w)*((p.disc+eq.g*(p.crra-1.0))/(p.disc+alg.fbg*(p.crra-1.0)))^(1.0/(p.crra-1.0));
  Zrat = eq.w/alg.fbw;

  info = sprintf(' %5.3f   %5.3f   %5.3f   %5.3f   %5.3f   %13.10f   %12.10f   %10.8f   %10.8f   %10.8f\n',p.asubs,p.bsubs,eq.w*eq.cz,p.massac,alg.bdole,wval,ceval,eq.g,eq.rndlab,1.0/Zrat);
  fprintf(1,info);

  if (~isnan(lname))
    unix(['mv logs/eqlog.txt policy/' lname]);
  end

  ca_inc = sum(eq.ca.*eq.phivec);
  ca_tot = ca_inc + eq.cent;
  cb_tot = sum(eq.cb.*eq.phivec);
  taxes = p.asubs*ca_tot + p.bsubs*cb_tot;
  wp = eq.w*(1.0-taxes);
  epb_tot = sum(eq.epb.*eq.psipdf);
  acbar = p.massac*(1.0+p.pb)*eq.xu;

  pol_fid = fopen('policy/summary.txt','w');
  fprintf(pol_fid,'asubs = %11.8f\n',p.asubs*100);
  fprintf(pol_fid,'bsubs = %11.8f\n',p.bsubs*100);
  fprintf(pol_fid,'academic gdp frac = %11.8f\n',p.That*100);
  fprintf(pol_fid,'university labs = %11.8f\n',p.massac*100);
  fprintf(pol_fid,'Bayh-Dole = %11.8f\n',alg.bdole*100);
  fprintf(pol_fid,'abar = %11.8f\n',eq.abar*100);
  fprintf(pol_fid,'atilde = %11.8f\n',eq.a0*100);
  fprintf(pol_fid,'dbar = %11.8f\n',acbar*100);
  fprintf(pol_fid,'bbar = %11.8f\n',eq.bbar*100);
  fprintf(pol_fid,'Lprod = %11.8f\n',eq.prodlab*100);
  fprintf(pol_fid,'Lapp_incumbent = %11.8f\n',ca_inc*100);
  fprintf(pol_fid,'Lapp_entrant = %11.8f\n',eq.cent*100);
  fprintf(pol_fid,'Lbas = %11.8f\n',cb_tot*100);
  fprintf(pol_fid,'Lac = %11.8f\n',eq.cz*100);
  fprintf(pol_fid,'psi = %11.8f\n',eq.hot*100);
  fprintf(pol_fid,'c-ratio = %11.8f\n',1.0/Zrat*100);
  fprintf(pol_fid,'g = %11.8f\n',eq.g*100);
  fprintf(pol_fid,'alpha = %11.8f\n',ceval*100);
  fprintf(pol_fid,'abeta = %11.8f\n',eq.abeta);
  fprintf(pol_fid,'bbeta = %11.8f\n',eq.bbeta);
  fclose(pol_fid);

end
