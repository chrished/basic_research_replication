function socplan
  global p walg alg

  alg = {};
  initalg();
  alg.fid = -1;

  p = {};
  [par0,names] = parse_params(alg.par_file);
  loadparams(par0);

  cev = load('logs/cevals.txt');
  walg.lastg = 0.02;
  walg.ce = 0;
  walg.cebw = cev(1);
  walg.cebg = cev(2);

  %abh_init = [0.09085315 0.03530928 -2.47792617]; % 0
  %abh_init = [0.06255120 0.02906347 -4.07219379]; % 15
  %abh_init = [0.06958428 0.02745587 -4.16500972]; % 16
  %abh_init = [0.07470080 0.02533399 -4.28694370]; % 17
  %abh_init = [0.09094270 0.03309137 -2.54189101]; % 18
  %abh_init = [0.07799962 0.02874590 -2.96816904]; % 20
  %abh_init = [0.08115361 0.03435960 -2.67961879]; % 21
  %abh_init = [0.08482022 0.03513487 -2.64194748]; % 22
  %abh_init = [0.08481649 0.03513395 -2.63823266]; % 23
  %abh_init = [0.08312293 0.03508161 -2.73321923]; % 25
  %abh_init = [0.07156417 0.02992827 -3.10807492]; % 37
  %abh_init = [0.11035934 0.03247037 -2.51359504]; % 38
  %abh_init = [0.10319643 0.03243764 -2.72791469]; % 39
  %abh_init = [0.10548388 0.03162096 -2.86694974]; % 40
  %abh_init = [0.10561482 0.03164931 -3.06093742]; % 41
  %abh_init = [0.12462828 0.03278447 -4.49844687]; % 45
  %abh_init = [0.12195949 0.01180320 -6.82260061] % 47
  %abh_init = [0.10666941 0.03122819 -4.45317314] % 49
  %abh_init = [0.12090844 0.03978253 -3.14817543] % 62
  abh_init = [0.12323884 0.03910684 -3.15840373] % 63

  opts = optimset('Display','iter','TolX',1e-4);
  [maxabh,mv] = fminsearch(@wobj,abh_init,opts);
  %maxabh = abh_init;

  fprintf(1,'abh_init = [%10.8f %10.8f %10.8f]\n',maxabh(1),maxabh(2),maxabh(3));

  walg.ce = 1;
  wobj(maxabh);

  % save these for the consumption equiv calculations
  optw = walg.optw;
  optg = walg.optg;
  save('logs/optvals.txt','optw','optg','-ascii','-double');
end

% a - applied intensity
% b - basic intensity
% hf - log of basic fixed cost cutoff
function wv = wobj(abh)
  global p eq qd walg

  a = abh(1);
  b = abh(2);
  hf = abh(3);

  dstar = (hf-p.dmu)/p.dsig;
  Pf = normcdf(dstar);
  cf = exp(p.dmu+0.5*p.dsig^2.0)*normcdf(dstar-p.dsig);
  cu = exp(p.dmu+0.5*p.dsig^2.0)*p.massac;

  abar = (1.0+p.pa)*(1.0+p.massout)*a;
  bbar = (1.0+p.pb)*(p.massac+Pf)*b;

  fprintf(1,'a = %f, b = %f, hf = %f, Pf = %f\n',a,b,hf,Pf);

  % solve for F0/g
  function gdiff = gobj(gin)
    eq.g = gin;
    eq.taua = abar;
    eq.taub = bbar;
    eq.taus = 0.0;

    eq.tau = eq.taua+eq.taub;
    eq.r = eq.g*p.crra + p.disc;

    qdistbin();
    growthrate();

    gdiff = gin - eq.newg;

    fprintf(1,'gin = %15.12f, gdiff = %15.12f\n',gin,gdiff);
  end

  opts = optimset('TolX',1e-10,'Display','none');
  %g = fzero(@gobj,[0.005 (walg.lastg+0.02)],opts);
  [g,gdiff,geval] = fsolve(@gobj,walg.lastg,opts);
  walg.lastg = g;

  if (geval <= 0)
    fprintf(1,'Could not find equilibrium g!\n');
  end

  crnd = (1.0+p.massout)*xcost(a,p.asigma_ns,p.agamma) + (p.massac+Pf)*xcost(b,p.bsigma_ns,p.bgamma) + (cf+cu);
  w = ((p.epsn-1.0)/p.epsn)/(1.0-crnd);

  wv = -log(p.disc+g*(p.crra-1.0)) + (p.crra-1.0)*log(w);

  walg.ceval = (walg.cebw/w)*((p.disc+g*(p.crra-1.0))/(p.disc+walg.cebg*(p.crra-1.0)))^(1.0/(p.crra-1.0));

  ceval0 = (walg.cebw/w)*((p.disc+g*(p.crra-1.0))/(p.disc+walg.cebg*(p.crra-1.0)))^(1.0/(p.crra-1.0));
  Zrat0 = w/walg.cebw;

  fprintf(1,'gfin = %13.10f\n',g);
  fprintf(1,'rndlab = %13.10f\n',crnd);
  fprintf(1,'wv = %13.10f\n',wv);
  fprintf(1,'eqce = %13.10f\n',walg.ceval);

  if (walg.ce == 1)
    %{}
    sp_fid = fopen('output/socplan.txt','w');
    fprintf(sp_fid,'abar = %11.8f\n',100*a);
    fprintf(sp_fid,'atilde = %11.8f\n',100*p.massout*a);
    fprintf(sp_fid,'dbar = %11.8f\n',100*p.massac*b);
    fprintf(sp_fid,'bbar = %11.8f\n',100*Pf*b);
    fprintf(sp_fid,'Pf = %11.8f\n',100*Pf);
    fprintf(sp_fid,'Lprod = %11.8f\n',100*(1.0-crnd));
    fprintf(sp_fid,'Lapp_incumbent = %11.8f\n',100*xcost(a,p.asigma_ns,p.agamma));
    fprintf(sp_fid,'Lapp_entrant = %11.8f\n',100*p.massout*xcost(a,p.asigma_ns,p.agamma));
    fprintf(sp_fid,'Lbas = %11.8f\n',100*(Pf*xcost(b,p.bsigma_ns,p.bgamma)+cf));
    fprintf(sp_fid,'Lac = %11.8f\n',100*(p.massac*xcost(b,p.bsigma_ns,p.bgamma)+cu));
    fprintf(sp_fid,'psi = %11.8f\n',100*eq.hot);
    fprintf(sp_fid,'c-ratio = %11.8f\n',100/Zrat0);
    fprintf(sp_fid,'g = %11.8f\n',100*g);
    fprintf(sp_fid,'alpha = %11.8f\n',100*ceval0);
    fprintf(sp_fid,'\n');

    fclose(sp_fid);

    %{
    % growth decomposition
    xout = p.massout*a;
    xinc = a;

    gbase = g*eq.F0*(p.qmin^(p.epsn-1.0))/(eq.eqhat1);
    gout = xout*(qd.binmids_eps*qd.step_a_pmf'-qd.binmids_eps*qd.even_pmf')/((p.epsn-1.0)*eq.eqhat1);
    gapp = xinc*(qd.binmids_eps*qd.step_a_pmf'-qd.binmids_eps*qd.even_pmf')/((p.epsn-1.0)*eq.eqhat1);
    gbas = bbar*(qd.binmids_eps*qd.step_b_pmf'-qd.binmids_eps*qd.even_pmf')/((p.epsn-1.0)*eq.eqhat1);
    fprintf(1,'type 1:\n');
    fprintf(1,'gout = %8.6f%% (%8.6f%%)\n',100*gout,100*gout/g);
    fprintf(1,'gapp = %8.6f%% (%8.6f%%)\n',100*gapp,100*gapp/g);
    fprintf(1,'gbas = %8.6f%% (%8.6f%%)\n',100*gbas,100*gbas/g);
    fprintf(1,'gbase = %8.6f%% (%8.6f%%)\n',100*gbase,100*gbase/g);

    b0 = (p.massac+Pf)*b;
    bs = p.p*(p.massac+Pf)*b;

    gbase = g*eq.F0*(p.qmin^(p.epsn-1.0))/(eq.eqhat1);
    ga = abar*(qd.binmids_eps*qd.pure_a_pmf'-qd.binmids_eps*qd.even_pmf')/((p.epsn-1.0)*eq.eqhat1);
    gb = b0*(qd.binmids_eps*qd.step_b_pmf'-qd.binmids_eps*qd.even_pmf')/((p.epsn-1.0)*eq.eqhat1);
    ghot = abar*(qd.binmids_eps*qd.step_a_pmf'-qd.binmids_eps*qd.pure_a_pmf')/((p.epsn-1.0)*eq.eqhat1);
    gspill = bs*(qd.binmids_eps*qd.step_b_pmf'-qd.binmids_eps*qd.even_pmf')/((p.epsn-1.0)*eq.eqhat1);
    fprintf(1,'type 2:\n');
    fprintf(1,'ga0 = %8.6f%% (%8.6f%%)\n',100*ga,100*ga/g);
    fprintf(1,'gb0 = %8.6f%% (%8.6f%%)\n',100*gb,100*gb/g);
    fprintf(1,'gbase = %8.6f%% (%8.6f%%)\n',100*gbase,100*gbase/g);
    fprintf(1,'ghot = %8.6f%% (%8.6f%%)\n',100*ghot,100*ghot/g);
    fprintf(1,'gspill = %8.6f%% (%8.6f%%)\n',100*gspill,100*gspill/g);
    fprintf(1,'\n');

    fprintf(1,'cebw = %11.8f\n',walg.cebw);
    fprintf(1,'cebg = %11.8f\n',walg.cebg);
    fprintf(1,'\n');
    %}

    Zrat = w/walg.cebw;
    fprintf(1,'Z*/ZSP = %11.8f\n',Zrat);

    %{
    walg.ceval = (w/walg.cebw)*((p.disc+walg.cebg*(p.crra-1.0))/(p.disc+g*(p.crra-1.0)))^(1.0/(p.crra-1.0));
    fprintf(1,'eqce = %10.8f\n',walg.ceval);
    %}

    walg.optw = w;
    walg.optg = g;
  end

end

function c = xcost(x,sig,gam)
  c = sig*x.^gam;
end

