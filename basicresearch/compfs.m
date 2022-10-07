function [score,m_vec_wl,m_dat_wl] = compfs(do_run)
  global m_vec m_desc m_pos alg p eq qd

  %load('workspace.mat');

  if (nargin == 0)
    do_run = 1;
  end

  % type dependent
  qbins = qd.binmids;
  norm_cdf_a = qd.step_a_cmf/qd.step_a_cmf(end);
  norm_cdf_b = qd.step_b_cmf/qd.step_b_cmf(end);
  norm_cdf_even = qd.even_cmf/qd.even_cmf(end);

  qdists = [norm_cdf_a norm_cdf_b norm_cdf_even];

  % spillovers
  dpupa = cumsum(eq.pupa, 1);
  dpupb = cumsum(eq.pupb, 1);

  % load in
  xat = eq.xa;
  xbt = eq.xb;
  xet = eq.mrate;
  eprt = p.eprob;
  tau = eq.tau;
  g = eq.g;
  r = eq.r;
  kappa = p.kappa;
  epsn = p.epsn;
  qmin = p.qmin;
  nBpow = qd.nBpow;

  % launch 'er
  R_PER_T = 100;

  if (do_run == 1)
    [fnprod,fmind,fage,fmpos,fqpow1,fnprod_zero,fmind_zero,fage_zero,fmpos_zero,fqpow1_zero,fexited] = firmsim(nBpow,qbins,qdists,xat,xbt,xet,eprt,dpupa,dpupb,tau,g,r,epsn,kappa,qmin,int32(R_PER_T));
    % save('firmpanel.mat');
  else
    load('firmpanel.mat');
  end

  dfnprod = double(fnprod)';
  dfmind = double(fmind);
  dfage = double(fage);
  dfmpos = double(fmpos);
  dfqpow1 = double(fqpow1)';
  dfnprod_zero = double(fnprod_zero)';
  dfmind_zero = double(fmind_zero);
  dfage_zero = double(fage_zero);
  dfmpos_zero = double(fmpos_zero);
  dfqpow1_zero = double(fqpow1_zero)';
  lfexited = logical(fexited);

  check = 1;

  nF = length(dfnprod);

  % derivative firm characteristics
  fexited = lfexited;
  fnot_exited = not(lfexited);
  n_exited = sum(fexited);
  n_not_exited = nF-n_exited;

  % age in years
  dfage_year = dfage/R_PER_T;

  % labor coeff
  lcf = (1.0/eq.w)*((p.epsn-1.0)/p.epsn); % = eq.prodlab

  % expected basic total cost
  ctb = eq.ccb + (eq.efb./eq.epb)/p.bfcoeff;
  ctb(eq.epb==0) = 0.0; % in case epb is truly zero

  % basic probabilities
  dfepb = eq.epb(dfmind);
  dfepb_zero = eq.epb(dfmind_zero);

  % initial employment
  dflabor_zero = lcf*dfqpow1_zero;

  dfrnd_empl_a_zero = eq.ca(dfmind_zero).*dfnprod_zero;
  dfrnd_empl_b_zero = ctb(dfmind_zero).*dfepb_zero.*dfnprod_zero;
  dfrnd_empl_tot_zero = dfrnd_empl_a_zero + dfrnd_empl_b_zero;

  dfempl_tot_zero = dflabor_zero + dfrnd_empl_tot_zero;

  % final employment
  dflabor = lcf*dfqpow1;

  dfrnd_empl_a = eq.ca(dfmind).*dfnprod;
  dfrnd_empl_b = ctb(dfmind).*dfepb_zero.*dfnprod;
  dfrnd_empl_tot = dfrnd_empl_a + dfrnd_empl_b;

  dfempl_tot = dflabor + dfrnd_empl_tot;

  % employment change
  dfempl_change = dfempl_tot - dfempl_tot_zero;

  % sales
  dfsales_zero = dfqpow1_zero;
  dfsales = dfqpow1;

  % sales change
  dfsales_change = dfsales - dfsales_zero;

  % profit
  dfprofit_zero = eq.pibar*dfqpow1_zero - eq.w*dfrnd_empl_tot_zero;
  dfros_zero = dfprofit_zero./dfsales_zero;

  % sales growth
  dfsales_growth = dfsales_change(fnot_exited)./dfsales_zero(fnot_exited)+eq.g;

  % employment growth
  dfempl_growth = dfempl_change(fnot_exited)./dfempl_tot_zero(fnot_exited);

  % R&D to sales ratio
  dfrnd_a_to_sales_zero = eq.w*dfrnd_empl_a_zero./dfsales_zero;
  dfrnd_b_to_sales_zero = eq.w*dfrnd_empl_b_zero./dfsales_zero;
  dfrnd_to_sales_zero = eq.w*dfrnd_empl_tot_zero./dfsales_zero;

  % R&D to nprod ratio
  dfrnd_to_nprod = eq.w*dfrnd_empl_tot_zero./dfnprod_zero;

  % industry expansion
  lfexpanded = dfmpos(fnot_exited) > dfmpos_zero(fnot_exited);

  % positive growth
  lfpos_growth = dfsales_growth > 0.0;

  % firm value
  dfbetaq_zero = eq.pibar*((dfqpow1_zero-p.qminp)/eq.denom1+p.qminp/eq.denom2);
  dfval_zero = dfbetaq_zero + dfnprod_zero.*eq.hvec(dfmind_zero);

  dfbetaq = eq.pibar*((dfqpow1-p.qminp)/eq.denom1+p.qminp/eq.denom2);
  dfval = dfbetaq + dfnprod.*eq.hvec(dfmind_zero);

  dfval_change = dfval-dfval_zero;
  dfval_growth = dfval_change(fnot_exited)./dfval_zero(fnot_exited);

  dfval_to_empl = dfval_zero./dfempl_tot_zero;

  % save for later
  save('sim_results.mat');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% calculate moments                                                   %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  wlo = 0.05;
  whi = 0.95;

  % init moments
  m_vec = [];
  m_desc = {};
  m_wgts = [];
  m_pos = 1;

  % sales growth - basic and no basic
  %b_ne_sel = (dfdo_basic(fnot_exited) == 1);
  %nb_ne_sel = ~b_ne_sel;

  %mean_sales_growth_b = mean(dfsales_growth(b_ne_sel));
  %mean_sales_growth_nb = mean(dfsales_growth(nb_ne_sel));
  mean_sales_growth_b = 0.0;
  mean_sales_growth_nb = 0.0;

  add_moment(mean_sales_growth_b,'sales_growth_b');
  add_moment(mean_sales_growth_nb,'sales_growth_nb');

  % spillover regression
  add_moment(0.0,'spill_reg');

  % basic research intensive margin
  fepb = eq.epb(dfmind_zero);
  for m=1:7
    d_epb(m) = mean(fepb(dfmind_zero == m));
    m_string = ['Basic Extensive $m = ' num2str(m) '$'];
    add_moment(d_epb(m),m_string);
  end
  d_epb(8) = mean(fepb(dfmind_zero >= 8));
  m_string = 'Basic Extensive $m >= 8$';
  add_moment(d_epb(8),m_string);

  % mean basic research intensity by m
  rnd_int = dfrnd_empl_b_zero./dfrnd_empl_a_zero;
  for m=1:7
    d_bri(m) = mean(rnd_int(dfmind_zero == m));
    m_string = ['Basic Intensive $m = ' num2str(m) '$'];
    add_moment(d_bri(m),m_string);
  end
  d_bri(8) = mean(rnd_int(dfmind_zero >= 8));
  m_string = 'Basic Intensive $m >= 8$';
  add_moment(d_bri(8),m_string);

  % industry presence
  mean_m = mean(dfmind_zero);
  mean_m2 = mean(dfmind_zero.^2);
  add_moment(mean_m,'Mean Industries');
  add_moment(mean_m2,'Mean Square Industries');

  % profitability
  pos_profit = dfros_zero > 0.0;
  %profit = mean(dfros_zero(pos_profit));
  profit = median(dfros_zero);
  %profit = mean_wins(dfros_zero,wlo,whi);
  add_moment(profit,'Return on Sales');

  % firm exit rate
  exit_rate = mean(fexited);
  add_moment(exit_rate,'Exit Rate');

  % R&D to sales ratio
  rnd_to_sales = mean_wins(dfrnd_to_sales_zero,wlo,whi);
  add_moment(rnd_to_sales,'R\&D/Sales');

  % applied R&D to sales ratio
  %rnd_a_to_sales = mean(dfrnd_empl_a./dfsales);
  %add_moment(rnd_a_to_sales,'rnd_a_to_sales',1.0);

  % basic R&D to sales ratio
  %rnd_b_to_sales = mean(dfrnd_empl_b./dfsales);
  %add_moment(rnd_b_to_sales,'rnd_b_to_sales',1.0);

  % dispersion of employment
  %empl_disp = mean(fempl)/std(dfqpow1);
  %add_moment(empl_disp,'empl_disp',1.0);

  % private to public basic research spending
  %pub_priv_ratio = mean(eq.cb(dfmind_zero))/eq.cz/p.M;
  pub_priv_ratio = sum(eq.cb.*eq.phivec)/eq.cz;
  add_moment(pub_priv_ratio,'Private/public basic');

  %add_moment(sales_growth,'sales_growth',1.0);
  empl_growth = mean_wins(dfempl_growth,wlo,whi);
  add_moment(empl_growth,'Employment Growth');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% rnd/sales breakdown                                                     %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % set up quartiles
  nquarts = 4;
  quart_size = ceil(n_not_exited/nquarts);
  qrange = get_quarts(nquarts,n_not_exited);

  dfrnd_empl_tot_zero_ne = dfrnd_empl_tot_zero(fnot_exited);
  dfsales_zero_ne = dfsales_zero(fnot_exited);
  dfrnd_to_sales_zero_ne = eq.w*dfrnd_empl_tot_zero_ne./dfsales_zero_ne;

  % get quart sets
  [ftab_rnd_sales,frnd_sales_rank] = sort(dfrnd_to_sales_zero_ne);
  quart_sets_rnd = get_quart_sets(frnd_sales_rank,qrange,nquarts);

  % growth rates
  q_sales_growth = zeros(1,nquarts);
  q_empl_growth = zeros(1,nquarts);
  for q=1:nquarts
    q_fsales_growth_sel = quart_sets_rnd(:,q);

    q_fsales_growth = dfsales_growth(q_fsales_growth_sel);
    q_sales_growth(q) = mean(q_fsales_growth);

    q_fempl_growth = dfempl_growth(q_fsales_growth_sel);
    q_empl_growth(q) = mean(q_fempl_growth);

    %add_moment(q_sales_growth(q),sprintf('sales_growth_q%i',q),1.0);
  end

  add_moment(q_empl_growth(2),'Empl growth R&D Q2');
  add_moment(q_empl_growth(3),'Empl growth R&D Q3');

  %add_moment(q_sales_growth(2),'Sales growth, R\&D quartile 2',0.5);
  %add_moment(q_sales_growth(3),'Sales growth, R\&D quartile 3',0.5);

  %disp(mean(q_sales_growth(1:2)));
  %disp(mean(q_sales_growth(3:4)));

  % aggregate growth rate
  add_moment(eq.g,'Aggregate Growth');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% firm size or industry presence breakdown                                %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % industry presence selection
  small_sel = dfmpos_zero <= alg.m_cutoff;
  large_sel = dfmpos_zero > alg.m_cutoff;

  % firm size selection
  %med_empl_zero = median(dfempl_tot_zero);
  %small_sel = dfempl_tot_zero <= med_empl_zero;
  %large_sel = dfempl_tot_zero > med_empl_zero;

  % total rnd selection
  %med_rnd_zero = median(dfrnd_empl_tot);
  %small_sel = dfrnd_empl_tot <= med_rnd_zero;
  %large_sel = dfrnd_empl_tot > med_rnd_zero;

  small_sel_ne = small_sel(fnot_exited);
  large_sel_ne = large_sel(fnot_exited);

  % sales growth
  sales_growth_low_m = mean_wins(dfsales_growth(small_sel_ne),wlo,whi);
  sales_growth_high_m = mean_wins(dfsales_growth(large_sel_ne),wlo,whi);

  add_moment(sales_growth_low_m,'Sales growth, $m \le 3$');
  add_moment(sales_growth_high_m,'Sales growth, $m > 3$');

  % employment growth
  empl_growth_low_m = mean_wins(dfempl_growth(small_sel_ne),wlo,whi);
  empl_growth_high_m = mean_wins(dfempl_growth(large_sel_ne),wlo,whi);

  add_moment(empl_growth_low_m,'Empl growth m <= 3');
  add_moment(empl_growth_high_m,'Empl growth m > 3');

  % return on sales
  ros_low_m = mean_wins(dfros_zero(small_sel),wlo,whi);
  ros_high_m = mean_wins(dfros_zero(large_sel),wlo,whi);

  add_moment(ros_low_m,'Return on sales $m \le 3$');
  add_moment(ros_high_m,'Return on sales $m > 3$');

  % exit rate
  exit_rate_low_m = mean(lfexited(small_sel));
  exit_rate_high_m = mean(lfexited(large_sel));

  add_moment(exit_rate_low_m,'Exit rate m <= 3');
  add_moment(exit_rate_high_m,'Exit rate m > 3');

  % R&D to sales
  rnd_to_sales_low_m = mean_wins(dfrnd_to_sales_zero(small_sel),wlo,whi);
  rnd_to_sales_high_m = mean_wins(dfrnd_to_sales_zero(large_sel),wlo,whi);

  add_moment(rnd_to_sales_low_m,'rnd_to_sales_low_m');
  add_moment(rnd_to_sales_high_m,'rnd_to_sales_high_m');

  % industry expansion
  expansion_low_m = mean(lfexpanded(small_sel_ne));
  expansion_high_m = mean(lfexpanded(large_sel_ne));

  add_moment(expansion_low_m,'Expansion, $m \le 3$');
  add_moment(expansion_high_m,'Expansion, $m > 3$');

  % positive growth
  pos_growth_low_m = mean(lfpos_growth(small_sel_ne));
  pos_growth_high_m = mean(lfpos_growth(large_sel_ne));
  pos_growth = mean(lfpos_growth);

  %add_moment(pos_growth_low_m,'pos_growth_low_m',1.0);
  %add_moment(pos_growth_high_m,'pos_growth_high_m',1.0);
  add_moment(pos_growth,'Fraction positive growth');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% some extras, value etc                                                  %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % value/employee ratio -> high_m to low_m
  %val_empl_basic = mean(dfval_to_empl(lfdo_basic));
  %val_empl_nobasic = mean(dfval_to_empl(not(lfdo_basic)));
  %val_ratio = val_empl_basic./val_empl_nobasic;
  val_ratio = 0.0;

  add_moment(val_ratio,'val_ratio');

  % correlation between val growth and bri
  valg_bri_corr_mat = corrcoef(dfval_growth,rnd_int(fnot_exited));
  valg_bri_corr = valg_bri_corr_mat(1,2);

  add_moment(valg_bri_corr,'valg_bri_corr');

  % industry expansion
  expansion_all = mean(lfexpanded);

  add_moment(expansion_all,'Expansion Rate');

  % within industry spillover
  spill_diff = (1.0/p.zeta)*(eq.taua/eq.tau);

  add_moment(spill_diff,'Spillover Differential');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% NEW MOMENTS - 2012/02/27                                                %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % employment groups
  empl_median = median(dfempl_tot_zero);
  empl_small = dfempl_tot_zero < empl_median;
  empl_large = dfempl_tot_zero >= empl_median;

  % age groups
  age_median = median(dfage_year);
  age_young = dfage_year < age_median;
  age_old = dfage_year >= age_median;

  % industry groups
  sel_m1 = dfmpos_zero == 1;
  sel_mp = dfmpos_zero > 1;

  % R&D intensity groups
  rnd_median = median(dfrnd_to_sales_zero);
  rnd_high = dfrnd_to_sales_zero >= rnd_median;
  rnd_low = dfrnd_to_sales_zero < rnd_median;

  % m by employment size
  mean_m_small = mean(dfmpos_zero(empl_small));
  mean_m_large = mean(dfmpos_zero(empl_large));
  add_moment(mean_m_small,'Industries, small firms');
  add_moment(mean_m_large,'Industries, large firms');

  % R&D expenditures or labor costs
  mean_rnd_labor = mean_wins(dfrnd_empl_tot_zero./dflabor_zero,0.0,0.95);
  % mean_rnd_labor = median(dfrnd_empl_tot_zero./dflabor_zero);
  add_moment(mean_rnd_labor,'R\&D/Labor');

  % Employment growth by basic/non-basic
  mean_empl_growth_nb = mean_wins(dfempl_growth.*(1.0-eq.epb(dfmind_zero(fnot_exited))),wlo,whi)/mean(1.0-eq.epb(dfmind_zero(fnot_exited)));
  mean_empl_growth_b = mean_wins(dfempl_growth.*eq.epb(dfmind_zero(fnot_exited)),wlo,whi)/mean(eq.epb(dfmind_zero(fnot_exited)));
  add_moment(mean_empl_growth_nb,'Employment growth, no basic');
  add_moment(mean_empl_growth_b,'Employment growth, basic');

  % Employment growth by m=1/m>1
  mean_empl_growth_m1 = mean_wins(dfempl_growth(sel_m1(fnot_exited)),wlo,whi);
  mean_empl_growth_mp = mean_wins(dfempl_growth(sel_mp(fnot_exited)),wlo,whi);
  add_moment(mean_empl_growth_m1,'Employment growth, m=1');
  add_moment(mean_empl_growth_mp,'Employment growth, m>1');

  % Exit by employment size
  mean_exit_small = mean(fexited(empl_small));
  mean_exit_large = mean(fexited(empl_large));
  add_moment(mean_exit_small,'Exit rate, small firms');
  add_moment(mean_exit_large,'Exit rate, large firms');

  % Age by employment size
  mean_age_small = mean(dfage_year(empl_small));
  mean_age_large = mean(dfage_year(empl_large));
  add_moment(mean_age_small,'Age, Small Firms');
  add_moment(mean_age_large,'Age, Large Firms');

  %% 2012/02/29 %%

  % Employment growth by firm size
  mean_empl_growth_small = mean_wins(dfempl_growth(empl_small(fnot_exited)),wlo,whi);
  mean_empl_growth_large = mean_wins(dfempl_growth(empl_large(fnot_exited)),wlo,whi);
  add_moment(mean_empl_growth_small,'Employment growth, small');
  add_moment(mean_empl_growth_large,'Employment growth, large');

  % Employment growth by firm age
  mean_empl_growth_young = mean_wins(dfempl_growth(age_young(fnot_exited)),wlo,whi);
  mean_empl_growth_old = mean_wins(dfempl_growth(age_old(fnot_exited)),wlo,whi);
  add_moment(mean_empl_growth_young,'Employment growth, young');
  add_moment(mean_empl_growth_old,'Employment growth, old');

  % Employment growth by R&D intensity
  mean_empl_growth_high = mean_wins(dfempl_growth(rnd_high(fnot_exited)),wlo,whi);
  mean_empl_growth_low = mean_wins(dfempl_growth(rnd_low(fnot_exited)),wlo,whi);
  add_moment(mean_empl_growth_high,'Employment growth, high');
  add_moment(mean_empl_growth_low,'Employment growth, low');

  % R&D intensity by firms size
  mean_rnd_small = mean_wins(dfrnd_to_sales_zero(empl_small),wlo,whi);
  mean_rnd_large = mean_wins(dfrnd_to_sales_zero(empl_large),wlo,whi);
  add_moment(mean_rnd_small,'R&D intensity, small');
  add_moment(mean_rnd_large,'R&D intensity, large');

  % Employment growth (central)
  empl_growth_central = dfempl_change(fnot_exited)./(0.5*(dfempl_tot_zero(fnot_exited)+dfempl_tot(fnot_exited)));
  mean_empl_growth_central = mean_wins(empl_growth_central,wlo,whi);
  add_moment(mean_empl_growth_central,'Employment growth (central)');

  % Sales growth (central)
  sales_growth_central = (dfsales_change(fnot_exited)+dfsales(fnot_exited)*eq.g)./(0.5*(dfsales(fnot_exited)*(1.0+eq.g)+dfsales_zero(fnot_exited)));
  mean_sales_growth_central = mean_wins(sales_growth_central,wlo,whi);
  add_moment(mean_sales_growth_central,'Sales growth (central)');

  % Sales growth
  sales_growth = dfsales_change(fnot_exited)./dfsales_zero(fnot_exited)+eq.g;
  mean_sales_growth = mean_wins(sales_growth,wlo,whi);
  add_moment(mean_sales_growth,'Sales growth');

  % total basic / total applied
  basic_applied_ratio = (sum(eq.cb.*eq.phivec)+eq.cz)/sum(eq.ca.*eq.phivec);
  add_moment(basic_applied_ratio,'Total basic/applied');

  % public basic / GDP
  public_basic_gdp = p.That;
  add_moment(public_basic_gdp,'Public basic / GDP');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % generality simulation   %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  zgeomh = eq.taua*p.citx*p.eta/eq.taub;
  dgeomh = 1.0/(1.0+zgeomh);

  zgeomc = eq.taua*p.citx*p.alpha/eq.taub;
  dgeomc = 1.0/(1.0+zgeomc);

  % generality
  [igenhm,egenhm,pgenhm] = gensim(p.pb,dgeomh,true);
  [igencm,egencm,pgencm] = gensim(p.pa,dgeomc,true);

  % aggregate over industry presence
  igenh = sum(eq.phivec.*igenhm);
  igenc = sum(eq.phivec.*igencm);
  egenh = sum(eq.phivec.*egenhm);
  egenc = sum(eq.phivec.*egencm);
  pgenh = sum(eq.phivec.*pgenhm);
  pgenc = sum(eq.phivec.*pgencm);

  % aggregate over hot/cold
  igena = eq.hot*igenh + (1.0-eq.hot)*igenc;
  igenb = igenh;
  egena = eq.hot*egenh + (1.0-eq.hot)*egenc;
  egenb = egenh;
  pgena = eq.hot*pgenh + (1.0-eq.hot)*pgenc;
  pgenb = pgenh;

  % aggregate over public/private
  igenprv = (eq.taua*igena+eq.taub*igenb)/eq.tau;
  igenpub = igenb;
  egenprv = (eq.taua*egena+eq.taub*egenb)/eq.tau;
  egenpub = egenb;
  pgenprv = (eq.taua*pgena+eq.taub*pgenb)/eq.tau;
  pgenpub = pgenb;

  %{
  fprintf(1,'intensive:\n');
  fprintf(1,'%6.4f %6.4f\n',igenc,igenh);
  fprintf(1,'%6.4f %6.4f\n',igena,igenb);
  fprintf(1,'%6.4f %6.4f\n',igenprv,igenpub);
  fprintf(1,'extensive:\n');
  fprintf(1,'%6.4f %6.4f\n',egenc,egenh);
  fprintf(1,'%6.4f %6.4f\n',egena,egenb);
  fprintf(1,'%6.4f %6.4f\n',egenprv,egenpub);
  fprintf(1,'probability:\n');
  fprintf(1,'%6.4f %6.4f\n',pgenc,pgenh);
  fprintf(1,'%6.4f %6.4f\n',pgena,pgenb);
  fprintf(1,'%6.4f %6.4f\n',pgenprv,pgenpub);
  %}

  % store moments
  add_moment(igenprv,'Private generality (intensive)');
  add_moment(igenpub,'Public generality (intensive)');
  add_moment(egenprv,'Private generality (extensive)');
  add_moment(egenpub,'Public generality (extensive)');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % citation simulation     %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (alg.nobasic == 1)
    taubs = eq.taus;
  else
    taubs = eq.taub;
  end

  zgeomh = eq.taua*p.citx*p.eta/taubs;
  dgeomh = zgeomh/(1.0+zgeomh);

  zgeomc = eq.taua*p.citx*p.alpha/taubs;
  dgeomc = zgeomc/(1.0+zgeomc);

  % simple basic spillovers
  cith1 = dgeomh/(1.0-dgeomh);
  cith2 = dgeomh*(1.0+dgeomh)/(1.0-dgeomh)^2;

  citc1 = dgeomc/(1.0-dgeomc);
  citc2 = dgeomc*(1.0+dgeomc)/(1.0-dgeomc)^2;

  citb1 = cith1;
  citb2 = cith2;

  cita1 = eq.hot*cith1 + (1.0-eq.hot)*citc1;
  cita2 = eq.hot*cith2 + (1.0-eq.hot)*citc2;

  % aggregate to public/private
  citpub1 = citb1;
  citpub2 = citb2;
  citpri1 = (eq.taua*cita1+eq.taub*citb1)/eq.tau;
  citpri2 = (eq.taua*cita2+eq.taub*citb2)/eq.tau;

  add_moment(citpub1,'Public Citations Mean');
  add_moment(sqrt(citpub2),'Public Citations RMS');

  add_moment(citpri1,'Private Citations Mean');
  add_moment(sqrt(citpri2),'Private Citations RMS');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % store moment values     %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % number of moments
  nm = length(m_vec);

  % load in data moments
  [m_dat,m_name] = parse_params(alg.targ_file);
  [m_wgt,m_wnam] = parse_params(alg.wgtvec_file);
  m_diff = m_vec-m_dat;

  % moments we use
  wlist = m_wgt > 0.0;
  m_vec_wl = m_vec(wlist);
  m_dat_wl = m_dat(wlist);
  m_wgt_wl = m_wgt(wlist);
  m_diff_wl = m_diff(wlist);
  nm_wl = length(m_vec_wl);

  % weighting matrix
  if (alg.wgt_opt)
    wgtmat = load(alg.wgtmat_file); % optimal
  else
    wgtmat = diag(m_wgt_wl.*m_dat_wl.^(-2)); % proportional
  end

  % find criterion function
  score = m_diff_wl*wgtmat*m_diff_wl';
  %score = sum(abs(m_diff_wl)./m_dat_wl);

  % proportional moments errors
  m_err = m_diff./m_dat;

  if (check == 1)
    save(alg.moments_file,'m_vec_wl','-ascii','-double');

    % construct cell array
    mcell = cell(nm,7);
    mcell(:,1) = num2cell(m_vec);
    mcell(:,2) = num2cell(m_dat);
    mcell(:,3) = num2cell(1:nm);
    mcell(:,4) = m_desc;
    mcell(:,5) = m_name;
    mcell(:,6) = num2cell(m_wgt);
    mcell(:,7) = num2cell(m_err);

    % moments text file
    mtxt_fid = fopen('output/moments_format.txt','w');
    fprintf(mtxt_fid,'%10s %10s %5s %30s %30s %10s %10s\n','Model','Data','#','Description1','Description2','Weight','Error');
    for i=1:nm
      fprintf(mtxt_fid,'%10.6f %10.6f %5i %30s %30s %10.6f %10.6f\n',mcell{i,1},mcell{i,2},mcell{i,3},mcell{i,4},mcell{i,5},mcell{i,6},mcell{i,7});
    end
    fprintf(mtxt_fid,'\n');
    fprintf(mtxt_fid,'score = %14.6f\n',score);
    fprintf(mtxt_fid,'\n');
    fclose(mtxt_fid);

    mtxt_fid = fopen('output/moments_used.txt','w');
    fprintf(mtxt_fid,'%10s %10s %5s %30s %10s %10s\n','Model','Data','#','Description','Weight','Error');
    for i=find(wlist)
      fprintf(mtxt_fid,'%10.6f %10.6f %5i %30s %10.6f %10.6f\n',mcell{i,1},mcell{i,2},i,mcell{i,4},mcell{i,6},mcell{i,7});
    end
    fprintf(mtxt_fid,'\n');
    fprintf(mtxt_fid,'score = %14.6f\n',score);
    fprintf(mtxt_fid,'\n');
    fclose(mtxt_fid);

    % moments csv file
    writetable(cell2table(mcell(wlist,[4,1,2,6,7]),'VariableNames',{'Description','Model','Data','Weight','Error'}),'output/moments_used.csv');

    % raw params file
    np = 21;
    pcell = {};
    pcell(1,:) = {'sigma',p.crra};
    pcell(2,:) = {'epsilon',p.epsn};
    pcell(3,:) = {'pa',p.pa};
    pcell(4,:) = {'pb',p.pb};
    pcell(5,:) = {'eta',p.eta};
    pcell(6,:) = {'alpha',p.alpha};
    pcell(7,:) = {'massout',p.massout};
    pcell(8,:) = {'massac',p.massac};
    pcell(9,:) = {'agamma',p.agamma-1.0};
    pcell(10,:) = {'bgamma',p.bgamma-1.0};
    pcell(11,:) = {'asigma',p.asigma_ns};
    pcell(12,:) = {'bsigma',p.bsigma_ns};
    pcell(13,:) = {'kappa',p.kappa};
    pcell(14,:) = {'dmu',-p.dmu};
    pcell(15,:) = {'dsig',p.dsig};
    pcell(16,:) = {'zeta',p.zeta};
    pcell(17,:) = {'qmin',p.qmin};
    pcell(18,:) = {'disc',p.disc};
    pcell(19,:) = {'merge',p.merge};
    pcell(20,:) = {'ufrac',p.ufrac};
    pcell(21,:) = {'x',p.citx};

    ptxt_fid = fopen('output/params_raw.txt','w');
    for i=1:np
      fprintf(ptxt_fid,'%20s : %15.12f\n',pcell{i,1},pcell{i,2});
    end
    fclose(ptxt_fid);

    ne = 6;
    ecell = {};
    ecell(1,:) = {'taua',eq.taua};
    ecell(2,:) = {'taub',eq.taub};
    ecell(3,:) = {'taus',eq.taus};
    ecell(4,:) = {'mrate',eq.mrate};
    ecell(5,:) = {'w',eq.w};
    ecell(6,:) = {'g',eq.g};

    etxt_fid = fopen('output/eqvars_raw.txt','w');
    for i=1:ne
      fprintf(etxt_fid,'%20s : %15.12f\n',ecell{i,1},ecell{i,2});
    end
    fclose(etxt_fid);

    %disp(score);
    %disp(n_used);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot some stuff         %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % bri
  t_epb = m_dat(4:11);
  t_bri = m_dat(12:19);

  eq.sim_epb = d_epb;
  eq.dat_epb = t_epb;
  eq.sim_bri = d_bri;
  eq.dat_bri = t_bri;

  % industry production distribution
  for m=1:7
    d_psipdf(m) = sum(dfmpos_zero == m);
  end
  d_psipdf(8) = sum(dfmpos_zero >= 8);
  d_psipdf = d_psipdf/sum(d_psipdf);

  z_psipdf(1:7) = eq.psipdf(1:7);
  z_psipdf(8) = sum(eq.psipdf(8:10));

  % true dists
  for m=1:p.M
    sel_m = (dfmpos_zero == m);
    base_psivec(m) = sum(sel_m);
    base_phivec(m) = sum(dfnprod_zero(sel_m));
    base_emplvec(m) = sum(dfempl_tot_zero(sel_m));
  end
  base_psivec = base_psivec/sum(base_psivec);
  base_phivec = base_phivec/sum(base_phivec);
  base_emplvec = base_emplvec/sum(base_emplvec);

  z_psivec(1:7) = base_psivec(1:7);
  z_psivec(8) = sum(base_psivec(8:10));

  z_phivec(1:7) = base_phivec(1:7);
  z_phivec(8) = sum(base_phivec(8:10));

  z_emplvec(1:7) = base_emplvec(1:7);
  z_emplvec(8) = sum(base_emplvec(8:10));

  %{
  function plot_ks(dat)
    [f,xi] = ksdensity(dat);
    plot(xi,f,'LineWidth',2);
    %hist(dat,50);
  end

  close all;
  figure;
  subplot(2,2,1), plot_ks(log(dfempl_tot)), title('employment');
  subplot(2,2,2), bar(d_psipdf), title('industry production');
  subplot(2,2,3), plot_ks(dfros_zero(dfros_zero>-0.5)), title('return on sales');
  subplot(2,2,4), plot_ks(dfrnd_to_sales_zero(dfrnd_to_sales_zero<1.0)), title('R&D to sales');
  %}

  %{
  figure;
  subplot(2,2,1), plot([m_bri' d_bri']);
  subplot(2,2,2), plot([m_epb' d_epb']);
  subplot(2,2,3), plot_ks(dfempl_tot), title('employment');
  subplot(2,2,4), bar([z_psipdf' d_psipdf']), title('industry presence');
  %}

  % basic research graphs

  % disp([d_bri; t_bri]');
  % disp([d_epb; t_epb]');

  %{
  set(0,'DefaultAxesColorOrder',[1 0 0; 0 0 1]);

  figure();
  plot([d_epb; t_epb]','LineWidth',2.0);
  legend('Model','Data','Location','NorthWest');
  xlim([0.5 8.5]);
  ylim([0.0 0.8]);
  xlabel('Number of Industries');
  ylabel('Fraction Positive Basic')

  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','inches');
  set(gcf,'PaperPosition',[0.0 0.0 3.5 3.0]);
  print('-depsc','latex/posbas.eps','-r400');

  figure();
  plot([d_bri; t_bri]','LineWidth',2.0);
  legend('Model','Data','Location','NorthWest');
  xlim([0.5 8.5]);
  ylim([0.0 0.15]);
  xlabel('Number of Industries');
  ylabel('Basic Research Intensity');

  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','inches');
  set(gcf,'PaperPosition',[0.0 0.0 3.5 3.0]);
  print('-depsc','latex/basint.eps','-r400');
  %}

  % qhat distribution graphs
  %{
  figure;
  nbins = qd.nB;
  plot(qd.binmids(2:nbins),[qd.even_pmf(2:nbins)' qd.step_a_pmf(2:nbins)' qd.step_b_pmf(2:nbins)'],'LineWidth',2);
  legend('Production','Applied','Basic','Location','NorthEast');
  xlim([0.3 2.1]);

  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','inches');
  set(gcf,'PaperPosition',[0.0 0.0 3.5 3.0]);
  print('-depsc','latex/qdists.eps','-r400');

  % distribution over m or n?
  figure;
  bar([z_psivec' z_emplvec']);
  legend('Firm Share','Employment Share');
  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','inches');
  set(gcf,'PaperPosition',[0.0 0.0 3.5 3.0]);
  print('-depsc','latex/m_dist.eps','-r400');

  figure;
  bar(eq.npdf);
  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','inches');
  set(gcf,'PaperPosition',[0.0 0.0 3.5 3.0]);
  print('-depsc','latex/n_dist.eps','-r400');

  % log employment dist
  figure;
  [vals,bins] = hist(log(dfempl_tot),50);
  vals = vals/sum(vals);
  bar(bins,vals);
  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','inches');
  set(gcf,'PaperPosition',[0.0 0.0 3.5 3.0]);
  print('-depsc','latex/empl_dist.eps','-r400');

  % R&D intensity dist
  figure;
  hist(log(dfrnd_to_sales_zero),50);
  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','inches');
  set(gcf,'PaperPosition',[0.0 0.0 3.5 3.0]);
  print('-depsc','latex/rndi_dist.eps','-r400');

  % beta function graphs
  betavec = eq.pibar*((qd.binmids_eps-p.qminp)/eq.denom1+p.qminp/eq.denom2);
  plot(qd.binmids,betavec);
  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','inches');
  set(gcf,'PaperPosition',[0.0 0.0 3.5 3.0]);
  print('-depsc','latex/beta_func.eps','-r400');
  %}

  %{
  fprintf(1,'corr(bri,ros) = %f\n',corr(rnd_int',dfros_zero'));
  fprintf(1,'corr(bri,empl) = %f\n',corr(rnd_int',dfempl_tot_zero'));
  fprintf(1,'corr(bri,lempl) = %f\n',corr(rnd_int',log(dfempl_tot_zero')));
  fprintf(1,'corr(ind,empl) = %f\n',corr(dfmpos_zero,dfempl_tot_zero'));
  fprintf(1,'corr(ind,lempl) = %f\n',corr(dfmpos_zero,log(dfempl_tot_zero')));
  fprintf(1,'corr(bri,sales) = %f\n',corr(rnd_int',dfsales_zero'));
  fprintf(1,'corr(bri,lsales) = %f\n',corr(rnd_int',log(dfsales_zero')));
  fprintf(1,'corr(ind,sales) = %f\n',corr(dfmpos_zero,dfsales_zero'));
  fprintf(1,'corr(ind,lsales) = %f\n',corr(dfmpos_zero,log(dfsales_zero')));
  fprintf(1,'skewness(empl) = %f\n',skewness(dfempl_tot_zero));
  %}

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% utility functions       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qrout = get_quarts(nq,nel)
  qdivs = floor((1:nq)/nq*nel);
  qrout = zeros(nq,2);
  qrout(1,1) = 1;
  for q=1:(nq-1)
    qrout(q,2) = qdivs(q);
    qrout(q+1,1) = qrout(q,2)+1;
  end
  qrout(nq,2) = nel;
end

function qsets = get_quart_sets(rank_vec,qrange,nquarts)
  qsets = logical(zeros(length(rank_vec),nquarts));
  for q=1:nquarts
    qsets(rank_vec(qrange(q,1):qrange(q,2)),q) = 1;
  end
end

function add_moment(m_val,m_str)
  global m_vec m_desc m_pos

  m_vec(m_pos) = m_val;
  m_desc{m_pos} = m_str;
  m_pos = m_pos + 1;
end

function fprint_vec(fid,vec)
  if (length(vec) > 1)
    fprintf(fid,'%f ',vec(1:end-1));
  end
  fprintf(fid,'%f\n',vec(end));
end

function wvec = mean_wins(vec,qlo,qhi)
  vlo = quantile(vec,qlo);
  vhi = quantile(vec,qhi);
  wvec = mean(max(vlo,min(vhi,vec)));
end
