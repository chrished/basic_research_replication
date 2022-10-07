function [score,moments] = targs_simple
  global alg p eq

  n_mmt = 25;

  brnd_ext(1:7) = eq.epb(1:7);
  brnd_ext(8) = sum(eq.epb(8:end).*eq.psipdf(8:end))/sum(eq.psipdf(8:end));

  brndi = eq.cb./eq.ca;
  brnd_int(1:7) = brndi(1:7);
  brnd_int(8) = sum(brndi(8:end).*eq.psipdf(8:end))/sum(eq.psipdf(8:end));

  mean_m = (1:p.M)*eq.psipdf';
  mean_m2 = ((1:p.M).^2)*eq.psipdf';

  mean_ros = 1.0/p.epsn - eq.w*eq.cinc;

  exit_rate = eq.tau*eq.npdf(1);

  rnd_labor = eq.cinc/eq.prodlab;

  priv_public = sum(eq.cb.*eq.phivec)./eq.cz;

  agg_growth = eq.g;

  basic_applied_ratio = (sum(eq.cb.*eq.phivec)+eq.cz)/sum(eq.ca.*eq.phivec);

  public_basic_gdp = p.That;

  basic_applied = (sum(eq.cb.*eq.phivec)+eq.cz)/sum(eq.ca.*eq.phivec);

  sim_targs = [brnd_ext,brnd_int,mean_m,mean_m2,mean_ros,exit_rate,priv_public,agg_growth,basic_applied_ratio,public_basic_gdp,basic_applied];
  dat_targs = [0.24319389,0.23079438,0.27025862,0.34827862,0.41836735,0.45175439,0.55084746,0.68032787,0.0665785,0.0466688,0.0617295,0.0800356,0.0762766,0.0786363,0.1224092,0.1005365,2.203677,6.975635,0.032614,0.091857,0.115714,0.015000,0.640000,0.005000,0.640000];
  wgt_targs = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0];
  targ_names = {'Basic R&D extensive m=1','Basic R&D extensive m=2','Basic R&D extensive m=3','Basic R&D extensive m=4','Basic R&D extensive m=5','Basic R&D extensive m=6','Basic R&D extensive m=7','Basic R&D extensive m=8+','Basic R&D intensive m=1','Basic R&D intensive m=2','Basic R&D intensive m=3','Basic R&D intensive m=4','Basic R&D intensive m=5','Basic R&D intensive m=6','Basic R&D intensive m=7','Basic R&D intensive m=8+','Mean m','Mean m^2','Return on sales','Exit rate','Private/public basic','Aggregate growth','Total basic/applied','Public basic / GDP','Total basic / applied'};

  mmt_error = sim_targs - dat_targs;
  mmt_weight  = wgt_targs./abs(dat_targs);
  score = sum(mmt_weight.*abs(mmt_error));
  moments = sim_targs;

  mtxt_fid = fopen('output/moments_simple.txt','w+');
  fprintf(mtxt_fid,'%10s %10s %5s %30s %10s %10s\n','Model','Data','#','Description','Weight','Error');
  for i=1:n_mmt
    fprintf(mtxt_fid,'%10.6f %10.6f %5i %30s %10.6f %10.6f\n',sim_targs(i),dat_targs(i),i,targ_names{i},mmt_weight(i),mmt_error(i));
  end
  fprintf(mtxt_fid,'\n');
  fprintf(mtxt_fid,'score = %14.6f\n',score);
  fprintf(mtxt_fid,'\n');
  fclose(mtxt_fid);

end

