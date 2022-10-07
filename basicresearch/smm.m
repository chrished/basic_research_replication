function smm
  global p alg

  s = RandStream.create('mt19937ar','seed',sum(100*clock));
  RandStream.setGlobalStream(s);

  alg = {};
  initalg();

  [par0,names] = parse_params(alg.par_file);
  eqv0 = load(alg.eqv_file);

  alg.nparams = length(par0);
  alg.pscale = par0;
  alg.eqv0 = eqv0;

  alg.fpar = [3 20 17]; % baseline!
  % alg.fpar = [3 15 17 20]; % nobasic
  % alg.fpar = [3 4 5 8 10 12 14 15 16 20 21 22]; % noacad
  % alg.fpar = [1 3 20 17]; % sigma2
  alg.spar = setdiff(1:alg.nparams,alg.fpar);
  alg.nfpar = length(alg.fpar);
  alg.nspar = length(alg.spar);

  out_idx = 63;
  alg.efid = fopen(['evals/evals' num2str(out_idx) '.txt'],'a+');
  alg.mfid = fopen(['evals/mvals' num2str(out_idx) '.txt'],'a+');
  alg.disp_set = 'iter';
  alg.firstrun = 1;
  alg.besteval = Inf;

  start = ones(1,alg.nspar);

  alg.simann = 0; % otherwise use simplex
  if (alg.simann == 1)
    [parfin,scfin] = anneal0(@smmobj,start,0.25,1000000);
  else
    mopts = optimset('Display','iter','MaxFunEvals',1000000,'MaxIter',1000000);
    [parfin,scfin] = fminsearch(@smmobj,start,mopts);
  end

  params = alg.pscale;
  params(alg.spar) = parfin.*alg.pscale(alg.spar);

  fprintf(1,'%16.12f\n',params);
  fprintf(1,'\n');
  fprintf(1,'%16.12f\n',scfin);

  fclose(alg.efid);
  fclose(alg.mfid);

end

function score = smmobj(parin)
  global alg

  %fprintf(1,'\nSMMOBJ\n\n');

  params = alg.pscale;
  params(alg.spar) = parin.*alg.pscale(alg.spar);

  fprintf(alg.efid,'%15.10f ',params);
  fprintf(1,'%15.10f ',params);
  fprintf(1,'\n');

  [eqfin,eqerr] = eqstand(params);

  if (eqerr > 5000.0)
    score = 10000.0;
    fprintf(alg.mfid,'FAILED\n');
  else
    [score,moments] = compfs();

    fprintf(alg.mfid,'%15.10f ',moments);
    fprintf(alg.mfid,'\n');
  end

  fprintf(alg.efid,'%15.10f ',eqfin);
  fprintf(alg.efid,'%15.10f ',eqerr);
  fprintf(alg.efid,'%15.10f\n',score);
  fprintf(1,' %15.10f\n\n',score);

  if (isnan(score))
    score = 20000.0;
  end

  if (score < alg.besteval)
    alg.besteval = score;
    alg.eqv0 = alg.lasteq;
    unix('cp output/moments_used.txt output/moments_best.txt');
    unix('cp output/params_raw.txt output/params_best.txt');
    unix('cp output/eqvars_raw.txt output/eqvars_best.txt');
  end
end
