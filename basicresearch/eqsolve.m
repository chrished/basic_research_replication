function eqsolve(save_eqvars)
  global alg p eq qd

  if (nargin < 1)
    save_eqvars = 0;
  end

  % alg init
  alg = {};
  initalg();

  % set parameters
  p = {};
  [par0,names] = parse_params(alg.par_file);
  eqv0 = load(alg.eqv_file);
  loadparams(par0);

  % solver for w
  options = optimset('Display','iter','MaxFunEvals',1000,'TypicalX',eqv0,'DiffMinChange',1e-6);
  [eqfin,eqdiff,eflag] = fsolve(@eqfunc,eqv0,options);

  if (eflag <= 0)
    disp('Failed.');
    return
  end

  fprintf(1,'\neqvars:\n'); fprintf(1,'%15.10f\n',eqfin); fprintf(1,'\n');
  save('logs/eqvars.txt','eqfin','-ascii','-double');
  if (save_eqvars == 1)
    unix(['cp logs/eqvars.txt ' alg.eqv_file]);
  end

  alg.check = 1;
  eqout = eqfunc(eqfin);

  fprintf(1,'\nerr:\n'); fprintf(1,'%e\n',eqout); fprintf(1,'\n');

  err = mean(abs(eqdiff));
  if (err > 1e-7)
    fprintf(1,'Warning: err = %e\n\n',err);
  end

  % moments
  %[score,alg.mvec] = compfs();
  %fprintf(1,'score = %10.5f\n',score);

  % save these for the consumption equiv calculations
  cew = eq.w;
  ceg = eq.g;
  save('logs/cevals.txt','cew','ceg','-ascii','-double');

  % save for cons equiv later
  if (save_eqvars == 1)
    save('logs/optvals.txt','cew','ceg','-ascii','-double');
  end

  % display stats
  output();

  % save equilibrium stats
  % save('workspace.mat');
end
