function [eqfin,err] = eqstand(params)
  global alg p

  % alg init
  alg.fevals = 0;
  alg.solve_fid = fopen(alg.solve_file,'w');

  % load in params
  p = {};
  loadparams(params);

  % determine initial guess
  if (isfield(alg,'eqv0'))
    eqv0 = alg.eqv0;
  else
    eqv0 = load(alg.eqv_file);
  end

  % set display mode
  if (~isfield(alg,'disp_set'))
    alg.disp_set = 'iter';
  end

  % solver for w
  options = optimset('Display',alg.disp_set,'MaxFunEvals',200,'TypicalX',eqv0,'DiffMinChange',1e-6);
  [eqfin,eqdiff,eflag] = fsolve(@eqfunc,eqv0,options);

  fclose(alg.solve_fid);

  if (eflag <= 0)
    disp('Failed.');
    eqfin = zeros(size(eqv0));
    err = 10000.0;
    return
  end

  eqout = eqfunc(eqfin);

  err = mean(abs(eqout));
  if (err > 1e-6)
    disp(['Failed. err = ' num2str(err)]);
    eqfin = zeros(size(eqv0));
    err = 10000.0;
    return
  end

  alg.lasteq = eqfin;

end

