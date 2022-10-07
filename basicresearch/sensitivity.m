function sensitivity
  global alg

  % step size
  step = 0.05;

  % init alg
  alg = {};
  initalg();
  alg.disp_set = 'none';

  % initial eq guess
  [pvec0,names] = parse_params(alg.par_file);
  alg.eqv0 = load(alg.eqv_file);

  % get baseline values
  [eqfin0,errs0] = eqstand(pvec0);
  [score0,mvec0] = compfs();

  % allocate matrix
  nparm = length(pvec0);
  nmomt = length(mvec0);
  smat = zeros(nparm,nmomt);

  for i=1:nparm
    fprintf(1,'PARAM %i\n',i);
    parp = pvec0;
    vstep = step*parp(i);
    parp(i) = parp(i) + vstep;
    [eqfin,errs] = eqstand(parp);
    [score,mvec] = compfs();
    smat(i,:) = 100*(mvec-mvec0)./mvec0;
  end

  save('output/sensitivity.txt','smat','-ascii','-double');
end
