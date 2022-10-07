function gradmat
  global alg p eq

  % alg init
  alg = {};
  initalg();

  [par0,names] = parse_params(alg.par_file);
  alg.nparams = length(par0);
  eqv0 = load(alg.eqv_file);

  wgt = load(alg.wgtmat_file);
  alg.ntargs = length(wgt);

  alg.firstrun = 1;
  alg.disp_set = 'iter';
  alg.eqv0 = eqv0;

  pstep = 1e-4;
  stepvec = pstep*par0;

  gradmat = zeros(alg.ntargs,alg.nparams);
  for pnum=1:alg.nparams
    fprintf(1,'pnum = %i\n',pnum);

    parp = par0;
    parp(pnum) = par0(pnum) + stepvec(pnum);
    [eqfin,eqerr] = eqstand(parp);
    [score,eqdp1] = compfs();

    parp = par0;
    parp(pnum) = par0(pnum) - stepvec(pnum);
    [eqfin,eqerr] = eqstand(parp);
    [score,eqdp2] = compfs();

    gradmat(:,pnum) = (eqdp1-eqdp2)/(2.0*stepvec(pnum));
  end

  save(alg.gradmat_file,'gradmat','-ascii','-double');
end

