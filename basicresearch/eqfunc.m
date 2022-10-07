function eqcond = eqfunc(eqvec)
  global p alg eq qd

  if (any(eqvec <= 0.0))
    eqcond = -100000.0;
    return
  end

  alg.hfail = 0;
  alg.fid = fopen('logs/eqlog-tmp.txt','w+');
  fprintf(alg.fid,'asubs = %5.3f, bsubs = %5.3f, That = %5.3f\n',p.asubs,p.bsubs,p.That);

  % direct vars
  eq.taua = eqvec(1);
  eq.taub = eqvec(2);
  eq.taus = eqvec(3);
  eq.mrate = eqvec(4);
  eq.w = eqvec(5);
  eq.g = eqvec(6);

  if (alg.nobasic == 1)
    eq.taub = 0.0;
  end
  if ((alg.nobasic == 1) & (alg.noacad == 1))
    eq.taus = 0.0;
  end

  % indirect vars
  eq.tau = eq.taua + eq.taub;
  eq.r = eq.g*p.crra + p.disc;

  fprintf(alg.fid,'eq.taua = %15.10f\n',eq.taua);
  fprintf(alg.fid,'eq.taub = %15.10f\n',eq.taub);
  fprintf(alg.fid,'eq.taus = %15.10f\n',eq.taus);
  fprintf(alg.fid,'eq.mrate = %15.10f\n',eq.mrate);
  fprintf(alg.fid,'eq.w = %15.10f\n',eq.w);
  fprintf(alg.fid,'eq.g = %15.10f\n',eq.g);
  fprintf(alg.fid,'eq.r = %15.10f\n',eq.r);
  fprintf(alg.fid,'eq.tau = %15.10f\n',eq.tau);
  fprintf(alg.fid,'tau_heat = %15.10f\n',eq.taub+eq.taus);

  qdistbin(); % calculate distribution over q and hot/cold
  findbetas(); % find production values
  academic(); % academic innovation

  % solve value function - check for failure
  findh();
  if (alg.hfail == 1)
    unix('mv logs/eqlog-tmp.txt logs/eqlog.txt');
    eqcond = -100000.0;
    return
  end

  innovation(); % find innovation intensities
  nmdist(); % calculate the (n,m) distribution - as well as psivec and phivec
  entrants(); % solve for entrants' innovation rate
  labor(); % lobor demand
  growthrate(); % growth rate

  % aggregate rates
  eq.a0 = p.massout*eq.xout;
  eq.abar = sum((1.0+p.rhoa).*eq.xa.*eq.phivec);
  eq.bbar = sum((1.0+p.rhob).*eq.xb.*eq.phivec);
  eq.aspill = sum((p.pa-p.rhoa).*eq.xa.*eq.phivec);
  eq.bspill = sum((p.pb-p.rhob).*eq.xb.*eq.phivec);

  % equilibrium equations - (taua,taub,taus,labor,growth)
  eqcond(1) = eq.taua - eq.a0 - eq.abar;
  eqcond(2) = eq.taub - eq.bbar - eq.freebase;
  eqcond(3) = eq.taus - eq.aspill - (1-alg.bloss)*eq.bspill - eq.acspill;
  eqcond(4) = eq.fsize*eq.mrate - p.merge*eq.a0;
  eqcond(5) = 1.0 - eq.prodlab - eq.rndlab;
  eqcond(6) = 10.0*(eq.g - eq.newg);

  % fix wage at a certain level
  if (~isnan(alg.fixwage))
    eqcond(5) = eq.w - alg.fixwage;
  end
  if (alg.nobasic == 1)
    eqcond(2) = eqvec(2) - 0.1;
  end
  if ((alg.nobasic == 1) & (alg.noacad == 1))
    eqcond(3) = eqvec(3) - 0.1;
  end

  % solver log
  err = sum(eqcond.*eqcond);
  alg.fevals = alg.fevals + 1;
  %fprintf(alg.solve_fid,'%15.10f ',eqvec'); fprintf(alg.solve_fid,'\n');
  %fprintf(alg.solve_fid,'%15.10f ',eqcond); fprintf(alg.solve_fid,'\n');
  %fprintf(alg.solve_fid,'err = %11.6g, feval = %3i\n\n',err,alg.fevals);

  if (any(imag(eqcond)~=0))
    disp(eqcond');
  end

  fclose(alg.fid);
  unix('mv logs/eqlog-tmp.txt logs/eqlog.txt');

end
