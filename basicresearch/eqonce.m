function eqonce()
  global alg p eq qd

  % algo init
  alg = {};
  initalg();

  % load params
  p = {};
  [par0,names] = parse_params(alg.par_file);
  loadparams(par0);

  % load eqvars
  eqvin = load(alg.eqv_file);
  fprintf(1,'eqvars: '); fprintf(1,'%10.8f ',eqvin); fprintf(1,'\n');

  % find errors
  eqout = eqfunc(eqvin);
  fprintf(1,'errors: '); fprintf(1,'%10.8f ',eqout); fprintf(1,'\n');
end
