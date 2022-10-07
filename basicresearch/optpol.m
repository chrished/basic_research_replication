function polfin = optpol(polinit,ptype)
  global goodx

  if (strcmp(ptype,'uniform'))
    weltype = @welobj_uniform;
  elseif (strcmp(ptype,'targeted'))
    weltype = @welobj_targeted;
  elseif (strcmp(ptype,'academic'))
    weltype = @welobj_academic;
  elseif (strcmp(ptype,'academic_full'))
    weltype = @welobj_academic_full;
  elseif (strcmp(ptype,'academic_labs'))
    weltype = @welobj_academic_labs;
  elseif (strcmp(ptype,'uniform_academic'))
    weltype = @welobj_uniform_academic;
  elseif (strcmp(ptype,'targeted_academic'))
    weltype = @welobj_targeted_academic;
  elseif (strcmp(ptype,'uniform_noacademic'))
    weltype = @welobj_uniform_noacademic;
  elseif (strcmp(ptype,'basic'))
    weltype = @welobj_basic;
  elseif (strcmp(ptype,'academic_zero'))
    weltype = @welobj_academic_zero;
  else
    disp('Invalid policy type.');
    return;
  end

  goodx = nan;
  weltype_recurse = recursify(weltype);
  polfin = fminsearch(weltype_recurse,polinit);
  weltype_recurse(polfin);
end

function ret = welobj_uniform(pol)
  asubs = pol(1);
  bsubs = pol(1);
  acfrac = nan;
  massac = nan;

  [wval,ceval,info] = welfare(asubs,bsubs,acfrac,massac,'');
  ret = -ceval;
end

function ret = welobj_targeted(pol)
  asubs = pol(1);
  bsubs = pol(2);
  acfrac = nan;
  massac = nan;

  [wval,ceval,info] = welfare(asubs,bsubs,acfrac,massac,'');
  ret = -ceval;
end

function ret = welobj_academic(pol)
  asubs = nan;
  bsubs = nan;
  acfrac = exp(pol(1));
  massac = nan;

  if (acfrac < 0.0)
    ret = Inf;
  else
    [wval,ceval,info] = welfare(asubs,bsubs,acfrac,massac,'');
    ret = -ceval;
  end
end

function ret = welobj_academic_full(pol)
  asubs = nan;
  bsubs = nan;
  acfrac = exp(pol(1));
  massac = exp(pol(2));

  if ((acfrac < 0.0) || (massac < 0.0))
    ret = Inf;
  else
    [wval,ceval,info] = welfare(asubs,bsubs,acfrac,massac,'');
    ret = -ceval;
  end
end

function ret = welobj_academic_labs(pol)
  asubs = nan;
  bsubs = nan;
  acfrac = nan;
  massac = exp(pol(1));

  if (massac < 0.0)
    ret = Inf;
  else
    [wval,ceval,info] = welfare(asubs,bsubs,acfrac,massac,'');
    ret = -ceval;
  end
end

function ret = welobj_uniform_academic(pol)
  asubs = pol(1);
  bsubs = pol(1);
  acfrac = exp(pol(2));
  massac = nan;

  if (acfrac < 0.0)
    ret = Inf;
  else
    [wval,ceval,info] = welfare(asubs,bsubs,acfrac,massac,'');
    ret = -ceval;
  end
end

function ret = welobj_targeted_academic(pol)
  asubs = pol(1);
  bsubs = pol(2);
  acfrac = exp(pol(3));
  massac = nan;

  if (acfrac < 0.0)
    ret = Inf;
  else
    [wval,ceval,info] = welfare(asubs,bsubs,acfrac,massac,'');
    ret = -ceval;
  end
end

function ret = welobj_uniform_noacademic(pol)
  asubs = pol(1);
  bsubs = pol(1);
  acfrac = 0.0;
  massac = nan;

  [wval,ceval,info] = welfare(asubs,bsubs,acfrac,massac,'');
  ret = -ceval;
end

function ret = welobj_basic(pol)
  z = 0.5; % fraction of applied that gets reported as basic
  % z = 11.04970949/48.64275441; % test a precisely zbar

  asubs = z*pol(1);
  bsubs = pol(1);
  acfrac = nan;
  massac = nan;

  [wval,ceval,info] = welfare(asubs,bsubs,acfrac,massac,'');
  ret = -ceval;
end

function ret = welobj_academic_zero(pol)
  asubs = 0.0;
  bsubs = 0.0;
  acfrac = exp(pol(1));
  massac = nan;

  if (acfrac < 0.0)
    ret = Inf;
  else
    [wval,ceval,info] = welfare(asubs,bsubs,acfrac,massac,'');
    ret = -ceval;
  end
end

% down and dirty
function ret = recurse_objective(func,lastx,newx,rlevel)
  fprintf('ATTEMPTING:');
  fprintf(' %f',newx);
  fprintf('\n');
  ret = func(newx);
  if (~any(isnan(lastx)) && isnan(ret))
    midx = 0.5*(lastx+newx);
    fprintf('OBJECTIVE FAILURE AT RECURSION LEVEL %i:',rlevel);
    fprintf(' %f',newx);
    fprintf('\n');
    ret = recurse_objective(func,lastx,midx,rlevel+1);
  end
end

function ret = recurse_objective_toplevel(func,x)
  global goodx
  if (length(goodx) == 0)
    goodx = nan;
  end
  ret = recurse_objective(func,goodx,x,0);
  goodx = x;
end

function func_new = recursify(func)
  func_new = @(x) recurse_objective_toplevel(func,x);
end
