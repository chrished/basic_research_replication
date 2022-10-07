function pick_best(evnum,ptag)
  global alg

  alg = {};
  initalg();
  [params,names] = parse_params(alg.par_file);
  npar = length(params);

  evname = ['evals/evals' num2str(evnum) '.txt'];

  ev_mat = importdata(evname);
  ev_sort = sortrows(ev_mat,size(ev_mat,2));

  brow = ev_sort(1,:);
  bpar = brow(1:npar)';
  beqv = brow(npar+1:end-2)';

  pname = ['params/params' num2str(ptag) '.txt'];
  eqname = ['eqvars/eqvars' num2str(ptag) '.txt'];

  save_params(pname,bpar,names);
  save(eqname,'beqv','-ascii','-double');

end

