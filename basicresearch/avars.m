function avars
  global alg p

  alg = {};
  initalg();

  [params,names] = parse_params(alg.par_file);
  loadparams(params);
  moments = load(alg.moments_file);
  grad = load(alg.gradmat_file);
  wgt = load(alg.wgtmat_file);
  cov = load(alg.covmat_file);

  %nobs = 13705; % num firm-years
  nobs = 6775; % num firms

  N = nobs;

  npar = length(params);
  nmom = length(wgt);

  % don't work
  spar = setdiff(1:npar,[]);
  smom = setdiff(1:nmom,[21]); % no R&D/sales

  G = grad(smom,spar);
  W = wgt(smom,smom);
  S = cov(smom,smom);
  
  P = inv(G'*W*G);
  disp(det(G'*W*G));

  avmat = P*G'*W*S*W*G*P;
  
  avar = zeros(p.nparm,1);
  avar(spar) = sqrt(diag(avmat)/N);
  %avar(spar) = sqrt(diag(avmat));
  avarFrac = avar./params; 
  avarPct = 100.0*avarFrac;

  disp([(1:npar)' params avar avarPct]);

  % params text file
  pcell = cell(p.nparm,5);
  pcell(:,1) = num2cell(1:p.nparm);
  pcell(:,2) = p.pnames;
  pcell(:,3) = num2cell(params);
  pcell(:,4) = num2cell(avar);
  pcell(:,5) = num2cell(avarPct);

  ptxt_fid = fopen('logs/avars_format.txt','w');
  fprintf(ptxt_fid,'%5s %30s %10s %12s %12s\n','#','Description','Value','Avar','Avar Pct');
  for i=1:p.nparm
    fprintf(ptxt_fid,'%5i %30s %10.6f %12.8f %12.8f\n',pcell{i,:});
  end
  fclose(ptxt_fid);

  % params latex
  ptab_fid = fopen('latex/params_table.tex','w');
  fprintf(ptab_fid,'\\documentclass{article}\n');
  fprintf(ptab_fid,'\\begin{document}\n');
  fprintf(ptab_fid,'\\begin{tabular}{|c|c|l|c|}\n');
  fprintf(ptab_fid,'\\hline\n');
  fprintf(ptab_fid,'\\textbf{Value} & \\textbf{Std. Err.} & \\textbf{Description} & \\textbf{Equation} \\\\ \\hline\\hline\n');
  for i=1:p.nparm
    fprintf(ptab_fid,'%8.5f & %8.5f & %s & \\\\ \\hline\n',pcell{i,3},pcell{i,4},pcell{i,2});
  end
  fprintf(ptab_fid,'\\end{tabular}\n');
  fprintf(ptab_fid,'\\end{document}\n');
  fprintf(ptab_fid,'\n');
  fclose(ptab_fid);

  % sensitivity matrix
  %sens = grad.*repmat(params',nmom,1)./repmat(moments(smom)',1,npar);
  %save('logs/sens.txt','sens','-ascii','-double');

end

