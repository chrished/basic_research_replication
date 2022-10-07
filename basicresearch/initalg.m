function initalg
  global alg

  % alg init
  alg.hmax = 5.0;
  alg.prodmax = 64;
  alg.M = 10;
  alg.glimfact = 0.99; % in [0,1), increase if lim=1
  alg.rlim = 1; % limit total innovation intensities
  alg.citK = 10000;
  alg.check = 0;
  alg.noexp = 0;
  alg.fixwage = nan; % 0.968940077536654;
  alg.acmult = nan; % 5.0
  alg.nobasic = 0; % 0
  alg.noacad = 0; % 0
  alg.bintern = 0; % is private basic internalized
  alg.bloss = 0.0; % fraction of private basic kept secret
  alg.bdole = 0.0; % academic applicability

  % warnings
  if (alg.noexp == 1)
    fprintf(1,'EXPANSION SHUT DOWN\n');
  end
  if (~isnan(alg.fixwage))
    fprintf(1,'WAGE FIXED AT %10.8f\n',alg.fixwage);
  end
  if (~isnan(alg.acmult))
    fprintf(1,'AC MULT AT %4.2f\n',alg.acmult);
  end
  if (alg.nobasic == 1)
    fprintf(1,'BASIC RESEARCH SHUT DOWN\n');
  end
  if (alg.noacad == 1)
    fprintf(1,'ACADEMIC RESEARCH SHUT DOWN\n');
  end
  if (alg.bintern == 1)
    fprintf(1,'PRIVATE BASIC IS INTERNALIZED\n');
  end
  if (alg.bloss > 0)
    fprintf(1,'BLOSS SECRET SAUCE IS NON-ZERO\n');
  end
  if (alg.bdole > 0)
    fprintf(1,'BOB DOLE IS IN THE HOUSE\n');
  end

  % state info
  alg.fevals = 0;
  alg.solve_file = 'logs/solvelog.txt';
  alg.solve_fid = 1;

  % options
  alg.par_idx = '1_3'; % param index
  alg.targ_idx = '15'; % '15' % data target index
  alg.wgt_idx = '3'; % '3' % weighting vec used
  alg.m_cutoff = 3; % cutoff for m for split moments
  alg.wgt_opt = 0; % use optimal weighting matrix?

  % filenames
  alg.par_file = ['params/params' alg.par_idx '.txt'];
  alg.eqv_file = ['eqvars/eqvars' alg.par_idx '.txt'];
  alg.moments_file = ['logs/moments' alg.par_idx '.txt'];
  alg.targ_file = ['targets/target' num2str(alg.targ_idx) '_m' num2str(alg.m_cutoff) '.txt'];
  alg.wgtvec_file = ['targets/wgtvec' num2str(alg.targ_idx) '_m' alg.wgt_idx '.txt'];
  alg.covmat_file = ['targets/covmat' num2str(alg.targ_idx) '_m' num2str(alg.m_cutoff) '.txt'];
  alg.wgtmat_file = ['targets/wgtmat' num2str(alg.targ_idx) '_m' alg.wgt_idx '.txt'];
  alg.gradmat_file = ['logs/gradmat' alg.par_idx '.txt'];

  % create directories
  if (~isdir('logs'))
    mkdir('logs');
  end
  if (~isdir('output'))
    mkdir('output');
  end
  if (~isdir('evals'))
    mkdir('evals');
  end
  if (~isdir('policy'))
    mkdir('policy');
  end
  if (~isdir('latex'))
    mkdir('latex');
  end

end
