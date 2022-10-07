function makeweights
  global alg

  alg = {};
  cd ..
  initalg();
  cd targets

  wgtmat_file = ['../' alg.wgtmat_file];
  covmat_file = ['../' alg.covmat_file];
  wgtvec_file = ['../' alg.wgtvec_file];
  target_file = ['../' alg.targ_file];

  % find moments used
  wgtvec = load(wgtvec_file);
  wlist = wgtvec > 0.0;

  % load targets
  mvec = load(target_file);
  mvec_wl = mvec(wlist);

  % load covmat
  covmat = load(covmat_file);
  n_tot = length(covmat);

  % exit rate
  %covmat(20,:) = 0.0;
  %covmat(:,20) = 0.0;
  %covmat(20,20) = 0.0001766^2;

  % blasted R&D/sales
  %covmat(21,:) = 0.0;
  %covmat(:,21) = 0.0;
  %covmat(21,21) = (0.05*mvec_wl(21)).^2;

  % other missing covs (22-pubpriv,24-agggrowth,25-spilldiff)
  covmat(22,22) = (0.02*mvec_wl(22))^2; % priv/pub basic
  covmat(24,24) = (0.02*mvec_wl(24))^2; % agg growth
  covmat(25,25) = (0.02*mvec_wl(25))^2; % spillover diff

  % invert to find weighting matrix
  wgtmat = inv(covmat);

  % kill R&D/sales
  wgtmat(21,:) = 0.0;
  wgtmat(:,21) = 0.0;

  % double up basic research moments
  %for i=1:16
  %  wgtmat(i,i) = 3.0*wgtmat(i,i);
  %end
  wgtmat(1:16,:) = 2*wgtmat(1:16,:);
  wgtmat(17:end,1:16) = 2*wgtmat(17:end,1:16);

  % save to files
  save(wgtmat_file,'wgtmat','-ascii','-double');

  % find the contribution of each element
  cont = mvec_wl.*diag(wgtmat).*mvec_wl/10000;
  disp([1:n_tot; cont']');

end
