function loadparams(params)
  global alg p

  % policies
  if (~isfield(p,'asubs'))
    p.asubs = 0.10;
  end
  if (~isfield(p,'bsubs'))
    p.bsubs = 0.10;
  end

  % fixed
  p.M = 10;

  % loaded
  p.crra = params(1);
  p.epsn = params(2);
  p.pa = params(3);
  p.pb = params(4);
  p.eta = params(5);
  p.alpha = params(6);
  p.massout = params(7);
  p.massac = params(8);
  p.agamma = params(9)+1.0;
  p.bgamma = params(10)+1.0;
  p.ugamma = params(10)+1.0;
  p.asigma_ns = params(11);
  p.asigma = params(11)*(1.0-p.asubs);
  p.bsigma_ns = params(12);
  p.bsigma = params(12)*(1.0-p.bsubs);
  p.usigma_ns = params(12);
  p.usigma = params(12);
  p.kappa = params(13);
  p.dmu = -params(14);
  p.dsig = params(15);
  p.zeta = params(16);
  p.qmin = params(17);
  p.disc = params(18);
  p.merge = params(19);
  p.ufrac = params(20);
  p.citx = params(21);

  if (alg.noexp == 1)
    p.merge = 0.0;
  end

  if (~isnan(alg.acmult))
    p.usigma_ns = p.usigma_ns/alg.acmult;
    p.usigma = p.usigma/alg.acmult;
  end

  % auxiliaries
  mvec = 1:p.M;
  p.rhoa = p.pa*mvec./p.M;
  p.rhob = p.pb*mvec./p.M;
  p.eprob = ((p.M-1):-1:0)/p.M;
  p.epsfrac = (p.epsn-1.0)/p.epsn;
  p.qminp = p.qmin^(p.epsn-1.0);
  p.bfcoeff = 1.0-p.bsubs;
  p.acfixed = exp(p.dmu+0.5*p.dsig^2);

  % academic stuff
  if (alg.noacad)
    p.That = 0.0;
  else
    if (isfield(p,'polzout'))
      p.That = p.polzout;
    else
      p.That = p.ufrac;
    end
    if (isfield(p,'polmac'))
      p.massac = p.polmac;
    end
  end

  % fully internalized basic
  if (alg.bintern)
    p.rhoa(:) = 0.99*p.pa;
    p.rhob(:) = 0.99*p.pb;
  end

  % param identifiers
  p.nparm = length(params);
  p.pnames = cell(p.nparm,1);
  p.pnames{1} = 'CRRA Utility Parameter';
  p.pnames{2} = 'Elasticity of Subsitution';
  p.pnames{3} = 'Cross-industry Spillover';
  p.pnames{4} = 'Basic Step Size';
  p.pnames{5} = 'Applied Step Size';
  p.pnames{6} = 'Mass of Entrants';
  p.pnames{7} = 'Mass of Academic Labs';
  p.pnames{8} = 'Applied Cost Curvature';
  p.pnames{9} = 'Basic Cost Curvature';
  p.pnames{10} = 'Applied Cost Scale';
  p.pnames{11} = 'Basic Cost Scale';
  p.pnames{12} = 'Exogenous Exit Rate';
  p.pnames{13} = 'Basic Fixed Mean';
  p.pnames{14} = 'Basic Fixed Stdev';
  p.pnames{15} = 'Product Cooldown Rate';
  p.pnames{16} = 'Minimum Quality Threshold';
  p.pnames{17} = 'Time Discount Rate';
  p.pnames{18} = 'Probability of Merger Offer';
  p.pnames{19} = 'Academic Fraction of GDP';

end
