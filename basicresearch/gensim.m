function [igen,egen,pgen] = gensim(sprob,cgeom,cond)
  global alg p

  for m=1:p.M
    am = sprob*m/(sprob*m+p.M);
    sm = 1 - am;

    % accumulate
    in = 0;
    en = 0;
    pn = 0;
    igensum = 0;
    egensum = 0;
    pgensum = 0;

    % simulate
    for j=1:alg.citK
      k = 1 + geornd(sm); % number of spillovers
      mv = randi([1,p.M],k,1); % particular industry draw
      nv = geornd(cgeom,k,1); % number of cites per spillover
      nb = zeros(1,p.M); % number of cites per industry
      for i=1:k
        mi = mv(i);
        nb(mi) = nb(mi) + nv(i);
      end

      mb = nb > 0;
      tn = sum(nb);
      tm = sum(mb);

      % whether to condition on m>1
      if (cond)
        inc = (tm > 1) && (tn > 0);
      else
        inc = true;
      end

      if (inc)
        if (tn > 0)
          shr = nb/tn;
          ig = 1 - sum(shr.^2);
        else
          ig = 0.0;
        end
        in = in + 1;
        igensum = igensum + ig;

        eg = tm;
        en = en + 1;
        egensum = egensum + eg;

        pg = tm > 1;
        pn = pn + 1;
        pgensum = pgensum + pg;
      end
    end

    if (cond)
      enull = 1.0;
      pnull = 1.0;
    else
      enull = 0.0;
      pnull = 0.0;
    end

    % return average
    if (in == 0)
      igen(m) = 0.0;
    else
      igen(m) = igensum/in;
    end
    if (en == 0)
      egen(m) = enull;
    else
      egen(m) = egensum/en;
    end
    if (pn == 0)
      pgen(m) = pnull;
    else
      pgen(m) = pgensum/pn;
    end
  end
end
