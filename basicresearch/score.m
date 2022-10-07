function sc = score(moments)
  global alg p eq

  err = alg.targs-moments;

  % blacklist - no hope on variance of profitability
  for b=alg.blist
    err(b) = 0.0;
  end

  nerr = err./alg.targs;
  sc = nerr*nerr';

  if (isnan(sc))
    sc = 1000.0;
  end
  if (isinf(sc))
    sc = 1000.0;
  end
  
  bind = setdiff(1:alg.ntargs,alg.blist);
  reltarg = eq.moments(bind);
  
  if (alg.mfid ~= -1)
    fprintf(alg.mfid,'%20.15f ',reltarg); fprintf(alg.mfid,'\n');
  end
end

