function academic
  global alg p eq

  % academic innovation
  eq.ufund = max(0.0,p.That/(eq.w)-p.massac*p.acfixed)/p.massac;
  eq.xu = (eq.ufund/p.usigma_ns)^(1.0/p.ugamma);

  eq.freebase = (1.0+p.pb)*alg.bdole*p.massac*eq.xu;
  eq.acspill = (1.0+p.pb)*(1.0-alg.bdole)*p.massac*eq.xu;

  if (alg.fid ~= -1)
    fprintf(alg.fid,'xu = %7.5f, acspill = %7.5f, freebase = %7.5f\n',eq.xu,eq.freebase,eq.acspill);
  end

end
