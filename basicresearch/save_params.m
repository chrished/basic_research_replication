function save_params(fname,pvec,pnames)

  np = length(pvec);

  fpid = fopen(fname,'w');
  for i=1:np
    fprintf(fpid,'%20s : %18.15f\n',pnames{i},pvec(i));
  end

  fclose(fpid);

end

