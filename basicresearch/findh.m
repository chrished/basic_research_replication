function findh
  global p alg eq

  eq.denomh = eq.r-eq.g+eq.tau;

  for m=(p.M-(0:(p.M-1)))
    alg.m = m;

    try
      eq.hvec(m) = fzero(@hfunc,[0.0 alg.hmax]);
    catch
      try
        alg.rlim = 0;
        eq.hvec(m) = fminbnd(@(h)(-hfunc(h)),0.0,alg.hmax);
        alg.rlim = 1;
      catch
        fprintf(1,'hvec failed\n');
        eq.hvec(m) = 0.0;
        alg.hfail = 1;

        %{}
        hvec = 0.0:0.001:alg.hmax;
        for i=1:length(hvec)
          hval(i) = hfunc(hvec(i));
        end
        plot(hvec,hval); grid on;
        pause;
        %}
      end

    end

    %fprintf(alg.fid,'hvec(%i) = %15.10f\n',m,eq.hvec(m));
  end
end
