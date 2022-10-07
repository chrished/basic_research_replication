function run_objective
    global alg p eq

    alg = {};
    initalg();
    [par0,names] = parse_params(alg.par_file);

    fid = fopen('output/objective_applied_spill.csv', 'w+');
    fprintf(fid,'pa,taua,w,obj\n');

    for pa=0.0:0.005:0.08
        params = par0;
        params(3) = pa;

        [alg.eqv0,err] = eqstand(params);
        [score,moments] = compfs();
        fprintf(1,'%f -> %f\n',pa,score);

        fprintf(fid,'%f,%f,%f,%f\n',pa,eq.taua,eq.w,score);
    end

    fclose(fid);

end
