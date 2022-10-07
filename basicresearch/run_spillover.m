function run_objective
    global alg p eq

    alg = {};
    initalg();
    [par0,names] = parse_params(alg.par_file);

    fid = fopen('output/identify_spillover.csv', 'w+');
    fprintf(fid,'pb,epb1,epb2,epb3,epb4,epb5,epb6,epb7,epb8,bri1,bri2,bri3,bri4,bri5,bri6,bri7,bri8,\n');

    for pb=[0.0,0.05,0.1,0.15]
        params = par0;
        params(4) = pb;

        fprintf(1,'pb = %f\n',pb);
        [alg.eqv0,err] = eqstand(params);
        [score,moments,data] = compfs();

        fmt = strcat(strjoin(repmat({'%f'},1,17),','),'\n');
        fprintf(fid,fmt,pb,moments(1:16));
    end

    fclose(fid);

end
