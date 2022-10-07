function run_academic
    global eq

    fid = fopen('output/academic_gdp_fixw.csv', 'w+');
    fprintf(fid,'ac_gdp,taua,Lapp,wage,hot\n');

    for ac=0.0:0.001:0.02
        welfare(nan,nan,ac,nan,nan);
        fprintf(fid,'%f,%f,%f,%f,%f\n',ac,eq.taua,sum(eq.ca.*eq.phivec),eq.w,eq.hot);
    end

    fclose(fid);

end
