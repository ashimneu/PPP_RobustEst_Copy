function saveoutlier(database_path,outliervec,outlierbin,outlierparam)
    outlierdb.vec   = outliervec;
    outlierdb.bin   = outlierbin;
    outlierdb.param = outlierparam;    
    save(database_path,'outlierdb');    
end

