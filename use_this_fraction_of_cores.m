function use_this_fraction_of_cores(fraction)
    physical_core_count = feature('numcores') ;
    maximum_core_count_desired = round(fraction * physical_core_count) ;
    poolobj = gcp('nocreate');  % If no pool, do not create new one.
    if isempty(poolobj) ,
        parpool([1 maximum_core_count_desired]) ;
    end
    poolobj = gcp('nocreate');  % If no pool, do not create new one.
    core_count = poolobj.NumWorkers ;
    fprintf('Using %d cores.\n', core_count) ;
end
