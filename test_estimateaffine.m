load estimateaffine_inputs_2019_10_31.mat

[scopeparams_good,validthis_good] = util.estimateaffine_good(paireddescriptor,neighbors,scopeloc,params,curvemodel,old) ;

%[scopeparams_bad,validthis_bad] = util.estimateaffine_bad(paireddescriptor,neighbors,scopeloc,params,curvemodel,old) ;

[scopeparams,validthis] = util.estimateaffine(paireddescriptor,neighbors,scopeloc,params,curvemodel,old)

scopeparams_check = isequal(scopeparams, scopeparams_good) 

validthis_check = isequal(validthis, validthis_good)
