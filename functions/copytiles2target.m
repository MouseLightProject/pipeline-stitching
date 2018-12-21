function tilestruct = copytiles2target(targetfolder,scopeloc,inds)
% copies raw files to a target location. useful for cropped renders and
% comparison to other stitching software
mkdir(targetfolder)
ntiles = length(inds);
tilestruct = [];
for ix = 1:ntiles
    it = inds(ix);
    [folderpath,tilepath] = fileparts(scopeloc.filepath{it});
    mkdir(fullfile(targetfolder,scopeloc.relativepaths{it}))
    unix(sprintf('cp -a %s/* %s',folderpath,fullfile(targetfolder,scopeloc.relativepaths{it})));
    
    % targetstructure
    scopeacqparams = readScopeFile(fileparts(scopeloc.filepath{it}));
    tilestruct(ix).scopeacqparams = scopeacqparams;
    tilestruct(ix).relativepaths = scopeloc.relativepaths{it};
    tilestruct(ix).loc = scopeloc.loc(it,:);
    tilestruct(ix).gridix = scopeloc.gridix(it,:);
end
