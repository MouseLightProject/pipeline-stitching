function scope = getScopeCoordinates(stitching_output_folder_path, tile_folder_path, is_sample_post_2016_04_04)
%GETSCOPECOORDINATES Summary of this function goes here
%
% [OUTPUTARGS] = GETSCOPECOORDINATES(INPUTARGS) Explain usage here
%
% Inputs:
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2016/09/21 16:33:49 $	$Revision: 0.1 $
% Copyright: HHMI 2016

% Recurse through the input tiles folder, collecting the names of the
% .acquisition files.  Store those in scopeacquisitionlist.txt
args = struct() ;
args.level = 3;
args.ext = 'acquisition';
args.skip = {''};
args.keep = {''};
args.pattern = '\d';
opt.seqtemp = fullfile(stitching_output_folder_path, 'scopeacquisitionlist.txt') ;
opt.inputfolder = tile_folder_path;
% if exist(opt.seqtemp, 'file') == 2
%     % load file directly
% else
args.fid = fopen(opt.seqtemp,'w');
recdir(opt.inputfolder,args)
% end

% Now read the file we just wrote, and make an in-memory list of the .acquisition files
fid=fopen(opt.seqtemp,'r');
inputfiles = textscan(fid,'%s');
inputfiles = inputfiles{1};
fclose(fid);
input_file_count = length(inputfiles) ;

gridix = cell(1,input_file_count) ;
loc = cell(1,input_file_count) ;
if is_sample_post_2016_04_04 ,
    parfor_progress(input_file_count);
    for ifile = 1:input_file_count
        parfor_progress();
        scvals = util.scopeparser(inputfiles{ifile});
        gridix{ifile} = [scvals.x scvals.y scvals.z scvals.cut_count];
        loc{ifile} = [scvals.x_mm scvals.y_mm scvals.z_mm];
    end
    parfor_progress(0);
    raw_grids = cat(1,gridix{:});
    grids = shift_zs(raw_grids) ;  % sometimes there are z=0 or a z=-1 tiles, so we shift up to accomodate Matlab indexing
    locs = cat(1,loc{:});
else
    parfor_progress(input_file_count);
    parfor ifile = 1:input_file_count
        parfor_progress;
        scvals = scopeparser(inputfiles{ifile});
        loc{ifile} = [scvals.x_mm scvals.y_mm scvals.z_mm];
    end
    parfor_progress(0);
    locs = cat(1,loc{:});
    locs_ = cat(1,loc{:});
    mins = min(locs_);
    locs_ = locs_-ones(size(locs_,1),1)*mins;
    
    
    difss = diff(locs_);
    r(1) = median(difss(:,1));
    r(2) = median(difss(difss(:,2)>eps,2));
    clear x y
    [x,y] = deal(zeros(size(locs_,1),1));
    for ii=1:size(locs_)
        x(ii) = fix((locs_(ii,1)+r(1)/2)./r(1))+1;
        y(ii) = fix((locs_(ii,2)+r(2)/2)./r(2))+1;
    end
    ran = max([y x],[],1);
    oldind = 0;
    currz = 1;
    grids = zeros(size(locs_));
    for ii=1:size(locs_)
        currind = sub2ind(ran([2 1]),x(ii),y(ii));
        if currind<oldind
            % append
            currz = currz+1;
        end
        grids(ii,:) = [x(ii),y(ii),currz];
        oldind = currind;
    end
    
end

% Make a list of all the relative paths to the .acquisition files, relative
% to the input tile folder
relativepaths = cell(length(inputfiles),1);
for ii=1:length(inputfiles)
    relativepaths{ii} = fileparts(inputfiles{ii}(length(tile_folder_path)+1:end));
end

% Store everything in the scope structure for return to caller
scope = struct() ;
scope.gridix = grids ;
scope.loc = locs;
scope.filepath = inputfiles;
scope.relativepaths = relativepaths;
end


function result = shift_zs(grids)
    % sometimes there are z=0 or a z=-1 tiles, so we shift up to accomodate Matlab indexing
    zs = grids(:,3) ;
    min_z = min(zs) ;
    if min_z<1 ,
        new_zs = zs-min_z+1 ;
        result = grids ;
        result(:,3) = new_zs ;
    else
        result = grids ;
    end    
end
