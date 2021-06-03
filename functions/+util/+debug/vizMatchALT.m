function [Ic_m,It_m] = vizMatchALT(scopeloc,neigs,descriptors,ineig,pixshift,iadj,X,Y,X_,Y_)
    % scopeloc: Info about all the raw tiles, including file paths and location in the lattice
    % neigs: soemthing x 4 double array, first col is tile index, 2nd col is x+1 tile, 3rd col is y+1 col,
    %        4th col is z+1 col
    % descriptors: 1 x tile_count cell array, each element containing a point_count x 5 double array.  1st three cols
    %              are integer xyz voxel coords of fiducial points (not sure if zero- or one-based), 4th col is DoG
    %              filter value at the point, 5th col is the foreground p-map value at the point.
    % ineig: index of the row in neigs we're going to use to get the central and adjacent tile index.
    %        Seems like the central til index is almost always equal to this, but I guess maybe not always?
    % imsize_um: 1 x 3, xyz dimensions of each tile, computed in the (not optimal IMHO) (n-1)*spacing way
    % idj: index of dimension along which the tile differ.  1==x, 2==y, 3==z.
    % X: something x 3, xyz coords of fiducials in central tile that are likely to be near the overlap
    %    region between the two tiles.
    % Y: something x 3, xyz coords of fiducials in other tile that are likely to be near the overlap
    %    region between the two tiles.
    % X_: match_count x 3, xyz coords of fiducials in central tile for which matches have been found
    %     with other tile.
    % Y_: match_count x 3, xyz coords of fiducials in other tile for which matches have been found
    %     with central tile.  (I assume ordered to match up with points in X_.)
    
    persistent fig
    
    idxcent = neigs(ineig,1);
    idxadj = neigs(ineig,iadj+1);
    
    channel_index = 1 ;
    Ic = util.getTilefromId(scopeloc,idxcent, channel_index);
    It = util.getTilefromId(scopeloc,idxadj, channel_index);
    Ic_m = rot90(max(Ic,[],3),2);
    It_m = rot90(max(It,[],3),2);
    dims = size(Ic);dims=dims([2 1 3]);
    %%
    %idxcent = neigs(ineig,1);
    descent = descriptors{idxcent};
    descent = double(descent(:,1:3));
    descent = util.correctTiles(descent,dims); % flip dimensions
    
    descadj = descriptors{idxadj};
    descadj = double(descadj(:,1:3)); % descadj has x-y-z-w1-w2 format
    descadj = util.correctTiles(descadj,dims); % flip dimensions
    
    %stgshift = 1000*(scopeloc.loc(idxadj,:)-scopeloc.loc(idxcent,:));
    %pixshift = round(stgshift.*(dims-1)./(imsize_um));
    descadj = descadj + ones(size(descadj,1),1)*pixshift; % shift with initial guess based on stage coordinate
    
    %%
    if isempty(fig) || ~isvalid(fig),
        fig = figure('Color', 'w', 'Name', mfilename()) ;
    end
    clf(fig) ;
    ax = axes(fig) ;
    %figure(43) ;
    %cla(ax) ;
    RA = imref2d(size(Ic_m),[1 dims(1)],[1 dims(2)]);
    if iadj==1
        RB = imref2d(size(It_m),[1 dims(1)]+pixshift(1),[1 dims(2)]);
    else
        RB = imref2d(size(It_m),[1 dims(1)],[1 dims(2)]+pixshift(2));
    end
    imshowpair(imadjust(Ic_m),RA,imadjust(It_m),RB,'Parent',ax,'falsecolor','Scaling','joint','ColorChannels','green-magenta')
    hold(ax,'on') ;
    myplot3pp(ax, descent-1, 'bo', 'MarkerSize', 6, 'LineWidth', 1) ;
    myplot3pp(ax, descadj-1, 'yo', 'MarkerSize', 6, 'LineWidth', 1) ;
    myplot3pp(ax, X-1, 'bo', 'MarkerSize', 12, 'LineWidth', 1) ;
    myplot3pp(ax, Y-1, 'yo', 'MarkerSize', 12, 'LineWidth', 1) ;
    Y_2 = Y_;
    Y_2(:,iadj) = Y_2(:,iadj) + pixshift(iadj);
    XX = [X_(:,1),Y_2(:,1),nan*X_(:,1)]'-1;
    YY = [X_(:,2),Y_2(:,2),nan*X_(:,2)]'-1;
    plot(ax, XX, YY, 'r') ;
    hold(ax,'off') ;
    drawnow
end
