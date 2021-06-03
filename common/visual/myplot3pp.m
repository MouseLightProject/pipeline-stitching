function myplot3pp(varargin)
    if ishghandle(varargin{1}) ,
        ax = varargin{1} ;
        pts = varargin{2} ;
        rest_of_args = varargin(3:end) ;
        plot3(ax, pts(:,1), pts(:,2), pts(:,3), rest_of_args{:}) ;
    else
        pts = varargin{1} ;
        rest_of_args = varargin(2:end) ;
        plot3(pts(:,1), pts(:,2), pts(:,3), rest_of_args{:}) ;
    end
end
