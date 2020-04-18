function result  = is_scattered_interpolant_degenerate(interpolant)
    lastwarn('') ;  % clear the lastwarn state
    interpolant.Points ;  % this will provoke the warning if the triangulation is degenerate
    [~, id] = lastwarn() ;
    result = isequal(id, 'MATLAB:scatteredInterpolant:InterpEmptyTri3DWarnId') ;
end
