function [fixed_point_xyz_from_match_index, moving_point_xyz_from_match_index] = ...
        descriptorMatch4XY(xyz_from_fixed_point_index, xyz_from_moving_point_index, params)
    % Descriptor match for XY directions.  Uses tighter constraints as there is
    % no physical cut and initialization can be reliably estimated from beads

    % "fp" == fixed_point
    % "mp" == moving_point
    % "fpi" = fixed_point_index
    % "mpi" = moving_point_index
    % "vfp" = viable fixed point
    % "vmp" = viable moving point
    % "vfpi" = viable fixed point index
    % "vmpi" = viable moving point index
    
    % Get some dimensions
    moving_point_count = size(xyz_from_moving_point_index,1) ;  %#ok<NASGU>
    fixed_point_count = size(xyz_from_fixed_point_index,1) ;  %#ok<NASGU>
    
    % Set parameters
    cpd_options = params.opt ;
    initial_distance_threshold = params.initial_distance_threshold ;
    registered_distance_threshold = params.registered_distance_threshold ;
    kernel_bandwith_cutoff = 20 ;
    
    % Eliminate fixed/moving landamrks that are too far from the nearest
    % moving/fixed landmark.  The remaining ones are "viable".
    distance_from_fixed_point_index_from_moving_point_index = pdist2(xyz_from_fixed_point_index, xyz_from_moving_point_index) ;
      % fixed_point_count x moving_point_count
    distance_to_nearest_fixed_point_from_moving_point_index = min(distance_from_fixed_point_index_from_moving_point_index, [], 1) ;
      % 1 x moving_point_count
    distance_to_nearest_moving_point_from_fixed_point_index = min(distance_from_fixed_point_index_from_moving_point_index, [], 2) ;
      % fixed_point_count x 1
    is_near_moving_point_from_fixed_point_index = distance_to_nearest_moving_point_from_fixed_point_index<initial_distance_threshold ;
    xyz_from_viable_fixed_point_index = xyz_from_fixed_point_index(is_near_moving_point_from_fixed_point_index,:) ;
    is_near_fixed_point_from_moving_point_index = distance_to_nearest_fixed_point_from_moving_point_index<initial_distance_threshold ;
    xyz_from_viable_moving_point_index = xyz_from_moving_point_index(is_near_fixed_point_from_moving_point_index,:) ;
    
    % Count how many are viable
    viable_fixed_point_count = size(xyz_from_viable_fixed_point_index, 1) ;
    viable_moving_point_count = size(xyz_from_viable_moving_point_index, 1) ;
    
    % If not enough viable points, fail out now
    if viable_fixed_point_count<3 || viable_moving_point_count<3 ,
        fixed_point_xyz_from_match_index = zeros(0,3) ;
        moving_point_xyz_from_match_index = zeros(0,3) ;
        return
    end

    % Run CPD to register points in the moving image
    diffs = distance_from_fixed_point_index_from_moving_point_index(distance_from_fixed_point_index_from_moving_point_index<kernel_bandwith_cutoff) ;  
    std_data = min(kernel_bandwith_cutoff,max(2,sqrt(sum(diffs.^2)/(numel(diffs)-1)))) ; 
      % force kernel bandwidth to [2 10]. These are pixel diffences, very unlikely to have a matching feature more than cuttoff 
    h = std_data*(4/3/numel(diffs))^(1/5) ; % silvermans rule
    cpd_options.beta = h ;  % (default 2) Gaussian smoothing filter size. Forces rigidity.
    Transform = cpd_register(xyz_from_viable_fixed_point_index, xyz_from_viable_moving_point_index, cpd_options) ;
    registered_xyz_from_viable_moving_point_index = Transform.Y ;
    
    % Check if match is found
    registered_distance_from_vfpi_from_vmpi = pdist2(xyz_from_viable_fixed_point_index, registered_xyz_from_viable_moving_point_index) ;
      % viable_fixed_point_count x viable_moving_point_count
    [registered_distance_to_nearest_vfp_from_vmpi, best_match_vfpi_from_vmpi] = min(registered_distance_from_vfpi_from_vmpi, [], 1) ;
      % 1 x viable_moving_point_count
    [registered_distance_to_nearest_vmp_from_vfpi, best_match_vmpi_from_vfpi] = min(registered_distance_from_vfpi_from_vmpi, [], 2) ;  %#ok<ASGLU>
      % viable_fixed_point_count x 1
    do_both_agree_to_match_from_vmpi = (1:viable_moving_point_count)'==best_match_vmpi_from_vfpi(best_match_vfpi_from_vmpi) ;
    vmpi_from_putative_match_index = find(do_both_agree_to_match_from_vmpi) ;
    vfpi_from_putative_match_index = best_match_vfpi_from_vmpi(vmpi_from_putative_match_index)' ;

    % For a putative match to be a final match, have to be close enough together
    % when registered
    registered_distance_from_putative_match_index = registered_distance_to_nearest_vfp_from_vmpi(vmpi_from_putative_match_index)' ;
    is_close_enough_from_putative_match_index = registered_distance_from_putative_match_index < registered_distance_threshold ;
    vfpi_from_match_index = vfpi_from_putative_match_index(is_close_enough_from_putative_match_index) ;
    vmpi_from_match_index = vmpi_from_putative_match_index(is_close_enough_from_putative_match_index) ;   
    
    % Compute the final outputs
    fixed_point_xyz_from_match_index = xyz_from_viable_fixed_point_index(vfpi_from_match_index,:) ;
    moving_point_xyz_from_match_index = xyz_from_viable_moving_point_index(vmpi_from_match_index,:) ;
end
