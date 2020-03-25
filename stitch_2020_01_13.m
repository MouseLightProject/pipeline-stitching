function stitch_2020_01_13()  
    addpath(genpath('/groups/mousebrainmicro/mousebrainmicro/Software/pipeline-stitching-shared/'))
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2020-01-13/Tiling' %#ok<NOPRT>
    pipeline_output_folder = '/nrs/mouselight/pipeline_output/2020-01-13' %#ok<NOPRT>
    %stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2020-01-13/stitching-output_TF'  %#ok<NOPRT>
    stitching_output_folder_path = '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2020-01-13/stitching-output'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
