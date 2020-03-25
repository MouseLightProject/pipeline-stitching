function stitch_2019_12_25()  
    addpath(genpath('/groups/mousebrainmicro/mousebrainmicro/Software/pipeline-stitching-shared/'))
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2019-12-25/Tiling' %#ok<NOPRT>
    pipeline_output_folder = '/nrs/mouselight/pipeline_output/2019-12-25' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-12-25/stitching-output_TF2'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
