function stitch_2020_02_24_round2()  
    addpath(genpath('/groups/mousebrainmicro/mousebrainmicro/Software/pipeline-stitching-shared/'))
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2020-02-24/Tiling' %#ok<NOPRT>
    pipeline_output_folder = '/nrs/mouselight/pipeline_output/2020-02-24' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2020-02-24/stitching-output_TF3'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
