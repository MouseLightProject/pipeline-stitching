function make_videos_2019_03_26() 
    
    inputfolder = '2019-03-26'

    brain = inputfolder;
    inputfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/acquisition/%s',brain);
    pipelineoutputfolder = sprintf('/nrs/mouselight/pipeline_output/%s',brain);
    arch = lower(computer('arch'));
    if arch(1:2) == 'pc'
        error('windows machine, set the input using input arguments')
    else
        experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s-%s',brain,getenv('USER'));
    end


    matfolder = fullfile(experimentfolder,'matfiles/');

    unix('umask u=rwx,g=rwx,o=rx');

    scopefile = fullfile(matfolder,'scopeloc.mat');


    fprintf('Running descriptorMatchQuality stage...\n') ;
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
    load(fullfile(matfolder,'scopeparams_pertile'),'scopeparams')
    load(fullfile(matfolder,'regpts'),'regpts')
    mkdir('./videos')
    videofile = sprintf('./videos/%s-1stiter-ch1-%s',brain,date())
    descriptorMatchQuality(regpts,scopeparams{end},scopeloc,videofile)
end
