# pipeline-stitching
computational pipeline that creates yml file that will be consumed by the render

1. Create a folder and clone the repo
```
cd pipe
clone https://github.com/erhanbas/pipeline-stitching.git 
```
2. run main.m in matlab
```
main(brain,[inputfolder,pipelineoutputfolder])
```
## inputs/outputs
    brain : sample name, e.g. '2018-06-14';
    inputfolder [OPTIONAL]: acqusition folder, e.g. sprintf('/groups/mousebrainmicro/mousebrainmicro/data/acquisition/%s',brain);
    pipelineoutputfolder [OPTIONAL]: outputfolder e.g. sprintf('/nrs/mouselight/cluster/sandbox2/%s',brain)

## Sample usage
```matlab
>> main('2018-06-14')
```
This will create tilebase.cache.yml file in the `/nrs/mouselight/cluster/classifierOutputs/2018-06-14` directory. Note that the render expects an hard-wired _tilebase.cache.yml_ filename.


## Sample render usage

1. Copy a sample parameter file next to yml file, e.g.
```sh
cp /groups/mousebrainmicro/mousebrainmicro/Software/barycentric4/src/render/parameters.jl /nrs/mouselight/cluster/classifierOutputs/2018-06-14/render_parameters.jl
```

2. Modify the paths in _render_parameters.jl_ as follows:
```julia
const user="<yourId>@hhmi.org"
const source="/nrs/mouselight/cluster/classifierOutputs/2018-06-14"  # path to tilebase.cache.yml
const output="/nrs/mouselight/SAMPLES/2018-06-14" # destination path of octree. Note the absence of trailing slash. Existing files will be overritten
```

3. Choose interpolation and downsampling_function

4. Start the render as follows:
```sh
ssh login1
cd /groups/mousebrainmicro/mousebrainmicro/Software/barycentric#/src/render
./render <path_to_render_parameters.jl>
```

### Rendering Notes
For an initial look at the data, one can speed up rendering times using:
```julia
const voxelsize_um=[1., 1, 1]  # desired pixel size
const interpolation = "nearest"  # "nearest" or "linear"
downsampling_function(arg::Array{UInt16,3}) = (@inbounds return arg[1,1,1])
```
For _daytime_ rendering and minimize the workload on the nrs filesystem, reduce to half the maximum number of jobs used to render leafs and the maximum number of nodes used to downsample the octree
```julia
const throttle_leaf_njobs = 64  # maximum number of jobs to use to render leafs
const throttle_octree_njobs = 128  # maximum number of compute nodes to use to downsample octree
```

Render documentation resides at `/groups/mousebrainmicro/mousebrainmicro/Software/barycentric#/src/render/README.md`
