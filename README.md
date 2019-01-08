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
const source="/nrs/mouselight/cluster/classifierOutputs/2018-06-14"  # path to directory holding tilebase.cache.yml
const output="/nrs/mouselight/SAMPLES/2018-06-14" # destination path of octree. Note the absence of trailing slash. Existing files will be overritten

const shared_scratch="/nrs/mouselight/scratch/<yourId>"
const logfile_scratch="/groups/mousebrainmicro/mousebrainmicro/scratch/<yourId>"  # should be on /groups
const delete_scratch="as-you-go"   # "never", "at-end", or "as-you-go"
```
Note that if the scratch paths are not unique, they will be wiped and overriten.

3. Choose interpolation and downsampling_function
For a 'regular' sample we typically use 2nd brightest. For registration, median, i.e., the 5th brighest non-zero of the 8 pixels. The way the parameters file is is set, we allow only one implementationm commenting out alternatives not in use. E.g., for 2nd brightest:

```julia
# 1. the simplest and fastest
#downsampling_function(arg::Array{UInt16,3}) = (@inbounds return arg[1,1,1])

# 2. equivalent to mean(arg) but 30x faster and half the memory
#downsampling_function(arg::Array{UInt16,3}) = UInt16(sum(arg)>>3)

# 3. 2nd brightest of the 8 pixels
# equivalent to sort(vec(arg))[7] but half the time and a third the memory usage
function downsampling_function(arg::Array{UInt16,3})
  m0::UInt16 = 0x0000
  m1::UInt16 = 0x0000
  for i = 1:8
    @inbounds tmp::UInt16 = arg[i]
    if tmp>m0
      m1=m0
      m0=tmp
    elseif tmp>m1
      m1=tmp
    end
  end
  m1
end

# 4. Nth brightest non-zero of the 8 pixels
#function downsampling_function(arg::Array{UInt16,3})
#  n=5
#  m = fill(0x0000,n)
#  for i = 1:8
#    @inbounds tmp::UInt16 = arg[i]
#    for i=1:n
#      if tmp>m[i]
#        m[i+1:n]=m[i:n-1]
#        m[i]=tmp
#        break
#      end
#    end
#  end
#  for i=n:-1:1
#    m[i]==0 || return m[i]
#  end
#  return 0x0000
#end
```

4. Start the render as follows:
```sh
ssh mluser@login1
cd /groups/mousebrainmicro/mousebrainmicro/Software/barycentric#/src/render/src
./render <path_to_render_parameters.jl>
```

### Rendering Notes
- As of January 2019, one needs to blacklist some mis-behaving nodes:
```julia
const bad_nodes = ["e10u13"]  # e.g. ["h09u20"]
```

- For an initial look at the data, one can speed up rendering times using:
```julia
const voxelsize_um=[1., 1, 1]  # desired pixel size
const interpolation = "nearest"  # "nearest" or "linear"
downsampling_function(arg::Array{UInt16,3}) = (@inbounds return arg[1,1,1])
```
- For _daytime_ rendering and minimize the workload on the nrs filesystem, reduce to half the maximum number of jobs used to render leafs and the maximum number of nodes used to downsample the octree
```julia
const throttle_leaf_njobs = 64  # maximum number of jobs to use to render leafs
const throttle_octree_njobs = 128  # maximum number of compute nodes to use to downsample octree
```

- Render documentation resides at `/groups/mousebrainmicro/mousebrainmicro/Software/barycentric#/src/render/README.md`

- For triggering the render remotely (e.g., over VPN) under mluser, one can first ssh to a local machine (e.g., vega) as oneself, so that the login is performed using the keys in the .ssh directory of one's home directory

