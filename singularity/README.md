You can build a [tracy](https://github.com/gear-genomics/tracy) singularity container (SIF file) using

`sudo singularity build tracy.sif tracy.def`

Once you have built the container you can run analysis using

`singularity exec tracy.sif tracy basecall input.ab1`
