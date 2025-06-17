## Benchmarks

Cf. folder "examples".
For our tool, just run the `.jl` file.
For `LoopInvGen`, use the `.sl` file.

"rotation": very fast (`LoopInvGen`: N/A because reals)<br>
"dai2020_SecIV_A": very fast (`LoopInvGen`: N/A because reals)<br>
"LIA_simple1D": very fast (`LoopInvGen`: very fast)<br>
"LIA_askew2D": very fast (`LoopInvGen`: very fast)<br>
"LIA_consensus2D": medium fast (`LoopInvGen`: easy->fast, medium->slow, hard->timeout)<br>
"liu2022_Fig1": medium fast (`LoopInvGen`: timeout)<br>

Instructions for `LoopInvGen`: https://github.com/SaswatPadhi/LoopInvGen
Instructions for `Interproc`: https://github.com/Edivad99/interproc-docker

## Illustrations

Method<br>
![GUI](https://github.com/guberger/CEGISPolyhedralBarrier.jl/blob/main/animation_rotating.gif)

Invariant for "examples/rotation.jl"<br>
![GUI](https://github.com/guberger/CEGISPolyhedralBarrier.jl/blob/main/fig_rotation_full.png)

Invariant for "examples/liu2022_Fig1.jl"<br>
![GUI](https://github.com/guberger/CEGISPolyhedralBarrier.jl/blob/main/fig_liu2022_Fig1.png)
