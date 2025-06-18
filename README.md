## Benchmarks

Cf. folder `examples`.<br>
For our tool, run the `.jl` file.
For `LoopInvGen`, use the `.sl` file.
For `Interproc`, use the `.ip` file.<br>
**N.B.:** For `LoopInvGen` and `Interproc`, when the base case fails, we sometimes consider simpler versions (mainly to get a confidence that the program description is correct).
We occasionally also consider harder versions.

| Case | Our | `LoopInvGen` | `Interproc` |
| --- | --- | --- | --- |
| rotation | very fast | N/A (reals) | very fast |
| dai2020_SecIV_A | very fast | N/A (reals) | very fast
| LIA_simple1D | very fast | very fast | to do
| LIA_askew2D | very fast | very fast | to do
| LIA_consensus2D | fast | timeout | timeout
||| simpler1: very fast | simpler: very fast
||| simpler2: sometimes fast
| liu2022_Fig1 | fast | very fast | very fast |
||| hard: very fast | hard: very fast |


Instructions for `LoopInvGen`: https://github.com/SaswatPadhi/LoopInvGen<br>
Instructions for `Interproc`: https://github.com/Edivad99/interproc-docker

## Illustrations

Method<br>
![](https://github.com/guberger/CEGISPolyhedralBarrier.jl/blob/main/animation_rotating.gif)

Invariant for "examples/rotation.jl"<br>
<img src="https://github.com/guberger/CEGISPolyhedralBarrier.jl/blob/main/fig_rotation_full.png" width="600">

Invariant for "examples/liu2022_Fig1.jl"<br>
<img src="https://github.com/guberger/CEGISPolyhedralBarrier.jl/blob/main/fig_liu2022_Fig1.png" width="600">
