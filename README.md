```julia
> ] activate ./examples
> ] instantiate
> include("./examples/name_of_example.jl")
```
## Benchmarks

Cf. folder `examples`.<br>
For our tool, run the `.jl` file.
For `LoopInvGen`, use the `.sl` file.
For `Interproc`, use the `.ip` file.<br>
**N.B.:** For `LoopInvGen` and `Interproc`, when the base case fails, we sometimes consider simpler versions (mainly to get a confidence that the program description is correct).
We occasionally also consider harder versions.
$n$ vars = (number of real variables, number of integer variables, number of locations)

| Case | $n$ vars | Our | `LoopInvGen` | `Interproc` |
| --- | --- | --- | --- | --- |
| rotation | (2, 0, 1) | very fast | N/A (reals) | very fast |
| dai2020_SecIV_A | (2, 0, 1) | very fast | N/A (reals) | very fast
| LIA_incr1D | (0, 1, 1) | very fast | very fast | to do
| LIA_incr2D | (0, 2, 1) | very fast | very fast | to do
| LIA_consensus2D | (0, 3, 2) | fast | timeout | timeout
|| same || simpler1: very fast |
|| same || simpler2: very fast |
|| same ||| simpler3: very fast |
| LIA_liu2022_Fig1 | (0, 2, 1) | fast (failed) | very fast | very fast |
|| same | simpler : fast |||
| thermostat | (2, 0, 4) | fast | N/A (reals) | very fast |
|| (0, 2, 4) || simpler1: very fast |
|| (0, 2, 4) || simpler2: timeout |
| roux2015_Eq3 | (4, 0, 1) | ~5min | N/A (reals) | to do |
| roux2015_Eq4 | (4, 0, 1) | ~5min | N/A (reals) | to do |
| train_speed | (4, 0, 2) | ~30sec | N/A (reals) | timeout |
|| (5, 0, 2) | ~10min | N/A (reals) | timeout |
|| (6, 0, 2) | ~10h | N/A (reals) | timeout |
|| (3, 1, 2) ||| 3-simpler: timeout |
| train_distance | (4, 0, 1) | ~10sec | N/A (reals) | to do |
|| (6, 0, 1) | ~5min | N/A (reals) | to do|

Instructions for `LoopInvGen`: https://github.com/SaswatPadhi/LoopInvGen<br>
Instructions for `Interproc`: https://github.com/Edivad99/interproc-docker

## Illustrations

Method<br>
![](https://github.com/guberger/CEGISPolyhedralBarrier.jl/blob/main/animation_rotating.gif)

Invariant for "examples/rotation.jl" (red=unsafe region, yellow=invariant, dots and lines=sample points and trajectories)<br>
<img src="https://github.com/guberger/CEGISPolyhedralBarrier.jl/blob/main/fig_rotation_full.png" width="600">

Invariant for "examples/liu2022_Fig1.jl" (red=unsafe region, yellow=invariant, dots and lines=sample points and trajectories)<br>
<img src="https://github.com/guberger/CEGISPolyhedralBarrier.jl/blob/main/fig_liu2022_Fig1.png" width="600">

## TODOs

- [x] Separate safe with margin (easy: just compute distance of safe face with
  inside and image points).
- [x] Maximize margin outside for contraction verifier (then: is counterexample only
  if margin > 0).
- [ ] Do not separate points that are already excluded by previous separator
  i.e., remove redundant separators (only when reset).
