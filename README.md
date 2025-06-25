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

| Experiment                  | # Var Int | # Var Real | Our Time | Our # Iter | LoopInvGen Time | Interproc Time |
|-----------------------------|-----------|------------|----------|------------|-----------------|----------------|
| `LIA_incr1D`                | 1         | 0          | 1 sec    | 1          | 1 sec           | 1 sec          |
| `LIA_incr2D`                | 2         | 0          | 1 sec    | 1          | 1 sec           | 1 sec          |
| `LIA_consensus`             | 4         | 0          | 18 sec   | 100        | T/O             | T/O            |
| `\|- simpler1 or 2`          | 4         | 0          | 12 sec   | 101        | ~3 sec          | 1 sec          |
| `LIA_liu20222_fig1`         | 2         | 0          | 5 sec*   | 104        | 1 sec           | 1 sec          |
| `\|- simpler`                | 2         | 0          | 8 sec    | 154        | 1 sec           | 1 sec          |
| `LRA_rotation`             | 0         | 2          | 2 sec    | 24         | N/A             | 1 sec          |
| `LRA_dai2020_sec4a`        | 0         | 2          | 1 sec    | 4          | N/A             | 1 sec          |
| `LRA_roux2015_fig3`        | 0         | 4          | 550 sec  | 787        | N/A             | T/O            |
| `LRA_roux2015_fig4`        | 0         | 4          | 1093 sec | 802        | N/A             | T/O            |
| `LRIA_heater`              | 2         | 2          | 21 sec   | 160        | N/A             | 1 sec          |
| `\|- LIA_heater`            | 4         | 0          | 9 sec    | 157        | T/O             | 1 sec          |
| ` \|- simpler`             | 4         | 0          | 6 sec    | 162        | 1 sec           | 1 sec          |
| `LIA_train_speed` (n=3)    | 1         | 4          | 5 sec    | 94         | N/A             | ??             |
| `LIA_train_speed` (n=4)    | 1         | 5          | 51 sec   | 173        | N/A             | T/O            |
| `LIA_train_speed` (n=5)    | 1         | 6          | 10 hour  | 956        | N/A             | T/O            |
| `LIA_train_distance` (n=2) | 0         | 4          | 2 sec    | 21         | N/A             | ??             |
| `LIA_train_distance` (n=3) | 0         | 6          | 600 sec  | 260        | N/A             | ??             |

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
