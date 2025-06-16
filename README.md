## Illustration

![GUI](https://github.com/guberger/CEGISPolyhedralBarrier.jl/blob/main/animation_rotating.gif)

## Benchmarks

"rotation": very fast (LoopInvGen: N/A because reals)<br>
"dai2020_SecIV_A": very fast (LoopInvGen: N/A because reals)<br>
"LIA_simple1D": very fast (LoopInvGen: very fast)<br>
"LIA_askew2D": very fast (LoopInvGen: very fast)<br>
"LIA_consensus2D": medium fast (LoopInvGen: easy->fast, medium->slow, hard->timeout)<br>
"liu2022_Fig1": medium fast (LoopInvGen: timeout)<br>
