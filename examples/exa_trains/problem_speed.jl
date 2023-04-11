# hot(i, N)'*x[1:N] = x[i]
hot(i, N) = [j == i ? 1.0 : 0.0 for j = 1:N]

function build_problem(Nt, α, Tstab, lim_lo, lim_up, xinit, xsafes)
    N = Nt + 1
    M = 2

    # System
    pieces = Piece{_PT,Matrix{Float64},Vector{Float64}}[]
    for cases in Iterators.product([(-1, 0, 1) for i = 1:N]...)
        afs = [AffForm([+ones(Nt); 0.0], 0.0), AffForm([-ones(Nt); 0.0], 0.0)]
        # afs = _AT[]
        A = Matrix{Float64}(I, N, N)
        b = zeros(N)
        i1 = i2 = 1
        for (i, case) in enumerate(cases) # dvcurr = vnext - vcurr
            if i < N
                icurr = mod(i - 1, Nt) + 1
                iprev = mod(i - 2, Nt) + 1
                if case == -1 # dvcurr ≥ lim_up/α --> vcurr += lim_up
                    push!(afs, AffForm(-hot(icurr, N), +lim_up/α))
                    b[icurr] -= lim_up
                    b[iprev] += lim_up
                elseif case == 1 # dvcurr ≤ lim_lo/α --> vcurr += lim_lo
                    push!(afs, AffForm(+hot(icurr, N), -lim_lo/α))
                    b[icurr] -= lim_lo
                    b[iprev] += lim_lo
                else # lim_lo/α ≤ dvcurr ≤ lim_up/α --> vcurr += dvcurr*α
                    push!(afs, AffForm(-hot(icurr, N), +lim_lo/α))
                    push!(afs, AffForm(+hot(icurr, N), -lim_up/α))
                    A[icurr, icurr] -= α
                    A[iprev, icurr] += α
                end
            else
                if case == -1 # Trans & t ≤ Tstab - 1
                    push!(afs, AffForm(+hot(N, N), -float(Tstab - 1)))
                    push!(afs, AffForm(-hot(N, N), 0.0))
                    b[N] = 1.0
                elseif case == 1 # Trans & t ≥ Tstab - 1
                    push!(afs, AffForm(-hot(N, N), +float(Tstab - 1)))
                    push!(afs, AffForm(+hot(N, N), -float(Tstab - 1)))
                    b[N] = Tstab
                    A[N, N] = 0.0
                    i2 = 2
                else # Stab
                    push!(afs, AffForm(-hot(N, N), +float(Tstab)))
                    push!(afs, AffForm(+hot(N, N), -float(Tstab)))
                    b[N] = Tstab
                    A[N, N] = 0.0
                    i1 = i2 = 2
                end
            end
        end
        if all(i -> iszero(cases[i]), 1:Nt)
            # display(afs)
            # display([A b])
            # display(abs.(eigvals(A)))
            # display(map(af -> (af.a..., Inf, af.β), afs))
        end
        pf_dom = PolyFunc(afs)
        push!(pieces, Piece(pf_dom, i1, A, b, i2))        
    end
    sys = System(pieces)

    # Initial set
    Nd = Nt ÷ 2
    x_list = Vector{Float64}[]
    for sets in combinations(1:Nt, Nt - Nd)
        x = [fill(float(xinit), Nt); 0.0]
        for i in sets
            x[i] = -xinit
        end
        if Nt == 2*Nd
            push!(x_list, x)
        else
            for i in sets
                y = copy(x)
                y[i] = 0
                push!(x_list, y)
            end
        end
    end
    mlist_init = [x_list, Vector{Float64}[]]

    mpf_inv = MultiPolyFunc([PolyFunc(_AT[]) for loc = 1:M])

    afs_ = [_AT[] for loc = 1:M]
    for i = 1:Nt
        for loc = 1:M
            push!(afs_[loc], AffForm(+hot(i, N), -xsafes[loc]))
            push!(afs_[loc], AffForm(-hot(i, N), -xsafes[loc]))
        end
    end
    mpf_safe = MultiPolyFunc(PolyFunc.(afs_))
    
    return N, M, sys, mlist_init, mpf_inv, mpf_safe
end