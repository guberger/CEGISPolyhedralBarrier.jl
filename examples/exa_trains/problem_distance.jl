# hot(i, N)'*x[1:N] = x[i]
hot(i, N) = [j == i ? 1.0 : 0.0 for j = 1:N]

function build_problem(
        Nt, α, β, lim_lo, lim_up, vinit, xsafe_lo, xsafe_up, vsafe
    )
    N = 2*Nt
    M = 1

    # System
    pieces = Piece{_PT,Matrix{Float64},Vector{Float64}}[]
    for cases in Iterators.product([(-1, 0, 1) for i = 1:Nt]...)
        afs = [
            AffForm([+ones(Nt); zeros(Nt)], 0.0),
            AffForm([-ones(Nt); zeros(Nt)], 0.0),
            AffForm([zeros(Nt); +ones(Nt)], 0.0),
            AffForm([zeros(Nt); -ones(Nt)], 0.0)
        ]
        # afs = _AT[]
        A = Matrix{Float64}(I, N, N)
        for i = 1:Nt
            A[i, Nt + i] = 1.0
        end
        b = zeros(N)
        for (i, case) in enumerate(cases) # dxcurr = xnext - xcurr
            icurr = mod(i - 1, Nt) + 1
            iprev = mod(i - 2, Nt) + 1
            h = α*hot(icurr, N) + β*hot(Nt + icurr, N)
            if case == -1 # α*dxcurr + β*dvcurr ≥ lim_up --> vcurr += lim_up
                push!(afs, AffForm(-h, +lim_up))
                b[icurr] -= lim_up/2
                b[Nt + icurr] -= lim_up
                b[iprev] += lim_up/2
                b[Nt + iprev] += lim_up
            elseif case == 1 # α*dxcurr + β*dvcurr ≤ lim_lo --> vcurr += lim_lo
                push!(afs, AffForm(+h, -lim_lo))
                b[icurr] -= lim_lo/2
                b[Nt + icurr] -= lim_lo
                b[iprev] += lim_lo/2
                b[Nt + iprev] += lim_lo
            else # lim_lo ≤ α*dxcurr + β*dvcurr ≤ lim_up --> vcurr += α*dxcurr + β*dvcurr
                push!(afs, AffForm(-h, +lim_lo))
                push!(afs, AffForm(+h, -lim_up))
                A[icurr, icurr] -= α/2
                A[icurr, Nt + icurr] -= β/2
                A[Nt + icurr, icurr] -= α
                A[Nt + icurr, Nt + icurr] -= β
                A[iprev, icurr] += α/2
                A[iprev, Nt + icurr] += β/2
                A[Nt + iprev, icurr] += α
                A[Nt + iprev, Nt + icurr] += β
            end
        end
        if all(iszero, cases)
            # display(afs)
            display([A b])
            display(abs.(eigvals(A)))
            display(map(af -> (af.a..., NaN, af.β), afs))
        end
        pf_dom = PolyFunc(afs)
        push!(pieces, Piece(pf_dom, 1, A, b, 1))
    end
    sys = System(pieces)

    # Initial set
    Nd = Nt ÷ 2
    x_list = Vector{Float64}[]
    for sets in combinations(1:Nt, Nt - Nd)
        x = [zeros(Nt); fill(float(vinit), Nt)]
        for i in sets
            x[Nt + i] = -vinit
        end
        if Nt == 2*Nd
            push!(x_list, x)
        else
            for i in sets
                y = copy(x)
                y[Nt + i] = 0
                push!(x_list, y)
            end
        end
    end
    mlist_init = [x_list]

    mpf_inv = MultiPolyFunc([PolyFunc(_AT[]) for loc = 1:M])

    afs_ = [_AT[] for loc = 1:M]
    for i = 1:Nt
        for loc = 1:M
            push!(afs_[loc], AffForm(+hot(i, N), -xsafe_up))
            push!(afs_[loc], AffForm(-hot(i, N), +xsafe_lo))
            push!(afs_[loc], AffForm(+hot(Nt + i, N), -vsafe))
            push!(afs_[loc], AffForm(-hot(Nt + i, N), -vsafe))
        end
    end
    mpf_safe = MultiPolyFunc(PolyFunc.(afs_))
    
    return N, M, sys, mlist_init, mpf_inv, mpf_safe
end