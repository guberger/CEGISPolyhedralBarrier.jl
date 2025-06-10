using Combinatorics

# hot(i, N)'*x[1:N] = x[i]
hot(i, N) = [j == i ? 1.0 : 0.0 for j = 1:N]

function build_problem(Nt, α, β, lim_lo, lim_up, vinit,
                       xsafe_lo, xsafe_up, vsafe)
    N = 2*Nt

    # System
    pieces = Piece[]
    for cases in Iterators.product([(-1, 0, 1) for i = 1:Nt]...)
        afs_dom = [
            AffForm([+ones(Nt); zeros(Nt)], 0.0),
            AffForm([-ones(Nt); zeros(Nt)], 0.0),
            AffForm([zeros(Nt); +ones(Nt)], 0.0),
            AffForm([zeros(Nt); -ones(Nt)], 0.0)
        ]
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
                push!(afs_dom, AffForm(-h, +lim_up))
                b[icurr] -= lim_up/2
                b[Nt + icurr] -= lim_up
                b[iprev] += lim_up/2
                b[Nt + iprev] += lim_up
            elseif case == 1 # α*dxcurr + β*dvcurr ≤ lim_lo --> vcurr += lim_lo
                push!(afs_dom, AffForm(+h, -lim_lo))
                b[icurr] -= lim_lo/2
                b[Nt + icurr] -= lim_lo
                b[iprev] += lim_lo/2
                b[Nt + iprev] += lim_lo
            else # lim_lo ≤ α*dxcurr + β*dvcurr ≤ lim_up --> vcurr += α*dxcurr + β*dvcurr
                push!(afs_dom, AffForm(-h, +lim_lo))
                push!(afs_dom, AffForm(+h, -lim_up))
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
            # display([A b])
            # display(abs.(eigvals(A)))
            # display(map(af -> (af.a..., NaN, af.β), afs))
        end
        push!(pieces, Piece(afs_dom, 1, A, b, 1))
    end

    # Initial set
    Nd = Nt ÷ 2
    states_init = State[]
    for sets in combinations(1:Nt, Nt - Nd)
        x = [zeros(Nt); fill(float(vinit), Nt)]
        for i in sets
            x[Nt + i] = -vinit
        end
        if Nt == 2*Nd
            push!(states_init, State(1, x))
        else
            for i in sets
                y = copy(x)
                y[Nt + i] = 0
                push!(states_init, State(1, y))
            end
        end
    end

    gfs_safe = [
        [
            GenForm(1, AffForm(+hot(i, N), -xsafe_up)) for i = 1:Nt
        ]...,
        [
            GenForm(1, AffForm(-hot(i, N), +xsafe_lo)) for i = 1:Nt
        ]...,
        [
            GenForm(1, AffForm(+hot(Nt + i, N), -vsafe)) for i = 1:Nt
        ]...,
        [
            GenForm(1, AffForm(-hot(Nt + i, N), -vsafe)) for i = 1:Nt
        ]...
    ]
    
    return BarrierProblem(
        N, pieces,
        GenForm[], # gfs_inv
        gfs_safe,
        states_init,
        vsafe*α/5, # ϵ
        1e-8 # δ
    )
end