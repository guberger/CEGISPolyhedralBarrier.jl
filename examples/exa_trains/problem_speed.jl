using Combinatorics

# hot(i, N)'*x[1:N] = x[i]
hot(i, N) = [j == i ? 1.0 : 0.0 for j = 1:N]

function build_problem(Nt, α, Tstab, lim_lo, lim_up, vinit, vsafes)
    N = Nt + 1

    # System
    pieces = Piece[]
    for cases in Iterators.product([(-1, 0, 1) for i = 1:N]...)
        afs_dom = [
            AffForm([+ones(Nt); 0.0], 0.0),
            AffForm([-ones(Nt); 0.0], 0.0)
        ]
        A = Matrix{Float64}(I, N, N)
        b = zeros(N)
        i1 = i2 = 1
        for (i, case) in enumerate(cases) # dvcurr = vnext - vcurr
            if i < N
                icurr = mod(i - 1, Nt) + 1
                iprev = mod(i - 2, Nt) + 1
                if case == -1 # dvcurr ≥ lim_up/α --> vcurr += lim_up
                    push!(afs_dom, AffForm(-hot(icurr, N), +lim_up/α))
                    b[icurr] -= lim_up
                    b[iprev] += lim_up
                elseif case == 1 # dvcurr ≤ lim_lo/α --> vcurr += lim_lo
                    push!(afs_dom, AffForm(+hot(icurr, N), -lim_lo/α))
                    b[icurr] -= lim_lo
                    b[iprev] += lim_lo
                else # lim_lo/α ≤ dvcurr ≤ lim_up/α --> vcurr += dvcurr*α
                    push!(afs_dom, AffForm(-hot(icurr, N), +lim_lo/α))
                    push!(afs_dom, AffForm(+hot(icurr, N), -lim_up/α))
                    A[icurr, icurr] -= α
                    A[iprev, icurr] += α
                end
            else
                if case == -1 # Trans & t ≤ Tstab - 1
                    push!(afs_dom, AffForm(+hot(N, N), -float(Tstab - 2)))
                    push!(afs_dom, AffForm(-hot(N, N), 0.0))
                    b[N] = 1.0
                elseif case == 1 # Trans & t ≥ Tstab - 1
                    push!(afs_dom, AffForm(-hot(N, N), +float(Tstab - 1)))
                    push!(afs_dom, AffForm(+hot(N, N), -float(Tstab - 1)))
                    b[N] = Tstab
                    A[N, N] = 0.0
                    i2 = 2
                else # Stab
                    push!(afs_dom, AffForm(-hot(N, N), +float(Tstab)))
                    push!(afs_dom, AffForm(+hot(N, N), -float(Tstab)))
                    b[N] = Tstab
                    A[N, N] = 0.0
                    i1 = i2 = 2
                end
            end
        end
        if all(i -> iszero(cases[i]), 1:Nt)
            # display([A b])
            # display(abs.(eigvals(A)))
            # display(map(af -> (af.a..., NaN, af.β), afs_dom))
        end
        push!(pieces, Piece(afs_dom, i1, A, b, i2))
    end

    # Initial set
    Nd = Nt ÷ 2
    states_init = State[]
    for sets in combinations(1:Nt, Nt - Nd)
        x = [fill(float(vinit), Nt); 0.0]
        for i in sets
            x[i] = -vinit
        end
        if Nt == 2*Nd
            push!(states_init, State(1, x))
        else
            for i in sets
                y = copy(x)
                y[i] = 0
                push!(states_init, State(1, y))
            end
        end
    end

    gfs_safe = [
        [
            GenForm(loc, AffForm(+hot(i, N), -vsafes[loc]))
            for i = 1:Nt, loc = 1:2
        ]...,
        [
            GenForm(loc, AffForm(-hot(i, N), -vsafes[loc]))
            for i = 1:Nt, loc = 1:2
        ]...
    ]
    
    return BarrierProblem(
        N, pieces,
        GenForm[], # gfs_inv
        gfs_safe,
        states_init,
        abs(vsafes[2])*α/5, # ϵ
        1e-8 # δ
    )
end