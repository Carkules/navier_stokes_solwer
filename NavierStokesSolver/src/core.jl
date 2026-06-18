function navier_stokes_solver(U, L, rho, miu, nx, ny, dt, steps, itol = 1e-2; metoda = :euler)

    Re = (rho * U * L) / miu
    ω = 0.7

    x_prim = 1.0 / (nx - 1)
    y_prim = 1.0 / (ny - 1)
    dt_prim = (dt * U) / L

    p_prim = zeros(nx, ny)
    p_new = similar(p_prim)

    v_x_prim = zeros(nx + 1, ny)
    v_y_prim = zeros(nx, ny + 1)

    v_x_start = zeros(nx + 1, ny)
    v_y_start = zeros(nx, ny + 1)

    b = zeros(nx, ny)

    # BC init
    nalozenie_bc_init!(v_x_prim, v_y_prim, v_x_start, v_y_start)

    v_x_hist = Vector{Array{Float64,2}}()
    v_y_hist = Vector{Array{Float64,2}}()
    p_hist = Vector{Array{Float64,2}}()

    push!(v_x_hist, copy(v_x_prim))
    push!(v_y_hist, copy(v_y_prim))
    push!(p_hist, copy(p_prim))

    for step in 1:steps

        v_x_start .= v_x_prim
        v_y_start .= v_y_prim

        if metoda == :euler
            euler_step!(v_x_prim, v_y_prim, v_x_start, v_y_start,
                        x_prim, y_prim, nx, ny, dt_prim, Re)
        elseif metoda == :rk2
            rk2_step!(v_x_prim, v_y_prim, v_x_start, v_y_start,
                      x_prim, y_prim, nx, ny, dt_prim, Re)
        elseif metoda == :rk4
            rk4_step!(v_x_prim, v_y_prim, v_x_start, v_y_start,
                      x_prim, y_prim, nx, ny, dt_prim, Re)
        end

        # Utrwalamy brzegi dla wyznaczonych prędkości wstępnych
        nalozenie_bc!(v_x_start, v_y_start)

        licz_blad!(v_x_start, v_y_start, b, x_prim, y_prim, nx, ny, dt_prim)

        rozwiaz_p!(p_prim, p_new, b, x_prim, y_prim, nx, ny, itol)

        projection!(v_x_prim, v_y_prim, v_x_start, v_y_start,
                    p_prim, x_prim, y_prim, nx, ny, dt_prim)

        # Utrwalamy brzegi dla prędkości ostatecznych
        nalozenie_bc!(v_x_prim, v_y_prim)

        sprawdz_dywerencje(v_x_prim, v_y_prim, x_prim, y_prim, nx, ny, itol, step)

        push!(v_x_hist, copy(v_x_prim))
        push!(v_y_hist, copy(v_y_prim))
        push!(p_hist, copy(p_prim))
    end

    return v_x_hist, v_y_hist, p_hist
end