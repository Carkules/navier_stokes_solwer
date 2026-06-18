function rk2_step!(v_x_prim, v_y_prim, v_x_start, v_y_start,
                   x_prim, y_prim, nx, ny, dt_prim, Re)

    function rhs_x(v_x, v_y)
        f = zeros(nx + 1, ny) # Poprawiony rozmiar siatki x
        for i in 2:nx, j in 2:(ny-1)
            v_y_dim_x = 0.25 * (v_y[i-1, j] + v_y[i, j] + v_y[i-1, j+1] + v_y[i, j+1])
            du_dx = (v_x[i+1, j] - v_x[i-1, j]) / (2 * x_prim)
            du_dy = (v_x[i, j+1] - v_x[i, j-1]) / (2 * y_prim)
            adwekcja_x = v_x[i, j] * du_dx + v_y_dim_x * du_dy

            d2u_dx2 = (v_x[i+1, j] - 2v_x[i, j] + v_x[i-1, j]) / x_prim^2
            d2u_dy2 = (v_x[i, j+1] - 2v_x[i, j] + v_x[i, j-1]) / y_prim^2

            f[i, j] = -adwekcja_x + (d2u_dx2 + d2u_dy2) / Re
        end
        return f
    end

    function rhs_y(v_x, v_y)
        f = zeros(nx, ny + 1) # Poprawiony rozmiar siatki y
        for i in 2:(nx-1), j in 2:ny
            v_x_dim_y = 0.25 * (v_x[i, j-1] + v_x[i+1, j-1] + v_x[i, j] + v_x[i+1, j])
            dv_dx = (v_y[i+1, j] - v_y[i-1, j]) / (2 * x_prim)
            dv_dy = (v_y[i, j+1] - v_y[i, j-1]) / (2 * y_prim)
            adwekcja_y = v_x_dim_y * dv_dx + v_y[i, j] * dv_dy

            d2v_dx2 = (v_y[i+1, j] - 2v_y[i, j] + v_y[i-1, j]) / x_prim^2
            d2v_dy2 = (v_y[i, j+1] - 2v_y[i, j] + v_y[i, j-1]) / y_prim^2

            f[i, j] = -adwekcja_y + (d2v_dx2 + d2v_dy2) / Re
        end
        return f
    end

    # k1
    k1x = rhs_x(v_x_prim, v_y_prim)
    k1y = rhs_y(v_x_prim, v_y_prim)

    # Predyktor
    vx_temp = v_x_prim .+ dt_prim .* k1x
    vy_temp = v_y_prim .+ dt_prim .* k1y
    nalozenie_bc!(vx_temp, vy_temp) 
    
    # k2
    k2x = rhs_x(vx_temp, vy_temp)
    k2y = rhs_y(vx_temp, vy_temp)

    # Korektor (aktualizacja wnętrza)
    for i in 2:nx, j in 2:(ny-1)
        v_x_start[i, j] = v_x_prim[i, j] + (dt_prim / 2) * (k1x[i, j] + k2x[i, j])
    end

    for i in 2:(nx-1), j in 2:ny
        v_y_start[i, j] = v_y_prim[i, j] + (dt_prim / 2) * (k1y[i, j] + k2y[i, j])
    end
end