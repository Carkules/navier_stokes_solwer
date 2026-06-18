function euler_step!(v_x_prim, v_y_prim, v_x_start, v_y_start,
                     x_prim, y_prim, nx, ny, dt_prim, Re)

    for i in 2:nx, j in 2:(ny-1)

        v_y_dim_x = 0.25 * (v_y_prim[i-1, j] + v_y_prim[i, j] +
                            v_y_prim[i-1, j+1] + v_y_prim[i, j+1])

        du_dx = (v_x_prim[i+1, j] - v_x_prim[i-1, j]) / (2x_prim)
        du_dy = (v_x_prim[i, j+1] - v_x_prim[i, j-1]) / (2y_prim)

        adwekcja_x = v_x_prim[i, j] * du_dx + v_y_dim_x * du_dy

        d2u_dx2 = (v_x_prim[i+1, j] - 2v_x_prim[i, j] + v_x_prim[i-1, j]) / x_prim^2
        d2u_dy2 = (v_x_prim[i, j+1] - 2v_x_prim[i, j] + v_x_prim[i, j-1]) / y_prim^2

        v_x_start[i, j] = v_x_prim[i, j] +
            dt_prim * (-adwekcja_x + (d2u_dx2 + d2u_dy2)/Re)
    end

    for i in 2:(nx-1), j in 2:ny

        v_x_dim_y = 0.25 * (v_x_prim[i, j-1] + v_x_prim[i+1, j-1] +
                            v_x_prim[i, j] + v_x_prim[i+1, j])

        dv_dx = (v_y_prim[i+1, j] - v_y_prim[i-1, j]) / (2x_prim)
        dv_dy = (v_y_prim[i, j+1] - v_y_prim[i, j-1]) / (2y_prim)

        adwekcja_y = v_x_dim_y * dv_dx + v_y_prim[i, j] * dv_dy

        d2v_dx2 = (v_y_prim[i+1, j] - 2v_y_prim[i, j] + v_y_prim[i-1, j]) / x_prim^2
        d2v_dy2 = (v_y_prim[i, j+1] - 2v_y_prim[i, j] + v_y_prim[i, j-1]) / y_prim^2

        v_y_start[i, j] = v_y_prim[i, j] +
            dt_prim * (-adwekcja_y + (d2v_dx2 + d2v_dy2)/Re)
    end
end