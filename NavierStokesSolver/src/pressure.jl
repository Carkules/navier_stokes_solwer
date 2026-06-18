function licz_blad!(v_x_start, v_y_start, b, x_prim, y_prim, nx, ny, dt_prim)
    for i in 2:(nx-1), j in 2:(ny-1)
        blad = (v_x_start[i+1, j] - v_x_start[i, j]) / x_prim +
               (v_y_start[i, j+1] - v_y_start[i, j]) / y_prim

        b[i, j] = blad / dt_prim
    end
end

function rozwiaz_p!(p_prim, p_new, b, x_prim, y_prim, nx, ny, itol)
    for iteracja in 1:500
        p_new .= p_prim

        for i in 2:(nx-1), j in 2:(ny-1)
            p_new[i, j] =
                ((p_prim[i+1, j] + p_prim[i-1, j]) / x_prim^2 +
                 (p_prim[i, j+1] + p_prim[i, j-1]) / y_prim^2 -
                 b[i, j]) / (2 / x_prim^2 + 2 / y_prim^2)
        end

        p_new[1, :] .= p_new[2, :]
        p_new[end, :] .= 0.0
        p_new[:, 1] .= p_new[:, 2]
        p_new[:, end] .= p_new[:, end-1]

        roznica = maximum(abs.(p_new .- p_prim))
        p_prim .= p_new

        if roznica < itol
            break
        end
    end
end

function projection!(v_x_prim, v_y_prim, v_x_start, v_y_start,
                     p_prim, x_prim, y_prim, nx, ny, dt_prim)

    for i in 2:nx, j in 1:ny
        v_x_prim[i, j] =
            v_x_start[i, j] - dt_prim * (p_prim[i, j] - p_prim[i-1, j]) / x_prim
    end

    for i in 1:nx, j in 2:ny # Poprawiono z 1:ny na 1:nx
        v_y_prim[i, j] =
            v_y_start[i, j] - dt_prim * (p_prim[i, j] - p_prim[i, j-1]) / y_prim
    end
end