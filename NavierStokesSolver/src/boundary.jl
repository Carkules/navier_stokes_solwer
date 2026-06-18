function nalozenie_bc_init!(v_x_prim, v_y_prim, v_x_start, v_y_start)
    nalozenie_bc!(v_x_prim, v_y_prim)
    nalozenie_bc!(v_x_start, v_y_start)
end

function nalozenie_bc!(v_x, v_y)
    v_x[1, :] .= 1.0
    v_y[1, :] .= 0.0
    v_x[end, :] .= v_x[end-1, :]
    v_y[end, :] .= v_y[end-1, :]

    v_x[:, 1] .= 0.0
    v_x[:, end] .= 0.0
    v_y[:, 1] .= 0.0
    v_y[:, end] .= 0.0
end