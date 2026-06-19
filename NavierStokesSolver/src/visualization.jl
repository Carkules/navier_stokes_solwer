using GLMakie
function animate_flow(v_x_hist, v_y_hist, p_hist, L; fps=20, title="flow.gif")
    
    nx, ny = size(p_hist[1])
    x = range(0, L, length=nx)
    y = range(0, L, length=ny)

    fig = Figure()
    ax = Axis(fig[1, 1])

    X = repeat(x, 1, ny)
    Y = repeat(y', nx, 1)

    scale = 0.05

    function prep(u, v)
        mag = sqrt.(u.^2 .+ v.^2)
        m = maximum(mag)

        m == 0 && return u .* scale, v .* scale

        return (u ./ m) .* scale,
               (v ./ m) .* scale
    end

    function to_center(u, v)
        uc = 0.5 .* (u[1:end-1, :] .+ u[2:end, :])
        vc = 0.5 .* (v[:, 1:end-1] .+ v[:, 2:end])
        return uc, vc
    end

    p0 = p_hist[1]
    u0_raw = v_x_hist[1]
    v0_raw = v_y_hist[1]

    u0, v0 = to_center(u0_raw, v0_raw)
    u0, v0 = prep(u0, v0)

    hm = heatmap!(ax, x, y, p0; colormap=:viridis)

    qv = arrows!(
        ax,
        vec(X), vec(Y),
        vec(u0), vec(v0),
        arrowsize=8,
        linewidth=1,
        color=:white
    )

    limits!(ax, minimum(x), maximum(x), minimum(y), maximum(y))
    step_label = Label(fig[1, 1, Top()], "step = 1")

    record(fig, title, 1:length(p_hist); framerate=fps) do t

        hm[3] = p_hist[t]

        u_raw = v_x_hist[t]
        v_raw = v_y_hist[t]

        u, v = to_center(u_raw, v_raw)
        u, v = prep(u, v)

        qv[3] = vec(u)
        qv[4] = vec(v)

        step_label.text = "step = $t"
    end
end