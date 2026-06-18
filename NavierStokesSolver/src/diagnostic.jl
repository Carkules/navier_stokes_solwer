function sprawdz_dywerencje(v_x_prim, v_y_prim, x_prim, y_prim, nx, ny, itol, step)
    max_blad = 0.0

    for i in 2:(nx-1), j in 2:(ny-1)
        blad =
            (v_x_prim[i+1, j] - v_x_prim[i, j]) / x_prim +
            (v_y_prim[i, j+1] - v_y_prim[i, j]) / y_prim

        max_blad = max(max_blad, abs(blad))
    end

    if max_blad > itol
        println("Krok $step: dywergencja = $max_blad > itol = $itol")
    end

    return max_blad
end