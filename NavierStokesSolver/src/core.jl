"""
    navier_stokes_solver(U, L, rho, miu, nx, ny, dt, steps, itol = 1e-2; metoda = :euler)

Główna funkcja solvera rozwiązująca dwuwymiarowe, nieściśliwe równania Naviera-Stokesa w formie bezwymiarowej na siatce przesuniętej (staggered grid).

# Argumenty pozycyjne
- `U::Real`: Prędkość wpływu do kanału.
- `L::Real`: Długość boku kanału.
- `rho::Real`: Gęstość płynu (rho).
- `miu::Real`: Lepkość dynamiczna płynu (miu).
- `nx::Int`, `ny::Int`: Liczba węzłów siatki odpowiednio w kierunku X i Y.
- `dt::Real`: FKrok czasowy symulacji.
- `steps::Int`: Całkowita liczba kroków czasowych do wykonania.
- `itol::Real`: Tolerancja zbieżności dla solvera ciśnienia oraz testu bezdywergencyjności pola prędkości (domyślnie: `1e-2`).

# Argumenty nazwane
- `metoda::Symbol`: Schemat dyskretyzacji członów adwekcyjno-dyfuzyjnych w czasie. Dostępne opcje:
  - `:euler` - jawny schemat Eulera 1. rzędu (domyślny),
  - `:rk2` - schemat Rungego-Kutty 2. rzędu,
  - `:rk4` - schemat Rungego-Kutty 4. rzędu.

# Zwraca
Trzyelementową krotkę (tuple) zawierającą historię stanów pól w każdym kroku czasowym:
1. `v_x_hist::Vector{Matrix{Float64}}`: Historia rozkładu prędkości poziomej zlokalizowanej na pionowych krawędziach komórek.
2. `v_y_hist::Vector{Matrix{Float64}}`: Historia rozkładu prędkości pionowej zlokalizowanej na poziomych krawędziach komórek.
3. `p_hist::Vector{Matrix{Float64}}`: Historia rozkładu pola ciśnienia zlokalizowanego w centrach komórek.

# Szczegóły algorytmu
- Metoda rzutowania: W każdym kroku czasowym obliczana jest prędkość pomocnicza (wstępna) bez uwzględnienia gradientu ciśnienia, po czym rozwiązywane jest równanie Poissona dla ciśnienia w celu skorygowania pola prędkości do postaci spełniającej warunek ciągłości (zerowej dywergencji).
- Siatka typu Staggered: Zapobiega powstawaniu niefizycznych oscylacji ciśnienia poprzez obliczanie zmiennych skalarnych w środkach komórek, a wektorowych na ich krawędziach.

# Przykład użycia
```julia
# Symulacja dla liczby Reynoldsa Re = 100 przy użyciu metody RK4
v_x_hist, v_y_hist, p_hist = navier_stokes_solver(1.0, 1.0, 1.0, 0.01, 20, 20, 0.06, 500; metoda=:rk4)
"""

function navier_stokes_solver(U, L, rho, miu, nx, ny, dt, steps, itol = 1e-2; metoda = :euler)

    Re = (rho * U * L) / miu

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

        nalozenie_bc!(v_x_start, v_y_start)

        licz_blad!(v_x_start, v_y_start, b, x_prim, y_prim, nx, ny, dt_prim)

        rozwiaz_p!(p_prim, p_new, b, x_prim, y_prim, nx, ny, itol)

        projection!(v_x_prim, v_y_prim, v_x_start, v_y_start,
                    p_prim, x_prim, y_prim, nx, ny, dt_prim)

        nalozenie_bc!(v_x_prim, v_y_prim)

        sprawdz_dywerencje(v_x_prim, v_y_prim, x_prim, y_prim, nx, ny, itol, step)

        push!(v_x_hist, copy(v_x_prim))
        push!(v_y_hist, copy(v_y_prim))
        push!(p_hist, copy(p_prim))
    end

    return v_x_hist, v_y_hist, p_hist
end