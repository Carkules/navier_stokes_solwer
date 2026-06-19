# navier_stokes_solwer
NavierStokesSolver to pakiet funkcji napisany w języku Julia opisujący zagadnienie numerycznej mechaniki płynów.

# NavierStokesSolver.jl

Solver równań Naviera-Stokesa dla nieściśliwego przepływu 2D w kwadratowym kanale, zaimplementowany w języku Julia.

## Metoda numeryczna

Solver wykorzystuje metodę projekcyjną na siatce staggered MAC:

1. Predyktor — całkowanie równania pędu bez członu ciśnienia (adwekcja + dyfuzja)
2. Równanie Poissona — wyznaczenie ciśnienia z warunku nieściśliwości ($\nabla \cdot v = 0$)
3. Korekcja — poprawka prędkości o gradient ciśnienia

Równanie rozwiązywane w formie bezwymiarowej:

$$\frac{\partial \vec{v}}{\partial t} + (\vec{v} \cdot \nabla)\vec{v} = -\nabla p + \frac{1}{Re} \Delta \vec{v}$$

## Struktura pakietu

- NavierStokesSolver
    - src
        - NavierStokesSolver.jl — moduł główny
        - core.jl — główna pętla czasowa
        - boundary.jl — warunki brzegowe
        - euler.jl — schemat Eulera (1. rzędu)
        - rk2.jl — schemat Runge-Kutta 2. rzędu
        - rk4.jl — schemat Runge-Kutta 4. rzędu
        - pressure.jl — solver Poissona i korekcja prędkości
        - diagnostic.jl — sprawdzanie dywergencji
        - visualization.jl — animacja wyników (GLMakie)
    - Project.toml
    - Manifest.toml

## Instalacja
Żeby zainstalować pakiet NavierStokesSolver należy wykonać następujące kroki:
```julia
import Pkg
Pkg.add(url="https://github.com/Carkules/navier_stokes_solwer", subdir="NavierStokesSolver")
```

## Użycie

```julia
using NavierStokesSolver

v_x_hist, v_y_hist, p_hist = navier_stokes_solver(
    U,      # prędkość referencyjna [m/s]
    L,      # długość boku kanału [m]
    rho,    # gęstość płynu [kg/m³]
    miu,    # lepkość dynamiczna [Pa·s]
    nx,     # liczba węzłów siatki w kierunku X
    ny,     # liczba węzłów siatki w kierunku Y
    dt,     # krok czasowy [s]
    steps,  # liczba kroków czasowych
    itol;   # tolerancja zbieżności (domyślnie 1e-2)
    metoda = :euler  # schemat czasowy: :euler, :rk2, :rk4
)
```

### Przykład

```julia
v_x_hist, v_y_hist, p_hist = navier_stokes_solver(
    1.0, 1.0, 1.0, 0.01,
    20, 20,
    0.001, 200;
    metoda = :rk2
)
```

## Warunki brzegowe

Wlot (lewa ściana): płyn wpływa z prędkością bezwymiarową `v_x = 1`, `v_y = 0`, ciśnienie ustala się samo.

Wylot (prawa ściana): swobodny wypływ — prędkości z wnętrza, ciśnienie referencyjne `p = 0`.

**Górna i dolna ściana:** warunek no-slip — `v_x = 0`, `v_y = 0`, ciśnienie z środka.

## Zwracane wartości

Funkcja `navier_stokes_solver` zwraca historię symulacji jako trzy wektory macierzy:

- `v_x_hist` — historia prędkości poziomej, rozmiar każdej macierzy `(nx+1, ny)`
- `v_y_hist` — historia prędkości pionowej, rozmiar `(nx, ny+1)`
- `p_hist` — historia ciśnienia, rozmiar `(nx, ny)`

## Animacja
Do tworzenia animacji z otrzymanych wyników służy funkckcja `animate_flow`. Tworzy ona plik o rozszerzeniu .gif zawierający animację przepływu.

```julia
animate_flow(
    v_x_hist, # Historia rozkładu prędkości poziomej w kolejnych krokach czasowych.
    v_y_hist, # Historia rozkładu prędkości pionowej w kolejnych krokach czasowych.
    p_hist,   # Historia rozkładu pola ciśnienia w kolejnych krokach czasowych.
    L;        # Długość boku kanału.
    fps=20,   # Liczba klatek na sekundę generowanej animacji GIF.
    title="flow.gif")   # Nazwa lub ścieżka pliku wynikowego wraz z rozszerzeniem `.gif`.
```
Przykład użycia

```julia
animate_flow(v_x_hist, v_y_hist, p_hist, 1.0;
    fps=20, title="flow.gif")
```

## Wymagania
- Julia 1.12+
- GLMakie 0.13.11


















