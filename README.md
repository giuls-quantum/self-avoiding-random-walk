# Self-Avoiding Random Walk Simulation

This project simulates **self-avoiding random walks (SARW)** on a 2D square lattice. It analyzes walk lengths, dead-end occurrences, and distance from the origin over multiple simulations.

## Features

- Simulates **self-avoiding random walks** on a square grid of size `L x L`.
- Implements **periodic and standard boundary conditions**.
- Saves simulation data in `.dat` files:
  - `positions.dat` – positions of the first walk.
  - `averages.dat` – average distance from origin over multiple walks.
  - `distance.dat` – squared distance vs. time for short walks.
  - `probability.dat` – probability distribution of walk lengths.
- Generates **plots in PDF** using gnuplot:
  - `walk.pdf`, `averageR.pdf`, `distance.pdf`, `prob.pdf`.
- Counts **distinct walks** and **rotation-invariant walks**.

## Requirements

- C compiler (e.g., `gcc`)  
- `gnuplot` installed for plotting  

Optional: Python with `matplotlib` if you want to create plots without gnuplot.

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/self_avoiding_random_walk.git
cd self_avoiding_random_walk
```

Compile the program:

```bash
gcc -o selfa_random_walk src/selfa_random_walk.c -lm
```

## Usage

Run the simulation:

```bash
./selfa_random_walk
```

- The program will generate .dat files in the data/ folder.
- PDF plots will be created in the plots/ folder (if gnuplot is installed).

## Output

- positions.dat – trajectory of the first simulated walk.
- averages.dat – average squared distance vs. attempt number.
- distance.dat – short-time squared distance from origin.
- probability.dat – probability of walks ending at a certain length.
- PDF plots – visualizations of the data for quick inspection.

## Repository Structure

```bash
self_avoiding_random_walk/
│
├── src/           # C source code
│   └─ selfa_random_walk.c
├── data/          # Simulation output (.dat)
├── plots/         # Generated plots (.pdf)
├── README.md
└── LICENSE
```

## License

This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.

