# Computational Physics Exercises

Welcome to my Computational Physics Exercises repository! This repository contains all the exercises from the course 'Computational Physics' at Universidad de Alicante. The aim is to provide a comprehensive set of problems and solutions to showcase the computational problems solved in this course. I will translate the exercises into English, and add solutions to the exercises. In the future, I also plan to have a Spanish version to help students of the course, which name has changed to 'Introduction to modelization in Physics', but which contents remain the same.

## Table of Contents
- [Introduction](#introduction)
- [Course Outline](#course-outline)
- [Setup and Installation](#setup-and-installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Introduction
This repository contains exercises and solutions from the Computational Physics course. The exercises cover various topics, including numerical methods and simulations, with a focus on their application to physical problems.

## Course Outline

### Block I - Fundamental Numerical Methods

1. **Tridiagonal Matrix Solver**: Implement a function to solve tridiagonal systems of equations and compare the performance of different solving methods.

2. **Laplace Equation Solver**: Solve the Laplace equation with Dirichlet and Neumann boundary conditions using sparse matrix methods and iterative solvers.

3. **Damped Harmonic Oscillator**: Solve the differential equation for a damped harmonic oscillator using explicit, implicit, and semi-implicit methods. Compare numerical results with the exact solution.

4. **Poisson's Equation in 1D and 2D**: Implement numerical methods to solve Poisson's equation in one and two dimensions.

5. **Damped Waves and Telegraph Equation**: Solve wave equations with damping and source terms using explicit and implicit methods.

### Block II - Advanced Topics and Applications

1. **Crystalline Lattice Generation**: Create a program to generate different types of crystalline lattices and visualize them in 3D.

2. **Orbital Mechanics Simulation**: Calculate the Earth's orbit around the Sun using numerical methods to solve the equations of motion.

3. **Lennard-Jones Potential**: Calculate the total energy of a system using the Lennard-Jones potential and simulate molecular dynamics.

4. **Force Calculation with Pair Potentials**: Extend the previous exercise to calculate forces on particles and simulate their trajectories using the Verlet algorithm.

5. **Monte Carlo Integration**: Implement the Monte Carlo method to integrate functions and explore statistical errors.

6. **Metropolis Monte Carlo Simulation**: Implement the Metropolis Monte Carlo method to simulate an ideal gas and study its energy distribution at different temperatures.

## Setup and Installation
To get started with the exercises, follow these steps:

1. **Clone the repository:**
    ```bash
    git clone https://github.com/ptr-25/computational_physics.git
    cd computational_physics
    ```

2. **Install the required dependencies:**
    - Ensure you have Python installed. You can download it from [python.org](https://www.python.org/).
    - Install the necessary Python packages:
    ```bash
    pip install -r requirements.txt
    ```

## Usage
Each exercise is contained in its own directory with the following structure:
|-- exercise-1
| |-- README.md
| |-- solution.py
| |-- data
| |-- input_data.txt
| |-- output_data.txt
|-- exercise-2
| |-- README.md
| |-- solution.py
| |-- data
| |-- input_data.txt
| |-- output_data.txt

- Navigate to the directory of the exercise you wish to work on.
- Read the `README.md` file in the exercise directory for specific instructions.
- Run the provided solutions or implement your own.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
For any questions or suggestions, feel free to open an issue or contact me directly at [rtwerdy@gmail.com](mailto:rtwerdy@gmail.com).

---

Happy computing!


