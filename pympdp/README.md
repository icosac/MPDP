# PyMPDP -- Python Package from MPDP


# Running the examples without installing the package

To run the examples without installing the package, navigate to the `pympdp/examples` folder and run the desired example using Python. For example, to run the MPMD example, navigate to the `pympdp/examples/MPMD` folder and run:
```bash
python mpmd.py
```


## Installation 

From the main folder of the repository, run the following commands.

First create a virtual environment:
```bash
python3 -m venv venv
source venv/bin/activate
```
Then install the required dependencies:
```bash
python3 -m pip install -r pympdp/requirements.txt
```

Then install the package using pip:
```bash
python3 -m pip install .
```

To check that the installation was successful open a Python shell and run:
```python
import pympdp
```
It should not return any error.

### Optional dependencies

You may need to install Tinker for matplotlib to work properly. On Ubuntu, you can install it using:
```bash
sudo apt-get install python3-tk
```
On MacOS, you can install it using Homebrew:
```bash
brew install python-tk
```


## Running the examples

Once the package is installed, you can run the examples in the `pympdp/examples` folder. For example, to run the MPMD example, navigate to the `pympdp/examples/MPMD` folder and run:
```bash
python mpmd.py
```

## Usage 

The main class of the package is `DP`, which can be imported from the `pympdp.dp` module. Here is a simple example of how to use it:

```python
from pympdp.dp import DP


if __name__ == "__main__":
    # Define the points and which angles are allowed at each point
    points = [(0, 0), (1, 1), (2, 2)]
    fixed_angles = [True, False, True]
    # If some angles are fixed, define them here
    angles = [0.0, 0.0, 1.57]

    # Create a DP instance specifying the points, the angles, the number of discretizations and refinements, and the maximum curvature
    dp_instance = DP(points = points,
            fixed_angles = fixed_angles,
            def_thetas = angles,
            discretizations = 4,
            refinements = 1,
            k_max = 1.0
    )

    # Solve the problem
    optimal_angles = dp_instance.solve()

```

A similar example can be found in `pympdp/dp/dp.py` by simply running the file:
```bash
python pympdp/dp/dp.py
```

Once the problem is solved, you can also visualize the dynamic programming matrix by running:
```python
dp_instance.visualize_dp_matrix(show_optimal_path=True)
```
Which will open a web page with an interactive visualization.

If you run the faster C++ solver you can still reuse this visualization. Dump the matrix by
calling `DP::exportVisualizationData("dp_snapshot.json")` in C++, then run:
```bash
python3 -m pympdp.dp.visualize_json dp_snapshot.json
```
The command reads the JSON file and produces the same HTML dashboard without rerunning the
Python solver.

Alternative, you can run the following command for a static visualization in the terminal:
```python
dp_instance.print_dp_matrix()
```

