# PotentialSphere

**PotentialSphere** is an interactive Processing sketch that simulates the equilibrium configuration of `n` points on a sphere under repulsive forces. It visually demonstrates how particles distribute themselves to minimize pairwise potential energy while constrained to a spherical surface.

This is a local optimization process governed by gradient descent-like dynamics, with constraint projection to the sphere at each step.

## Features

- Evenly distributes points on a sphere using physical simulation
- Uses inverse-square repulsive forces to drive the system to a low-energy configuration
- Projects motion back onto the sphere surface after each update
- Dynamically draws lines from each node to its `k` nearest neighbors (`k = ceil(log₂(n))`)
- 3D interactive visualization with rotation, zoom, and node dragging
- Optional density grid drawn as semi-transparent blue polygons on the sphere
- Real-time display of system energy

## Mathematical Background

This simulation relates to the **Thomson Problem** — a classic question in mathematical physics that asks:

> *What is the minimum energy configuration of `n` mutually repelling particles constrained to the surface of a sphere?*

Each pair of particles repels according to an inverse-square law, and the system seeks to minimize the total potential energy. While exact solutions are known for small `n` (e.g., 2 to 5), finding the **global minimum configuration for arbitrary `n` is an open problem**. The number of local minima increases rapidly with `n`, and the energy landscape becomes highly non-convex.

**PotentialSphere** provides a numerical method to approximate such configurations by simulating gradient flow, constrained to the sphere. It finds **locally minimal energy arrangements**, which may or may not be globally optimal.

## Controls

- **Left Mouse Drag** — Rotate the view
- **Right Mouse Drag** — Drag a specific node (automatically selects closest node on press)
- **Scroll Wheel** — Zoom in and out

## Parameters

- `nNodes` — Number of particles distributed on the sphere
- `nNeighbors` — Automatically set to `ceil(log₂(nNodes))`, controls nearest-neighbor lines
- `sphereRadius` — Radius of the constraint sphere
- `forceStrength` — Strength of repulsion (more negative = stronger repulsion)
- `scaleFactor` — Affects zoom level for rendering

## Code Structure

- `Node` class — Represents particles with position, velocity, and color
- `TransformationMatrix` class — Encapsulates 3D view rotation logic
- `calculateForces()` — Computes pairwise inverse-square repulsive forces
- `calculateTotalEnergy()` — Computes system potential energy
- `drawDensityPlot()` — Draws a static mesh on the sphere (no scalar field computation)
- `findClosestNodes()` — Finds the `k` nearest neighbors for graph drawing
- Event functions handle interactive view manipulation and dragging

## Requirements

- [Processing](https://processing.org/) (Java mode)

## Notes

- The simulation finds a **local minimum** of the system's potential energy. It does not guarantee a globally optimal or symmetrical configuration.
- You can tweak the `nNodes` parameter to explore how equilibrium changes with system size.
- All nodes remain constrained to the sphere throughout the simulation via normalization.

## License

Open source and free to use for educational, scientific, and visualization purposes.

---

Explore the physics of curved surfaces and the beauty of emergent structure with **PotentialSphere**.
