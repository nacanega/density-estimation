# density-estimation
Atmospheric Density Estimation from Large Satellite Constellation Navigation Data

---

## Essential Features Needed Before Release
Before creating the initial release, a few features are needed.
* Implement SNC and DMC in LKF_RTSpre
* Implement passing of parameters into filter for state equations.
* Implement drag into state equations.
* Implement density estimation into state equations.

## Future Quality of Life Improvements
* Implement robust input validation.
* Implement classes for each object.
* Implement SRIF
* Implement Information-based smoother

## Desired Features
* Implement Wind
* Implement More Sophisticated Truth Generation

## Immediate Goals

TODO Finish basic LKF_RTSpre
- [x] Variable Outputs
- [x] Argument Validation
- [x] Function Validation
- [x] Implement Structures
- [x] Endstate Output
- [x] Divergence Check
- [x] Tolerance Check
- [x] Iteration Limit
- [x] Nonconstant Q
- [x] Nonconstant R
- [x] Nonconstant H
- [ ] SNC
- [ ] DMC
- [ ] De-Orbit Check
- [ ] Alternate version with taylor series
- [ ] Alternate version with looped interval integration

TODO Make a call script
- [x] Evaluate Each satellite individually
- [ ] Evaluate Each satellite combined N satellite with different model and drag

TODO
- [ ] Block integration function with sparse matrices 

TODO Make a basic visualizer
- [x] Show Plots of components
- [x] Save results to process later
- [] Make a visualizer app with 3d animated plots

TODO Process Noise
- [x] DPN - Diagonal Process Noise
- [ ] SNC - State Noise Compensation

TODO Implement Drag
- [ ] Basic Drag Model
- [ ] Fix Plane of Earth Rotation
- [ ] Complex Drag Model

TODO Implement Density
- [x] Exponential Density Model
- [ ] Harris-Priester Density Model

TODO Implement Complex Reference Trajectory
- [ ] Simple Nonlinear Model (Drag Only) 
- [ ] Simple Nonlinear Model (Drag + J2-J6)