# anesthesia-control-system
**Summer 2025**

This is an NSF-funded research project done under the mentorship of Dr. Hosam Fathy and Dr. Jin-Oh Hahn at the University of Maryland during my sophomore summer. 

I apply a novel "graceful" control barrier function (CBF) first introduced by Yejin Moon to the problem of automatic sedation, but with two barriers in a mirrored fashion, as sedation has both safety-critical upper and lower bounds. I also apply a deterministic interval observer that provides conservative bounds to the controller.  

This approach poses two novel benefits: (1) The CBF is designed to apply to systems of degree 2, or those systems in which the input variable $u(t)$ appears only in the derivative of the output. (2) It provides multi-layered safety regimes that mimics induction/maintenance phases in the existing literature and enforces more aggressive control near the barrier. 

The code is organized as follows: 
- `Paper.pdf` is the paper, to be published in the Modeling, Estimation and Control Conference (MECC 2026). 
- `IntervalObserver.m` computes the interval observer gains in the diagonalized basis once for all simulations. The results were manually copied into `PropofolSystem.m`. The interval observer, as explained in `Paper.pdf`, requires $\mathbf{A} - \mathbf{LC}$ to be both Hurwitz and Metzler. It turns out that this latter constraint quite a numerical thorn to this system, so I implement the observer in the diagonalized basis, where it will be trivially Metzler, and transform back to conservative bounds in the original basis. 
- `PropofolSystem.m` is a module containing functions to simulate the model, the controller, and the observer, as well as to plot results. 
- `runPropofolSystem.m` is the master-runner that feeds parameters and settings into `PropofolSystem.m` to simulate and plot various initial conditions. 
