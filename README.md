# DPCTutorial
Yalmip implementation of data-driven predictive controllers. 

Requires the Yalmip toolbox for Matlab (https://yalmip.github.io/) and a compatible QP solver (Matlab Quadprog will do).

Matlab files that build SPC, DeePC and offset-free SPC and DeePC, respectively, from data generated by a linear discrete-time system.

The PADC file generates the system dynamics and constraints for a power amplifier example. 

If you define your own linear system you do not need to call the PADC function.

The Hanke file plots Hanke-Raus and Tikhonov heuristics for tuning regularized DeePC or iDeePC in time scale.

The DeePCRegularizationPlot file is a function that can be called to plot Hanke-Raus and Tikhonov heursitics in time and logarithmic scales. 

Comments are provided for each data-driven predictive controller implementation on how to modify data generation, noise, disturbance and tuning.

If you use DPCTutorial implementations in your research, please cite the following paper (to appear in IEEE Proceedings after the CDC 2022 in December):

```bibtex
@inproceedings{DPCTutorial2022,
	author      = {Lazar, M. and Verheijen, P. C. N.},
	booktitle   = {IEEE Conference on Decision and Control},
	title       = {{Offset-free Data-Driven Predictive Control}},
	year        = {2022},
	address     = {Cancun, Mexic},
	month       = {Dec.},
}
```
