# 1D Biot's poroelasticity equation solver using MATLAB

Author: Naren Vohra, Malgorzata Peszynska

The 1D poroelasticity code is used to model deformation and flow in fully saturated porous media using the quasi-static Biot's equations.

## Simulation example

To run example 1:

``` [xn, xcc, t, U, P, Q] = BiotHHMM(10, 10, [0 0 0 0], 1);```

To run example 4:

``` [xn, xcc, t, U, P, Q] = BiotHHMM(100, 10, [1 0 0 1], 4); ```

More examples and the complete documentation can be found in the [pdf]().

## Physical parameters:

| Parameter           | Description                              | Unit                 |
|---------------------|------------------------------------------|----------------------|
| E                   | Young's modulus                          | [MPa]                |
| nu                  | Poisson's ratio                          | [-]                  |
| lam_const, mu_const | Lamé parameter calculated using E and nu | [MPa]                |
| betaf               | Fluid compressibility                    | [1/MPa]              |
| kappa               | Permeability                             | [m<sup>2</sup>]      |
| viscosity           | Fluid viscosity                          | [MPa hr]             |
| alpha               | Biot's constant                          | [-]                  |
| rhof                | Fluid density                             | [kg / m<sup>3</sup>] |
| rhos                | Solid density (eg. of soil grains)       | [kg / m<sup>3</sup>] |
| G                   | Acceleration due to gravity; 1.27290582 $\times 10^8$                 | [m / hr<sup>2</sup>]  |

### Conversion factors used:
1. 1[kg / m<sup>3</sup>] =  = (1/3600 * 1/3600 * 1e-6) [MPa hr / m<sup>2</sup>].
2. 1[m / s<sup>2</sup>] =  = (3600 * 3600 ) [m / hr<sup>2</sup>].

## Acknowledgement 

This research was partially supported by NSF DMS-1912938 "Modeling with Constraints and Phase Transitions in Porous Media" and NSF DMS-1522734 ``Phase transitions in porous media across multiple scales"; PI: Dr. Malgorzata Peszynska. 

## License

We use the [GNU GPL](https://www.gnu.org/licenses/licenses.en.html#GPL) license for our software. 

We ask that you respect the [Creative Commons CC BY-NC-ND 4.0 Attribution-NonCommercial-NoDerivatives 4.0 International license](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode).

See License.md for the full text.


















