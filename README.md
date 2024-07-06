# Code for: Interpretable data-driven turbulence modeling for separated flows using symbolic regression with unit constraints

Machine learning techniques have been applied to enhance turbulence modeling in recent years. However, the "black box" nature of most machine learning techniques poses significant interpretability challenges in improving turbulence models. This paper introduces a novel unit-constrained turbulence modeling framework using symbolic regression to overcome these challenges. The framework amends the constitutive equation of linear eddy viscosity models (LEVMs) by establishing explicit equations between the Reynolds stress deviation and mean flow quantities, thereby improving the LEVM model's predictive capability for large separated turbulence. Unit consistency constraints are applied to the symbolic expressions to ensure physical realizability. The framework’s effectiveness is demonstrated through its application to the separated flow over 2D periodic hills. Numerical simulations are conducted to validate the accuracy and generalization capabilities of the learned turbulence model. Compared to the standard k-ε model, the learned model exhibits significantly enhanced predictive accuracy for anisotropic Reynolds stresses and flow velocities. It also more accurately identifies flow separations and exhibits promising generalization capabilities across different flow geometries.

A sketch map of the symbolic regression turbulence modeling framework with unit constraint: 
![image](https://github.com/Jamie-Z/PhySO-RANS/assets/72621861/60dabd06-0b48-41f0-b182-fb3c4f920a9e)


Requirements: 
- Ubuntu 20.04
- OpenFOAM 7
