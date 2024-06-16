# Advection

## 1-Dim Navier Stokes equation. caseA form.

$u$ is the velocity of the gas, $p$ is the pressure, $\rho$ is the molar density, $R$ is the ideal gas constant, $T$ is the absolute temperature.

### Navier Stokes equation's caseA form.

$$
\frac{\partial \left( \rho u \right)}{\partial t} + \frac{\partial \left( \rho u^2 \right)}{\partial x} = - \frac{\partial p}{\partial x} + \mu \frac{\partial^2 u}{\partial x^2} \ \ \rm{(caseA).}
$$

continuity equation.

$$
\frac{\partial \rho}{\partial t} + \frac{\partial \left( \rho u \right)}{\partial x} = 0.
$$

Equation of states.

$$
p = \rho R T.
$$

In this repository, temperature and molar density is assumed to be a constant.

## $q$-$\rho$ representation.

$$
q := \rho u
$$

$$
\frac{\partial q}{\partial t} + \frac{\partial}{\partial x}\left( \frac{q^2}{\rho} \right) = - RT \frac{\partial \rho}{\partial x} + \mu \frac{\partial^2}{\partial x^2} \left( \frac{q}{\rho} \right)
$$

### continuity equation. caseA

$$
\frac{\partial \rho}{\partial t} + \frac{\partial q}{\partial x} = 0.
$$

### continuity equation. caseB

$$
\frac{\partial \rho}{\partial t} + \frac{q}{\rho} \frac{\partial \rho}{\partial x} + \rho \frac{\partial}{\partial x}\left( \frac{q}{\rho} \right) = 0.
$$

### Navier Stokes equation's caseB form.

Assuming the continuity equation, the Navier-Stokes equation can be transformed as follows.

$$
\rho \left(\frac{\partial u}{\partial t} + u \frac{\partial u }{\partial x} \right)= - \frac{\partial p}{\partial x} + \mu \frac{\partial^2 u}{\partial x^2} \ \ \rm{(caseB).}
$$
