
thermolib
=========

An open-source library for the calculation of fluid properties.

# Cubic EOS {Vdw, Rk, Srk, Pr}

| Flash Calculation | Get Corresponding Properties |
| :---: | :---: |
| `t_flash(Ts)` | `T_s()` <br> `p_s()` <br> `rho_v()` <br> `rho_l()` <br> |
| `tp_flash(T,p)` | `T()` <br> `p()` <br> `rho()` <br> |

## Flash Calculation

```rust
use thermolib::{Vdw, Rk, Srk, Pr};
let crit_t = 430.64; // critical temperature of sulfur dioxide // K
let crit_p = 7886600.0; // critical pressure of sulfur dioxide // Pa
let acentric_factor = 0.256;
let mut fluid_vdw = Vdw::new_fluid(crit_t, crit_p);
let mut fluid_rk = Rk::new_fluid(crit_t, crit_p);
let mut fluid_srk = Srk::new_fluid(crit_t, crit_p, acentric_factor);
let mut fluid_pr = Pr::new_fluid(crit_t, crit_p, acentric_factor);
```

```python
from thermolib import Vdw, Rk, Srk, Pr
TC = 430.64
PC = 7886600
OMEGA = 0.256
fluid_vdw = Vdw(TC, PC)
fluid_rk = Rk(TC, PC)
fluid_srk = Srk(TC, PC, OMEGA)
fluid_pr = Pr(TC, PC, OMEGA)
```

## Get Corresponding Properties

```rust
if let Ok(_) = fluid.t_flash(273.15) {
    println!("T_s={}", fluid.T_s().unwrap());
    println!("p_s={}", fluid.p_s().unwrap());
    println!("rho_v={}", fluid.rho_v().unwrap());
    println!("rho_l={}", fluid.rho_l().unwrap());
}
fluid.tp_flash(273.15, 0.1e6);
println!("T={}", fluid.T().unwrap());
println!("p={}", fluid.p().unwrap());
println!("rho={}", fluid.rho().unwrap());
```

```python
fluid.t_flash(273.15)
print("T_s =", fluid.T_s())
print("p_s =", fluid.p_s())
print("rho_v =", fluid.rho_v())
print("rho_l =", fluid.rho_l())
SO2.tp_flash(273.15, 0.1e6)
print("T =", fluid.T())
print("p =", fluid.p())
print("rho =", fluid.rho())
```

# Helmholtz

| Flash Calculation | Get Corresponding Properties |
| :---: | :---: |
| `t_flash(Ts)` | `T_s()` <br> `p_s()` <br> `rho_v()` <br> `rho_l()` <br> |
| `tp_flash(T,p)` | `T()` <br> `p()` <br> `rho()` <br> `cv()` <br> `cp()` <br> `w()` <br> |

```rust
use thermolib::Helmholtz;

let mut SO2 = Helmholtz::read_json("SO2.json").expect("no SO2.json");

if let Ok(_) = SO2.t_flash(273.15) {
    println!("T_s={}", SO2.T_s().unwrap()); // temperature = 273.15 K
    println!("p_s={}", SO2.p_s().unwrap()); // pressure = 0.15549e6 Pa
    println!("rho_v={}", SO2.rho_v().unwrap()); // vapor density = 71.106 mol/m3
    println!("rho_l={}", SO2.rho_l().unwrap()); // liquid density = 22403 mol/m3
}

if let Ok(_) = SO2.tp_flash(273.15, 0.1e6) {
    println!("T={}", SO2.T().unwrap()); // temperature = 273.15 K
    println!("p={}", SO2.p().unwrap()); // pressure = 0.1e6 Pa
    println!("rho={}", SO2.rho().unwrap()); // density = 45.093 mol/m3
    println!("cv={}", SO2.cv().unwrap()); // isochoric heat capacity = 31.953 J/mol/K
    println!("cp={}", SO2.cp().unwrap()); // isobaric heat capacity = 41.477 J/mol/K
    println!("w={}", SO2.w().unwrap()); // speed of sound = 209.41 m/s
}

```

```python
from thermolib import Helmholtz

SO2 = Helmholtz("SO2.json")

SO2.t_flash(273.15)
print("T_s =", SO2.T_s())  # temperature = 273.15 K
print("p_s =", SO2.p_s())  # pressure = 0.15549e6 Pa
print("rho_v =", SO2.rho_v())  # vapor density = 71.106 mol/m3
print("rho_l =", SO2.rho_l())  # liquid density = 22403 mol/m3

SO2.tp_flash(273.15, 0.1e6)
print("T =", SO2.T())  # temperature = 273.15 K
print("p =", SO2.p())  # pressure = 0.1e6 Pa
print("rho =", SO2.rho())  # density = 45.093 mol/m3
print("cv =", SO2.cv())  # isochoric heat capacity = 31.953 J/mol/K
print("cp =", SO2.cp())  # isobaric heat capacity = 41.477 J/mol/K
print("w =", SO2.w())  # speed of sound = 209.41 m/s

```

# LiquidMetal

| Function | Unit |
| :---: | :---: |
| `calc_rho(T)` | `kg/m^3` |
| `calc_eta(T)` | `mPa*s` |
| `calc_lambda(T)` | `W/m/K` |

```rust
use thermolib::LiquidMetal;

if let Ok(metal) = LiquidMetal::new_metal("Si") {
    println!("rho = {}", metal.calc_rho(1800.0).unwrap()); // 2520 kg/m3
    println!("eta = {}", metal.calc_eta(1800.0).unwrap()); // 0.541 mPa*s
    println!("lambda = {}", metal.calc_lambda(1800.0).unwrap()); // 54.88 W/m/K
}

```

```python
from thermolib import LiquidMetal

metal = LiquidMetal("Si")
print("rho =", metal.calc_rho(1800))  # 2520 kg/m3
print("eta =", metal.calc_eta(1800))  # 0.541 mPa*s
print("lambda =", metal.calc_lambda(1800))  # 54.88 W/m/K

```

# PcSaftPure

| Flash Calculation | Get Corresponding Properties |
| :---: | :---: |
| `c_flash()` | `T()` <br> `p()` <br> `rho()` <br> |
| `t_flash(Ts)` | `T_s()` <br> `p_s()` <br> `rho_v()` <br> `rho_l()` <br> |
| `tp_flash(T,p)` | `T()` <br> `p()` <br> `rho()` <br> |

```rust
use thermolib::PcSaftPure;

let m = 2.8611;
let sigma = 2.6826;
let epsilon = 205.35;
let mut SO2 = PcSaftPure::new_fluid(m, sigma, epsilon);

if let Ok(_) = SO2.c_flash() {
    println!("T_c={}", SO2.T().unwrap());
    println!("p_c={}", SO2.p().unwrap());
    println!("rho_c={}", SO2.rho().unwrap());
}

if let Ok(_) = SO2.t_flash(273.15) {
    println!("T_s={}", SO2.T_s().unwrap());
    println!("p_s={}", SO2.p_s().unwrap());
    println!("rho_v={}", SO2.rho_v().unwrap());
    println!("rho_l={}", SO2.rho_l().unwrap());
}

if let Ok(_) = SO2.tp_flash(273.15, 0.1e6) {
    println!("T={}", SO2.T().unwrap());
    println!("p={}", SO2.p().unwrap());
    println!("rho={}", SO2.rho().unwrap());
}

```

```python
from thermolib import PcSaftPure

M = 2.8611  # m
S = 2.6826  # sigma
E = 205.35  # epsilon
SO2 = PcSaftPure(M, S, E)

SO2.c_flash()
print("T_c =", SO2.T())
print("p_c =", SO2.p())
print("rho_c =", SO2.rho())

SO2.t_flash(273.15)
print("T_s =", SO2.T_s())
print("p_s =", SO2.p_s())
print("rho_v =", SO2.rho_v())
print("rho_l =", SO2.rho_l())

SO2.tp_flash(273.15, 0.1e6)
print("T =", SO2.T())
print("p =", SO2.p())
print("rho =", SO2.rho())

```

