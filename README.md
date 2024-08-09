
thermolib
=========

An open-source library for the calculation of fluid properties.

# Vdw

| Flash Calculation | Get Corresponding Properties |
| :---: | :---: |
| `t_flash(Ts)` | `T_s()` <br> `p_s()` <br> `rho_v()` <br> `rho_l()` <br> |
| `tp_flash(T,p)` | `T()` <br> `p()` <br> `rho()` <br> |

```rust
use thermolib::Vdw;

let Tc = 430.64; // K
let pc = 7886600.0; // Pa
let M = 0.064064; // kg/mol
let mut SO2 = Vdw::new_fluid(Tc, pc, M);
let _ = SO2.set_molar_unit();

if let Ok(_) = SO2.t_flash(273.15) {
    println!("T_s={}", SO2.T_s().unwrap());
    println!("p_s={}", SO2.p_s().unwrap());
    println!("rho_v={}", SO2.rho_v().unwrap());
    println!("rho_l={}", SO2.rho_l().unwrap());
}

SO2.tp_flash(273.15, 0.1e6);
println!("T={}", SO2.T().unwrap());
println!("p={}", SO2.p().unwrap());
println!("rho={}", SO2.rho().unwrap());

```

```python
from thermolib import Vdw

Tc = 430.64
pc = 7886600
M = 0.064064
SO2 = Vdw(Tc, pc, M)
SO2.set_molar_unit()

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

# Rk

| Flash Calculation | Get Corresponding Properties |
| :---: | :---: |
| `t_flash(Ts)` | `T_s()` <br> `p_s()` <br> `rho_v()` <br> `rho_l()` <br> |
| `tp_flash(T,p)` | `T()` <br> `p()` <br> `rho()` <br> |

```rust
use thermolib::Rk;

let Tc = 430.64; // K
let pc = 7886600.0; // Pa
let M = 0.064064; // kg/mol
let mut SO2 = Rk::new_fluid(Tc, pc, M);
let _ = SO2.set_molar_unit();

```

```python
from thermolib import Rk

Tc = 430.64
pc = 7886600
M = 0.064064
SO2 = Rk(Tc, pc, M)
SO2.set_molar_unit()

```

# Srk

| Flash Calculation | Get Corresponding Properties |
| :---: | :---: |
| `t_flash(Ts)` | `T_s()` <br> `p_s()` <br> `rho_v()` <br> `rho_l()` <br> |
| `tp_flash(T,p)` | `T()` <br> `p()` <br> `rho()` <br> |

```rust
use thermolib::Srk;

let Tc = 430.64; // K
let pc = 7886600.0; // Pa
let omega = 0.256;
let M = 0.064064; // kg/mol
let mut SO2 = Srk::new_fluid(Tc, pc, omega, M);
let _ = SO2.set_molar_unit();

```

```python
from thermolib import Srk

Tc = 430.64
pc = 7886600
omega = 0.256
M = 0.064064
SO2 = Srk(Tc, Pc, omega, M)
SO2.set_molar_unit()

```

# Pr

| Flash Calculation | Get Corresponding Properties |
| :---: | :---: |
| `t_flash(Ts)` | `T_s()` <br> `p_s()` <br> `rho_v()` <br> `rho_l()` <br> |
| `tp_flash(T,p)` | `T()` <br> `p()` <br> `rho()` <br> |

```rust
use thermolib::Pr;

let Tc = 430.64; // K
let pc = 7886600.0; // Pa
let omega = 0.256;
let M = 0.064064; // kg/mol
let mut SO2 = Pr::new_fluid(Tc, pc, omega, M);
let _ = SO2.set_molar_unit();

```

```python
from thermolib import Pr

Tc = 430.64
pc = 7886600
omega = 0.256
M = 0.064064
SO2 = Pr(Tc, Pc, omega, M)
SO2.set_molar_unit()

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
    println!("rho = {}", metal.calc_rho(1800.0).unwrap()); // 2528 kg/m3
    println!("eta = {}", metal.calc_eta(1800.0).unwrap()); // 0.541 mPa*s
    println!("lambda = {}", metal.calc_lambda(1800.0).unwrap()); // 54.88 W/m/K
}

```

```python
from thermolib import LiquidMetal

metal = LiquidMetal("Si")
print("rho =", metal.calc_rho(1800))  # 2528 kg/m3
print("eta =", metal.calc_eta(1800))  # 0.541 mPa*s
print("lambda =", metal.calc_lambda(1800))  # 54.88 W/m/K

```

