
thermolib
=========

An open-source library for the calculation of fluid properties.

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

