
thermolib
=========

An open-source library for the calculation of fluid properties.

# Helmholtz

```rust
use thermolib::Helmholtz;

let mut SO2 = Helmholtz::read_json("SO2.json").expect("no SO2.json");

if let Ok(_) = SO2.tp_flash(273.15, 0.1e6) {
    println!("T={}", SO2.T()); // temperature = 273.15 K
    println!("p={}", SO2.p()); // pressure = 0.1e6 Pa
    println!("rho={}", SO2.rho().unwrap()); // density = 45.093 mol/m3
    println!("cv={}", SO2.cv().unwrap()); // isochoric heat capacity = 31.953 J/mol/K
    println!("cp={}", SO2.cp().unwrap()); // isobaric heat capacity = 41.477 J/mol/K
    println!("w={}", SO2.w().unwrap()); // speed of sound = 209.41 m/s
}

if let Ok(_) = SO2.t_flash(273.15) {
    println!("T={}", SO2.T()); // temperature = 273.15 K
    println!("p={}", SO2.p()); // pressure = 0.15549e6 Pa
    println!("vd={}", SO2.rho_v().unwrap()); // vapor density = 71.106 mol/m3
    println!("ld={}", SO2.rho_l().unwrap()); // liquid density = 22403 mol/m3
}

```

```python
from thermolib import Helmholtz

SO2 = Helmholtz("SO2.json")

SO2.tp_flash(273.15, 0.1e6)
print("T =", SO2.T())  # temperature = 273.15 K
print("p =", SO2.p())  # pressure = 0.1e6 Pa
print("rho =", SO2.rho())  # density = 45.093 mol/m3
print("cv =", SO2.cv())  # isochoric heat capacity = 31.953 J/mol/K
print("cp =", SO2.cp())  # isobaric heat capacity = 41.477 J/mol/K
print("w =", SO2.w())  # speed of sound = 209.41 m/s

SO2.t_flash(273.15)
print("T =", SO2.T())  # temperature = 273.15 K
print("p =", SO2.p())  # pressure = 0.15549e6 Pa
print("vd =", SO2.rho_v())  # vapor density = 71.106 mol/m3
print("ld =", SO2.rho_l())  # liquid density = 22403 mol/m3

```

