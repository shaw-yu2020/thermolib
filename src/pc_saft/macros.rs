/// macro_rules! fn_vec
macro_rules! fn_vec {
    ($name:ty) => {
        #[cfg_attr(feature = "with_pyo3", pymethods)]
        impl $name {
            pub fn vec_t_flash_v(&mut self, temp: Vec<f64>) -> Vec<f64> {
                temp.into_iter()
                    .map(|t| {
                        if self.t_flash(t).is_err() {
                            println!("t_flash_v diverge in {} K", t);
                            f64::NAN
                        } else {
                            self.rhov_num
                        }
                    })
                    .collect()
            }
            pub fn vec_t_flash_l(&mut self, temp: Vec<f64>) -> Vec<f64> {
                temp.into_iter()
                    .map(|t| {
                        if self.t_flash(t).is_err() {
                            println!("t_flash_l diverge in {} K", t);
                            f64::NAN
                        } else {
                            self.rhol_num
                        }
                    })
                    .collect()
            }
            pub fn vec_p(&mut self, temp: Vec<f64>, dens_num: Vec<f64>) -> Vec<f64> {
                temp.into_iter()
                    .zip(dens_num)
                    .map(|(t, d)| self.calc_p(t, d))
                    .collect()
            }
            pub fn vec_tp_flash(&mut self, temp: Vec<f64>, pres: Vec<f64>) -> Vec<f64> {
                temp.into_iter()
                    .zip(pres)
                    .map(|(t, p)| {
                        if self.tp_flash(t, p).is_err() {
                            println!("tp_flash diverge in {} K {} Pa", t, p);
                            f64::NAN
                        } else {
                            self.rho_num
                        }
                    })
                    .collect()
            }
            pub fn vec_tp_flash_v(&mut self, temp: Vec<f64>, pres: Vec<f64>) -> Vec<f64> {
                temp.into_iter()
                    .zip(pres)
                    .map(|(t, p)| {
                        let d3 =
                            self.sigma3 * (1.0 - 0.12 * (-3.0 * self.epsilon / t).exp()).powi(3);
                        // Iteration from vapor phase: eta = 1E-10
                        self.calc_density(t, p, 1E-10 / (FRAC_PI_6 * self.m * d3))
                    })
                    .collect()
            }
            pub fn vec_tp_flash_l(&mut self, temp: Vec<f64>, pres: Vec<f64>) -> Vec<f64> {
                temp.into_iter()
                    .zip(pres)
                    .map(|(t, p)| {
                        let d3 =
                            self.sigma3 * (1.0 - 0.12 * (-3.0 * self.epsilon / t).exp()).powi(3);
                        // Iteration from liquid phase: eta = 0.5
                        self.calc_density(t, p, 0.5 / (FRAC_PI_6 * self.m * d3))
                    })
                    .collect()
            }
            pub fn vec_cv(&mut self, temp: Vec<f64>, dens_num: Vec<f64>) -> Vec<f64> {
                temp.into_iter()
                    .zip(dens_num)
                    .map(|(t, d)| self.calc_cv(t, d))
                    .collect()
            }
            pub fn vec_cp(&mut self, temp: Vec<f64>, dens_num: Vec<f64>) -> Vec<f64> {
                temp.into_iter()
                    .zip(dens_num)
                    .map(|(t, d)| self.calc_cp(t, d))
                    .collect()
            }
            pub fn vec_w(
                &mut self,
                temp: Vec<f64>,
                dens_num: Vec<f64>,
                molar_mass: f64,
            ) -> Vec<f64> {
                temp.into_iter()
                    .zip(dens_num)
                    .map(|(t, d)| (self.calc_w2(t, d) / molar_mass).sqrt()) // m/s
                    .collect()
            }
        }
    };
}
/// macro_rules! fn_c_flash
macro_rules! fn_c_flash {
    ($name:ty) => {
        #[cfg_attr(feature = "with_pyo3", pymethods)]
        impl $name {
            pub fn c_flash(&mut self) -> anyhow::Result<()> {
                // Iteration from temp_c = 1000 dens_c = 1E-3
                let (mut temp_c, mut dens_c) = (1E3, 1E-3);
                // Define variables
                let mut p_t0d1 = self.calc_p_t0d1(temp_c, dens_c);
                let mut p_t0d2 = self.calc_p_t0d2(temp_c, dens_c);
                let (mut p_t1d1, mut p_t1d2, mut p_t0d3);
                for _i in 1..100000 {
                    p_t1d1 = self.calc_p_t1d1(temp_c, dens_c);
                    p_t1d2 = self.calc_p_t1d2(temp_c, dens_c);
                    p_t0d3 = self.calc_p_t0d3(temp_c, dens_c);
                    temp_c -=
                        (p_t0d1 * p_t0d3 - p_t0d2 * p_t0d2) / (p_t1d1 * p_t0d3 - p_t1d2 * p_t0d2);
                    dens_c -=
                        (p_t0d1 * p_t1d2 - p_t0d2 * p_t1d1) / (p_t0d2 * p_t1d2 - p_t0d3 * p_t1d1);
                    p_t0d1 = self.calc_p_t0d1(temp_c, dens_c);
                    p_t0d2 = self.calc_p_t0d2(temp_c, dens_c);
                    if p_t0d1.abs() < 1E3 && p_t0d2.abs() < 1E6 {
                        self.is_single_phase = true;
                        return Ok(());
                    }
                }
                Err(anyhow!(PcSaftErr::NotConvForC))
            }
            fn calc_p_t1d1(&mut self, temp: f64, rho_num: f64) -> f64 {
                (FRAC_RE30_NA * temp)
                    * (1.0
                        + 2.0 * self.r_t0d1(temp, rho_num)
                        + self.r_t0d2(temp, rho_num)
                        + 2.0 * self.r_t1d1(temp, rho_num)
                        + self.r_t1d2(temp, rho_num))
            }
            fn calc_p_t1d2(&mut self, temp: f64, rho_num: f64) -> f64 {
                FRAC_RE30_NA * temp / rho_num
                    * (2.0 * self.r_t0d1(temp, rho_num)
                        + 4.0 * self.r_t0d2(temp, rho_num)
                        + self.r_t0d3(temp, rho_num)
                        + 2.0 * self.r_t1d1(temp, rho_num)
                        + 4.0 * self.r_t1d2(temp, rho_num)
                        + self.r_t1d3(temp, rho_num))
            }
            fn calc_p_t0d3(&mut self, temp: f64, rho_num: f64) -> f64 {
                FRAC_RE30_NA * temp / rho_num.powi(2)
                    * (6.0 * self.r_t0d2(temp, rho_num)
                        + 6.0 * self.r_t0d3(temp, rho_num)
                        + self.r_t0d4(temp, rho_num))
            }
        }
    };
}
/// macro_rules! fn_flash
macro_rules! fn_flash {
    ($name:ty) => {
        #[cfg_attr(feature = "with_pyo3", pymethods)]
        impl $name {
            pub fn t_flash(&mut self, temp: f64) -> anyhow::Result<()> {
                let d3 = self.sigma3 * (1.0 - 0.12 * (-3.0 * self.epsilon / temp).exp()).powi(3);
                // Vapor phase: eta = 1E-10
                let rhov_num_guess = 1E-10 / (FRAC_PI_6 * self.m * d3);
                let mut rhov_num = rhov_num_guess;
                let mut p_t0d1_v = self.calc_p_t0d1(temp, rhov_num);
                for _i in 1..100000 {
                    if p_t0d1_v.abs() < 10.0 {
                        break;
                    } else {
                        rhov_num -= p_t0d1_v / self.calc_p_t0d2(temp, rhov_num);
                        p_t0d1_v = self.calc_p_t0d1(temp, rhov_num);
                    }
                }
                let pv_limit = self.calc_p(temp, rhov_num);
                if pv_limit.is_sign_negative() {
                    return Err(anyhow!(PcSaftErr::NotConvForT));
                }
                // Liquid phase: eta = 0.5
                let rhol_num_guess = 0.5 / (FRAC_PI_6 * self.m * d3);
                let mut rhol_num = rhol_num_guess;
                let mut p_t0d1_l = self.calc_p_t0d1(temp, rhol_num);
                for _i in 1..100000 {
                    if p_t0d1_l.abs() < 10.0 {
                        break;
                    } else {
                        rhol_num -= p_t0d1_l / self.calc_p_t0d2(temp, rhol_num);
                        p_t0d1_l = self.calc_p_t0d1(temp, rhol_num);
                    }
                }
                let pl_limit = self.calc_p(temp, rhol_num);
                if pl_limit > pv_limit {
                    return Err(anyhow!(PcSaftErr::NotConvForT));
                }
                // Iteration for saturation state
                if rhol_num < rhov_num {
                    return Err(anyhow!(PcSaftErr::NotConvForT));
                }
                let rhov_max = rhov_num;
                let lnphi_diff = |p: f64| {
                    rhov_num = self.calc_density(temp, p, rhov_num_guess);
                    if rhov_num.is_nan() && (p / pv_limit - 1.0).abs() < 1E-2 {
                        rhov_num = rhov_max;
                    };
                    rhol_num = self.calc_density(temp, p, rhol_num_guess);
                    self.calc_lnphi(temp, rhov_num) - self.calc_lnphi(temp, rhol_num)
                };
                let ps = brent_zero(lnphi_diff, pv_limit - 1.0, pl_limit.max(1.0));
                if ps.is_nan() {
                    Err(anyhow!(PcSaftErr::NotConvForT))
                } else {
                    self.is_single_phase = false;
                    self.temp = temp;
                    self.rhov_num = rhov_num;
                    self.rhol_num = rhol_num;
                    Ok(())
                }
            }
            pub fn tp_flash(&mut self, temp: f64, pres: f64) -> anyhow::Result<()> {
                let d3 = self.sigma3 * (1.0 - 0.12 * (-3.0 * self.epsilon / temp).exp()).powi(3);
                // Iteration from vapor phase: eta = 1E-10
                let rhov_num = self.calc_density(temp, pres, 1E-10 / (FRAC_PI_6 * self.m * d3));
                let lnphi_v = if rhov_num.is_nan() {
                    f64::INFINITY
                } else {
                    self.calc_lnphi(temp, rhov_num)
                };
                // Iteration from liquid phase: eta = 0.5
                let rhol_num = self.calc_density(temp, pres, 0.5 / (FRAC_PI_6 * self.m * d3));
                let lnphi_l = if rhol_num.is_nan() {
                    f64::INFINITY
                } else {
                    self.calc_lnphi(temp, rhol_num)
                };
                // Select the correct output
                if lnphi_v.is_infinite() && lnphi_l.is_infinite() {
                    Err(anyhow!(PcSaftErr::NotConvForTP))
                } else if lnphi_v.is_infinite() {
                    self.eta_flash(temp, rhol_num);
                    self.is_single_phase = true;
                    Ok(())
                } else if lnphi_l.is_infinite() {
                    self.eta_flash(temp, rhov_num);
                    self.is_single_phase = true;
                    Ok(())
                } else {
                    if lnphi_v < lnphi_l {
                        self.eta_flash(temp, rhov_num);
                    } else {
                        self.eta_flash(temp, rhol_num);
                    }
                    self.is_single_phase = true;
                    Ok(())
                }
            }
        }
        impl $name {
            fn calc_density(&mut self, temp: f64, p: f64, rho_num_guess: f64) -> f64 {
                let mut rho_num = rho_num_guess;
                let (mut p_diff, mut val_p_t0d1, mut rho_num_diff);
                for _i in 1..10000 {
                    p_diff = self.calc_p(temp, rho_num) - p;
                    if p_diff.abs() < f64::EPSILON {
                        return rho_num;
                    }
                    val_p_t0d1 = self.calc_p_t0d1(temp, rho_num);
                    if val_p_t0d1.is_sign_negative() {
                        return f64::NAN;
                    }
                    rho_num_diff = p_diff / val_p_t0d1;
                    if rho_num_diff.abs() < f64::EPSILON {
                        return rho_num;
                    }
                    rho_num -= rho_num_diff;
                    if rho_num.is_sign_negative() {
                        return f64::NAN;
                    }
                }
                f64::NAN
            }
            fn calc_lnphi(&mut self, temp: f64, rho_num: f64) -> f64 {
                self.r_t0d0(temp, rho_num) + self.r_t0d1(temp, rho_num)
                    - (1.0 + self.r_t0d1(temp, rho_num)).ln()
            }
        }
    };
}
/// macro_rules! fn_single_prop
macro_rules! fn_single_prop {
    ($name:ty) => {
        #[cfg_attr(feature = "with_pyo3", pymethods)]
        impl $name {
            pub fn print_derivatives(&mut self) {
                self.check_derivatives(true);
            }
            #[allow(non_snake_case)]
            pub fn T(&self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Ok(self.temp)
                } else {
                    Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
                }
            }
            pub fn rho(&self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Ok(self.rho_num / FRAC_NA_1E30)
                } else {
                    Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
                }
            }
            pub fn p(&mut self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Ok(self.calc_p(self.temp, self.rho_num))
                } else {
                    Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
                }
            }
            pub fn s_res(&mut self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Ok(self.calc_s_res(self.temp, self.rho_num))
                } else {
                    Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
                }
            }
            pub fn u_res(&mut self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Ok(self.calc_u_res(self.temp, self.rho_num))
                } else {
                    Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
                }
            }
            pub fn h_res(&mut self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Ok(self.calc_h_res(self.temp, self.rho_num))
                } else {
                    Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
                }
            }
            pub fn cv_res(&mut self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Ok(self.calc_cv_res(self.temp, self.rho_num))
                } else {
                    Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
                }
            }
            pub fn cp_res(&mut self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Ok(self.calc_cp_res(self.temp, self.rho_num))
                } else {
                    Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
                }
            }
        }
    };
}
/// macro_rules! fn_double_prop
macro_rules! fn_double_prop {
    ($name:ty) => {
        #[cfg_attr(feature = "with_pyo3", pymethods)]
        impl $name {
            pub fn p_s(&mut self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Err(anyhow!(PcSaftErr::OnlyInDoublePhase))
                } else {
                    Ok(self.calc_p(self.temp, self.rhov_num) / 2.0
                        + self.calc_p(self.temp, self.rhol_num) / 2.0)
                }
            }
            #[allow(non_snake_case)]
            pub fn T_s(&self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Err(anyhow!(PcSaftErr::OnlyInDoublePhase))
                } else {
                    Ok(self.temp)
                }
            }
            pub fn rho_v(&self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Err(anyhow!(PcSaftErr::OnlyInDoublePhase))
                } else {
                    Ok(self.rhov_num / FRAC_NA_1E30)
                }
            }
            pub fn rho_l(&self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Err(anyhow!(PcSaftErr::OnlyInDoublePhase))
                } else {
                    Ok(self.rhol_num / FRAC_NA_1E30)
                }
            }
        }
    };
}
/// macro_rules! fn_virial_prop
macro_rules! fn_virial_prop {
    ($name:ty) => {
        #[cfg_attr(feature = "with_pyo3", pymethods)]
        #[allow(non_snake_case)] // For pymethods
        impl $name {
            pub fn B(&mut self, temp: f64) -> f64 {
                let mut rho_num: f64 = 1E-9;
                let mut val_old: f64 = self.r_t0d1(temp, rho_num) / rho_num;
                let mut val_new: f64 = 0.0;
                loop {
                    rho_num /= 10.0;
                    if rho_num < 1E-30 {
                        break;
                    }
                    val_new = self.r_t0d1(temp, rho_num) / rho_num;
                    if (val_new / val_old - 1.0).abs() < 1E-9 {
                        break;
                    }
                    val_old = val_new;
                }
                val_new * FRAC_NA_1E30
            }
            pub fn C(&mut self, temp: f64) -> f64 {
                let mut rho_num: f64 = 1E-9;
                let mut val_old: f64 = self.r_t0d2(temp, rho_num) / rho_num.powi(2);
                let mut val_new: f64 = 0.0;
                loop {
                    rho_num /= 10.0;
                    if rho_num < 1E-30 {
                        break;
                    }
                    val_new = self.r_t0d2(temp, rho_num) / rho_num.powi(2);
                    if (val_new / val_old - 1.0).abs() < 1E-9 {
                        break;
                    }
                    val_old = val_new;
                }
                val_new * FRAC_NA_1E30.powi(2)
            }
            pub fn D(&mut self, temp: f64) -> f64 {
                let mut rho_num: f64 = 1E-9;
                let mut val_old: f64 = self.r_t0d3(temp, rho_num) / rho_num.powi(3);
                let mut val_new: f64 = 0.0;
                loop {
                    rho_num /= 10.0;
                    if rho_num < 1E-30 {
                        break;
                    }
                    val_new = self.r_t0d3(temp, rho_num) / rho_num.powi(3);
                    if (val_new / val_old - 1.0).abs() < 1E-9 {
                        break;
                    }
                    val_old = val_new;
                }
                val_new * FRAC_NA_1E30.powi(3) / 2.0
            }
        }
    };
}
/// macro_rules! fn_aly_lee_cp0
macro_rules! fn_aly_lee_cp0 {
    ($name:ty) => {
        #[cfg_attr(feature = "with_pyo3", pymethods)]
        impl $name {
            #[allow(non_snake_case)]
            pub fn set_aly_lee_cp0(&mut self, B: f64, C: f64, D: f64, E: f64, F: f64) {
                self.cv_B = B / R - 1.0; // cv_B = B/R -1
                self.cv_C = C / R; // cv_C = C/R
                self.cv_D = D; // cv_D = D
                self.cv_E = E / R; // cv_E = E/R
                self.cv_F = F; // cv_F = F
            }
            pub fn cv(&mut self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Ok(self.calc_cv(self.temp, self.rho_num))
                } else {
                    Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
                }
            }
            pub fn cp(&mut self) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Ok(self.calc_cp(self.temp, self.rho_num))
                } else {
                    Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
                }
            }
            pub fn w(&mut self, molar_mass: f64) -> anyhow::Result<f64> {
                if self.is_single_phase {
                    Ok((self.calc_w2(self.temp, self.rho_num) / molar_mass).sqrt())
                } else {
                    Err(anyhow!(PcSaftErr::OnlyInSinglePhase))
                }
            }
        }
        impl $name {
            fn calc_ideal_cv(&mut self, temp: f64) -> f64 {
                self.cv_B
                    + self.cv_C * (self.cv_D / temp / (self.cv_D / temp).sinh()).powi(2)
                    + self.cv_E * (self.cv_F / temp / (self.cv_F / temp).cosh()).powi(2)
            }
            fn calc_cv(&mut self, temp: f64, rho_num: f64) -> f64 {
                R * self.calc_ideal_cv(temp) + self.calc_cv_res(temp, rho_num)
            }
            fn calc_cp(&mut self, temp: f64, rho_num: f64) -> f64 {
                R * (self.calc_ideal_cv(temp) + 1.0) + self.calc_cp_res(temp, rho_num)
            }
            fn calc_w2(&mut self, temp: f64, rho_num: f64) -> f64 {
                R * temp
                    * ((1.0 + 2.0 * self.r_t0d1(temp, rho_num) + self.r_t0d2(temp, rho_num))
                        + (1.0 + self.r_t0d1(temp, rho_num) + self.r_t1d1(temp, rho_num)).powi(2)
                            / (self.calc_ideal_cv(temp)
                                - 2.0 * self.r_t1d0(temp, rho_num)
                                - self.r_t2d0(temp, rho_num)))
            }
        }
    };
}
/// macro_rules! fn_calc_prop
macro_rules! fn_calc_prop {
    () => {
        fn calc_p(&mut self, temp: f64, rho_num: f64) -> f64 {
            FRAC_RE30_NA * temp * rho_num * (1.0 + self.r_t0d1(temp, rho_num))
        }
        fn calc_p_t0d1(&mut self, temp: f64, rho_num: f64) -> f64 {
            (FRAC_RE30_NA * temp)
                * (1.0 + 2.0 * self.r_t0d1(temp, rho_num) + self.r_t0d2(temp, rho_num))
        }
        fn calc_p_t0d2(&mut self, temp: f64, rho_num: f64) -> f64 {
            FRAC_RE30_NA * temp / rho_num
                * (2.0 * self.r_t0d1(temp, rho_num)
                    + 4.0 * self.r_t0d2(temp, rho_num)
                    + self.r_t0d3(temp, rho_num))
        }
        fn calc_cv_res(&mut self, temp: f64, rho_num: f64) -> f64 {
            -R * (2.0 * self.r_t1d0(temp, rho_num) + self.r_t2d0(temp, rho_num))
        }
        fn calc_cp_res(&mut self, temp: f64, rho_num: f64) -> f64 {
            R * ((-1.0 - 2.0 * self.r_t1d0(temp, rho_num) - self.r_t2d0(temp, rho_num))
                + (1.0 + self.r_t0d1(temp, rho_num) + self.r_t1d1(temp, rho_num)).powi(2)
                    / (1.0 + 2.0 * self.r_t0d1(temp, rho_num) + self.r_t0d2(temp, rho_num)))
        }
        fn calc_s_res(&mut self, temp: f64, rho_num: f64) -> f64 {
            R * temp * (self.r_t0d0(temp, rho_num) + self.r_t1d0(temp, rho_num))
        }
        fn calc_u_res(&mut self, temp: f64, rho_num: f64) -> f64 {
            -R * temp * self.r_t1d0(temp, rho_num)
        }
        fn calc_h_res(&mut self, temp: f64, rho_num: f64) -> f64 {
            R * temp * (-self.r_t1d0(temp, rho_num) + self.r_t0d1(temp, rho_num))
        }
    };
}
/// macro_rules! fn_check_derivatives
macro_rules! fn_check_derivatives {
    () => {
        pub fn check_derivatives(&mut self, print_val: bool) {
            let (t, d) = (self.temp, self.rho_num);
            if print_val {
                println!("[t0d0 == t0d0] t0d0 ={}", self.r_t0d0(t, d));
            }
            let compare_val = |val_calc: f64, val_diff: f64| {
                assert_eq!(
                    format!("{:.8e}", &val_calc.abs()),
                    format!("{:.8e}", &val_diff.abs())
                )
            };
            // derivative for density
            let val_calc = self.r_t0d1(t, d) / d;
            let val_diff = romberg_diff(|dx: f64| self.r_t0d0(t, dx), d);
            if print_val {
                println!("[t0d1 => t0d1] t0d1 ={:e}", val_calc);
                println!("[t0d0 -> t0d1] diff ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
            // derivative for density+density
            let val_calc = self.r_t0d2(t, d) / d.powi(2);
            let val_diff = romberg_diff(|dx: f64| self.r_t0d1(t, dx) / dx, d);
            if print_val {
                println!("[t0d2 == t0d2] t0d2 ={:e}", val_calc);
                println!("[t0d1 -> t0d2] diff ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
            // derivative for density+density+density
            let val_calc = self.r_t0d3(t, d) / d.powi(3);
            let val_diff = romberg_diff(|dx: f64| self.r_t0d2(t, dx) / dx.powi(2), d);
            if print_val {
                println!("[t0d3 == t0d3] t0d3 ={:e}", val_calc);
                println!("[t0d2 -> t0d3] diff ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
            // derivative for density+density+density+density
            let val_calc = self.r_t0d4(t, d) / d.powi(4);
            let val_diff = romberg_diff(|dx: f64| self.r_t0d3(t, dx) / dx.powi(3), d);
            if print_val {
                println!("[t0d4 == t0d4] t0d4 ={:e}", val_calc);
                println!("[t0d3 -> t0d4] diff ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
            // derivative for temperature
            let val_calc = self.r_t1d0(t, d) / t;
            let val_diff = romberg_diff(|tx: f64| self.r_t0d0(tx, d), t);
            if print_val {
                println!("[t1d0 == t1d0] t1d0 ={:e}", val_calc);
                println!("[t0d0 -> t1d0] diff ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
            // derivative for temperature+density
            let val_calc = self.r_t1d1(t, d) / t / d;
            let val_diff = romberg_diff(|dx: f64| self.r_t1d0(t, dx) / t, d);
            if print_val {
                println!("[t1d1 == t1d1] t1d1 ={:e}", val_calc);
                println!("[t1d0 -> t1d1] diff ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
            let val_diff = romberg_diff(|tx: f64| self.r_t0d1(tx, d) / d, t);
            if print_val {
                println!("[t1d1 == t1d1] t1d1 ={:e}", val_calc);
                println!("[t0d1 -> t1d1] diff ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
            // derivative for temperature+density+density
            let val_calc = self.r_t1d2(t, d) / t / d.powi(2);
            let val_diff = romberg_diff(|dx: f64| self.r_t1d1(t, dx) / t / dx, d);
            if print_val {
                println!("[t1d2 == t1d2] t1d2 ={:e}", val_calc);
                println!("[t1d1 -> t1d2] diff ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
            let val_diff = romberg_diff(|tx: f64| self.r_t0d2(tx, d) / d.powi(2), t);
            if print_val {
                println!("[t1d2 == t1d2] t1d2 ={:e}", val_calc);
                println!("[t0d2 -> t1d2] diff ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
            // derivative for temperature+density+density+density
            let val_calc = self.r_t1d3(t, d) / t / d.powi(3);
            let val_diff = romberg_diff(|dx: f64| self.r_t1d2(t, dx) / t / dx.powi(2), d);
            if print_val {
                println!("[t1d3 == t1d3] t1d3 ={:e}", val_calc);
                println!("[t1d2 -> t1d3] diff ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
            let val_diff = romberg_diff(|tx: f64| self.r_t0d3(tx, d) / d.powi(3), t);
            if print_val {
                println!("[t1d3 == t1d3] t1d3 ={:e}", val_calc);
                println!("[t0d3 -> t1d3] didf ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
            // derivative for temperature+temperature
            let val_calc = self.r_t2d0(t, d) / t.powi(2);
            let val_diff = romberg_diff(|tx: f64| self.r_t1d0(tx, d) / tx, t);
            if print_val {
                println!("[t2d0 == t2d0] t2d0 ={:e}", val_calc);
                println!("[t1d0 -> t2d0] diff ={:e}", val_diff);
            } else {
                compare_val(val_calc, val_diff);
            }
        }
    };
}
/// macro_rules! _fn_test
macro_rules! _fn_test {
    ($fluid:expr) => {
        $fluid.c_flash().unwrap();
        $fluid.check_derivatives(false);
        let temp_min = (0.6 * $fluid.T().unwrap()).floor() as u32;
        let temp_max = $fluid.T().unwrap().ceil() as u32;
        let (mut p_s, mut rho_v, mut rho_l): (f64, f64, f64) = (0.0, 0.0, f64::INFINITY);
        for temp in temp_min..temp_max {
            $fluid.t_flash(temp as f64).unwrap();
            if $fluid.p_s().unwrap() < p_s
                || $fluid.rho_v().unwrap() < rho_v
                || $fluid.rho_l().unwrap() > rho_l
            {
                panic!()
            } else {
                p_s = $fluid.p_s().unwrap();
                rho_v = $fluid.rho_v().unwrap();
                rho_l = $fluid.rho_l().unwrap();
            }
        }
    };
}
