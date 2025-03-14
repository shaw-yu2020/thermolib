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
                        (self.temp, self.rho_num, self.is_single_phase) = (temp_c, dens_c, true);
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
                // Vapor phase: eta = 1E-10
                let rhov_num_guess = 1E-10 / self.eta0_coef(temp);
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
                let rhol_num_guess = 0.5 / self.eta0_coef(temp);
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
                let eta0_coef = self.eta0_coef(temp);
                // Iteration from vapor phase: eta = 1E-10
                let rhov_num = self.calc_density(temp, pres, 1E-10 / eta0_coef);
                let lnphi_v = if rhov_num.is_nan() {
                    f64::INFINITY
                } else {
                    self.calc_lnphi(temp, rhov_num)
                };
                // Iteration from liquid phase: eta = 0.5
                let rhol_num = self.calc_density(temp, pres, 0.5 / eta0_coef);
                let lnphi_l = if rhol_num.is_nan() {
                    f64::INFINITY
                } else {
                    self.calc_lnphi(temp, rhol_num)
                };
                // Select the correct output
                if lnphi_v.is_infinite() && lnphi_l.is_infinite() {
                    Err(anyhow!(PcSaftErr::NotConvForTP))
                } else if lnphi_v.is_infinite() {
                    (self.temp, self.rho_num, self.is_single_phase) = (temp, rhol_num, true);
                    Ok(())
                } else if lnphi_l.is_infinite() {
                    (self.temp, self.rho_num, self.is_single_phase) = (temp, rhov_num, true);
                    Ok(())
                } else {
                    if lnphi_v < lnphi_l {
                        (self.temp, self.rho_num, self.is_single_phase) = (temp, rhov_num, true);
                    } else {
                        (self.temp, self.rho_num, self.is_single_phase) = (temp, rhol_num, true);
                    }
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
/// macro_rules! fn_assoc
macro_rules! fn_assoc {
    ($name:ty) => {
        impl $name {
            pub fn t0d0(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                self.XA_flash(temp, rho_num, eta);
                match self.assoc_type {
                    AssocType::Type1 => self.site_t0d0::<1>(self.XA),
                    AssocType::Type2B => 2.0 * self.site_t0d0::<1>(self.XA),
                    AssocType::Type3B => {
                        2.0 * self.site_t0d0::<1>(self.XA)
                            + self.site_t0d0::<2>(2.0 * self.XA - 1.0)
                    }
                    AssocType::Type4C => 4.0 * self.site_t0d0::<1>(self.XA),
                }
            }
            pub fn t0d1(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                self.XA_flash(temp, rho_num, eta);
                match self.assoc_type {
                    AssocType::Type1 => self.site_t0d1::<1>(self.XA, eta),
                    AssocType::Type2B => 2.0 * self.site_t0d1::<1>(self.XA, eta),
                    AssocType::Type3B => {
                        2.0 * self.site_t0d1::<1>(self.XA, eta)
                            + self.site_t0d1::<2>(2.0 * self.XA - 1.0, eta)
                    }
                    AssocType::Type4C => 4.0 * self.site_t0d1::<1>(self.XA, eta),
                }
            }
            pub fn t0d2(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                self.XA_flash(temp, rho_num, eta);
                match self.assoc_type {
                    AssocType::Type1 => self.site_t0d2::<1>(self.XA, eta),
                    AssocType::Type2B => 2.0 * self.site_t0d2::<1>(self.XA, eta),
                    AssocType::Type3B => {
                        2.0 * self.site_t0d2::<1>(self.XA, eta)
                            + self.site_t0d2::<2>(2.0 * self.XA - 1.0, eta)
                    }
                    AssocType::Type4C => 4.0 * self.site_t0d2::<1>(self.XA, eta),
                }
            }
            pub fn t0d3(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                self.XA_flash(temp, rho_num, eta);
                match self.assoc_type {
                    AssocType::Type1 => self.site_t0d3::<1>(self.XA, eta),
                    AssocType::Type2B => 2.0 * self.site_t0d3::<1>(self.XA, eta),
                    AssocType::Type3B => {
                        2.0 * self.site_t0d3::<1>(self.XA, eta)
                            + self.site_t0d3::<2>(2.0 * self.XA - 1.0, eta)
                    }
                    AssocType::Type4C => 4.0 * self.site_t0d3::<1>(self.XA, eta),
                }
            }
            pub fn t0d4(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                self.XA_flash(temp, rho_num, eta);
                match self.assoc_type {
                    AssocType::Type1 => self.site_t0d4::<1>(self.XA, eta),
                    AssocType::Type2B => 2.0 * self.site_t0d4::<1>(self.XA, eta),
                    AssocType::Type3B => {
                        2.0 * self.site_t0d4::<1>(self.XA, eta)
                            + self.site_t0d4::<2>(2.0 * self.XA - 1.0, eta)
                    }
                    AssocType::Type4C => 4.0 * self.site_t0d4::<1>(self.XA, eta),
                }
            }
            pub fn t1d0(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                self.XA_flash(temp, rho_num, eta);
                match self.assoc_type {
                    AssocType::Type1 => self.site_t1d0::<1>(self.XA, eta, eta1),
                    AssocType::Type2B => 2.0 * self.site_t1d0::<1>(self.XA, eta, eta1),
                    AssocType::Type3B => {
                        2.0 * self.site_t1d0::<1>(self.XA, eta, eta1)
                            + self.site_t1d0::<2>(2.0 * self.XA - 1.0, eta, eta1)
                    }
                    AssocType::Type4C => 4.0 * self.site_t1d0::<1>(self.XA, eta, eta1),
                }
            }
            pub fn t1d1(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                self.XA_flash(temp, rho_num, eta);
                match self.assoc_type {
                    AssocType::Type1 => self.site_t1d1::<1>(self.XA, eta, eta1),
                    AssocType::Type2B => 2.0 * self.site_t1d1::<1>(self.XA, eta, eta1),
                    AssocType::Type3B => {
                        2.0 * self.site_t1d1::<1>(self.XA, eta, eta1)
                            + self.site_t1d1::<2>(2.0 * self.XA - 1.0, eta, eta1)
                    }
                    AssocType::Type4C => 4.0 * self.site_t1d1::<1>(self.XA, eta, eta1),
                }
            }
            pub fn t1d2(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                self.XA_flash(temp, rho_num, eta);
                match self.assoc_type {
                    AssocType::Type1 => self.site_t1d2::<1>(self.XA, eta, eta1),
                    AssocType::Type2B => 2.0 * self.site_t1d2::<1>(self.XA, eta, eta1),
                    AssocType::Type3B => {
                        2.0 * self.site_t1d2::<1>(self.XA, eta, eta1)
                            + self.site_t1d2::<2>(2.0 * self.XA - 1.0, eta, eta1)
                    }
                    AssocType::Type4C => 4.0 * self.site_t1d2::<1>(self.XA, eta, eta1),
                }
            }
            pub fn t1d3(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                self.XA_flash(temp, rho_num, eta);
                match self.assoc_type {
                    AssocType::Type1 => self.site_t1d3::<1>(self.XA, eta, eta1),
                    AssocType::Type2B => 2.0 * self.site_t1d3::<1>(self.XA, eta, eta1),
                    AssocType::Type3B => {
                        2.0 * self.site_t1d3::<1>(self.XA, eta, eta1)
                            + self.site_t1d3::<2>(2.0 * self.XA - 1.0, eta, eta1)
                    }
                    AssocType::Type4C => 4.0 * self.site_t1d3::<1>(self.XA, eta, eta1),
                }
            }
            pub fn t2d0(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64, eta2: f64) -> f64 {
                self.XA_flash(temp, rho_num, eta);
                match self.assoc_type {
                    AssocType::Type1 => self.site_t2d0::<1>(self.XA, eta, eta1, eta2),
                    AssocType::Type2B => 2.0 * self.site_t2d0::<1>(self.XA, eta, eta1, eta2),
                    AssocType::Type3B => {
                        2.0 * self.site_t2d0::<1>(self.XA, eta, eta1, eta2)
                            + self.site_t2d0::<2>(2.0 * self.XA - 1.0, eta, eta1, eta2)
                    }
                    AssocType::Type4C => 4.0 * self.site_t2d0::<1>(self.XA, eta, eta1, eta2),
                }
            }
        }
        impl $name {
            fn site_t0d0<const C: i32>(&mut self, x: f64) -> f64 {
                x.ln() - x / 2.0 + 0.5
            }
            fn site_t0d1<const C: i32>(&mut self, x: f64, eta: f64) -> f64 {
                self.site_x1::<C>(x) * self.x_t0d1(eta)
            }
            fn site_t0d2<const C: i32>(&mut self, x: f64, eta: f64) -> f64 {
                self.site_x2::<C>(x) * self.x_t0d1(eta).powi(2)
                    + self.site_x1::<C>(x) * self.x_t0d2(eta)
            }
            fn site_t0d3<const C: i32>(&mut self, x: f64, eta: f64) -> f64 {
                self.site_x3::<C>(x) * self.x_t0d1(eta).powi(3)
                    + 3.0 * self.site_x2::<C>(x) * self.x_t0d1(eta) * self.x_t0d2(eta)
                    + self.site_x1::<C>(x) * self.x_t0d3(eta)
            }
            fn site_t0d4<const C: i32>(&mut self, x: f64, eta: f64) -> f64 {
                self.site_x4::<C>(x) * self.x_t0d1(eta).powi(4)
                    + 6.0 * self.site_x3::<C>(x) * self.x_t0d1(eta).powi(2) * self.x_t0d2(eta)
                    + 3.0 * self.site_x2::<C>(x) * self.x_t0d2(eta).powi(2)
                    + 4.0 * self.site_x2::<C>(x) * self.x_t0d1(eta) * self.x_t0d3(eta)
                    + self.site_x1::<C>(x) * self.x_t0d4(eta)
            }
            fn site_t1d0<const C: i32>(&mut self, x: f64, eta: f64, eta1: f64) -> f64 {
                self.site_x1::<C>(x) * self.x_t1d0(eta, eta1)
            }
            fn site_t1d1<const C: i32>(&mut self, x: f64, eta: f64, eta1: f64) -> f64 {
                self.site_x2::<C>(x) * self.x_t1d0(eta, eta1) * self.x_t0d1(eta)
                    + self.site_x1::<C>(x) * self.x_t1d1(eta, eta1)
            }
            fn site_t1d2<const C: i32>(&mut self, x: f64, eta: f64, eta1: f64) -> f64 {
                self.site_x3::<C>(x) * self.x_t1d0(eta, eta1) * self.x_t0d1(eta).powi(2)
                    + 2.0 * self.site_x2::<C>(x) * self.x_t1d1(eta, eta1) * self.x_t0d1(eta)
                    + self.site_x2::<C>(x) * self.x_t1d0(eta, eta1) * self.x_t0d2(eta)
                    + self.site_x1::<C>(x) * self.x_t1d2(eta, eta1)
            }
            fn site_t1d3<const C: i32>(&mut self, x: f64, eta: f64, eta1: f64) -> f64 {
                self.site_x4::<C>(x) * self.x_t1d0(eta, eta1) * self.x_t0d1(eta).powi(3)
                    + 3.0 * self.site_x3::<C>(x) * self.x_t1d1(eta, eta1) * self.x_t0d1(eta).powi(2)
                    + (3.0 * self.site_x3::<C>(x))
                        * (self.x_t1d0(eta, eta1) * self.x_t0d1(eta) * self.x_t0d2(eta))
                    + 3.0 * self.site_x2::<C>(x) * self.x_t1d2(eta, eta1) * self.x_t0d1(eta)
                    + 3.0 * self.site_x2::<C>(x) * self.x_t1d1(eta, eta1) * self.x_t0d2(eta)
                    + self.site_x2::<C>(x) * self.x_t1d0(eta, eta1) * self.x_t0d3(eta)
                    + self.site_x1::<C>(x) * self.x_t1d3(eta, eta1)
            }
            fn site_t2d0<const C: i32>(&mut self, x: f64, eta: f64, eta1: f64, eta2: f64) -> f64 {
                self.site_x2::<C>(x) * self.x_t1d0(eta, eta1).powi(2)
                    + self.site_x1::<C>(x) * self.x_t2d0(eta, eta1, eta2)
            }
        }
        impl $name {
            fn site_x1<const C: i32>(&self, x: f64) -> f64 {
                (1.0 / x - 0.5) * C as f64
            }
            fn site_x2<const C: i32>(&self, x: f64) -> f64 {
                -1.0 / x.powi(2) * C.pow(2) as f64
            }
            fn site_x3<const C: i32>(&self, x: f64) -> f64 {
                2.0 / x.powi(3) * C.pow(3) as f64
            }
            fn site_x4<const C: i32>(&self, x: f64) -> f64 {
                -6.0 / x.powi(4) * C.pow(4) as f64
            }
        }
        impl $name {
            fn x_t0d1(&mut self, eta: f64) -> f64 {
                if self.x_t0d1.0 != self.XA {
                    self.x_t0d1 = (self.XA, self.xt1() * self.t_t0d1(eta))
                }
                self.x_t0d1.1
            }
            fn x_t0d2(&mut self, eta: f64) -> f64 {
                if self.x_t0d2.0 != self.XA {
                    self.x_t0d2 = (
                        self.XA,
                        self.xt2() * self.t_t0d1(eta).powi(2) + self.xt1() * self.t_t0d2(eta),
                    )
                }
                self.x_t0d2.1
            }
            fn x_t0d3(&mut self, eta: f64) -> f64 {
                if self.x_t0d3.0 != self.XA {
                    self.x_t0d3 = (
                        self.XA,
                        self.xt3() * self.t_t0d1(eta).powi(3)
                            + 3.0 * self.xt2() * self.t_t0d1(eta) * self.t_t0d2(eta)
                            + self.xt1() * self.t_t0d3(eta),
                    )
                }
                self.x_t0d3.1
            }
            fn x_t0d4(&mut self, eta: f64) -> f64 {
                if self.x_t0d4.0 != self.XA {
                    self.x_t0d4 = (
                        self.XA,
                        self.xt4() * self.t_t0d1(eta).powi(4)
                            + 6.0 * self.xt3() * self.t_t0d1(eta).powi(2) * self.t_t0d2(eta)
                            + 3.0 * self.xt2() * self.t_t0d2(eta).powi(2)
                            + 4.0 * self.xt2() * self.t_t0d1(eta) * self.t_t0d3(eta)
                            + self.xt1() * self.t_t0d4(eta),
                    )
                }
                self.x_t0d4.1
            }
            fn x_t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
                if self.x_t1d0.0 != self.XA {
                    self.x_t1d0 = (self.XA, self.xt1() * self.t_t1d0(eta, eta1))
                }
                self.x_t1d0.1
            }
            fn x_t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
                if self.x_t1d1.0 != self.XA {
                    self.x_t1d1 = (
                        self.XA,
                        self.xt2() * self.t_t1d0(eta, eta1) * self.t_t0d1(eta)
                            + self.xt1() * self.t_t1d1(eta, eta1),
                    )
                }
                self.x_t1d1.1
            }
            fn x_t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
                if self.x_t1d2.0 != self.XA {
                    self.x_t1d2 = (
                        self.XA,
                        self.xt3() * self.t_t1d0(eta, eta1) * self.t_t0d1(eta).powi(2)
                            + 2.0 * self.xt2() * self.t_t1d1(eta, eta1) * self.t_t0d1(eta)
                            + self.xt2() * self.t_t1d0(eta, eta1) * self.t_t0d2(eta)
                            + self.xt1() * self.t_t1d2(eta, eta1),
                    )
                }
                self.x_t1d2.1
            }
            fn x_t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
                if self.x_t1d3.0 != self.XA {
                    self.x_t1d3 = (
                        self.XA,
                        self.xt4() * self.t_t1d0(eta, eta1) * self.t_t0d1(eta).powi(3)
                            + 3.0 * self.xt3() * self.t_t1d1(eta, eta1) * self.t_t0d1(eta).powi(2)
                            + 3.0
                                * self.xt3()
                                * self.t_t1d0(eta, eta1)
                                * self.t_t0d1(eta)
                                * self.t_t0d2(eta)
                            + 3.0 * self.xt2() * self.t_t1d2(eta, eta1) * self.t_t0d1(eta)
                            + 3.0 * self.xt2() * self.t_t1d1(eta, eta1) * self.t_t0d2(eta)
                            + self.xt2() * self.t_t1d0(eta, eta1) * self.t_t0d3(eta)
                            + self.xt1() * self.t_t1d3(eta, eta1),
                    )
                }
                self.x_t1d3.1
            }
            fn x_t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
                if self.x_t2d0.0 != self.XA {
                    self.x_t2d0 = (
                        self.XA,
                        self.xt2() * self.t_t1d0(eta, eta1).powi(2)
                            + self.xt1() * self.t_t2d0(eta, eta1, eta2),
                    )
                }
                self.x_t2d0.1
            }
        }
        impl $name {
            fn xt1(&mut self) -> f64 {
                if self.xt1.0 != self.XA {
                    self.xt1 = (
                        self.XA,
                        match self.assoc_type {
                            AssocType::Type1 | AssocType::Type2B => {
                                self.XA.powi(3) / (self.XA - 2.0)
                            }
                            AssocType::Type3B => {
                                (self.XA * (2.0 * self.XA - 1.0)).powi(2)
                                    / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0)
                            }
                            AssocType::Type4C => 2.0 * self.XA.powi(3) / (self.XA - 2.0),
                        },
                    )
                }
                self.xt1.1
            }
            fn xt2(&mut self) -> f64 {
                if self.xt2.0 != self.XA {
                    self.xt2 = (
                        self.XA,
                        2.0 * match self.assoc_type {
                            AssocType::Type1 | AssocType::Type2B => {
                                self.XA.powi(5) / (self.XA - 2.0).powi(3) * (self.XA - 3.0)
                            }
                            AssocType::Type3B => {
                                (self.XA * (2.0 * self.XA - 1.0)).powi(3)
                                    / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0).powi(3)
                                    * (4.0 * self.XA.powi(3) - 12.0 * self.XA.powi(2)
                                        + 6.0 * self.XA
                                        - 1.0)
                            }
                            AssocType::Type4C => {
                                4.0 * self.XA.powi(5) / (self.XA - 2.0).powi(3) * (self.XA - 3.0)
                            }
                        },
                    )
                }
                self.xt2.1
            }
            fn xt3(&mut self) -> f64 {
                if self.xt3.0 != self.XA {
                    self.xt3 = (
                        self.XA,
                        6.0 * match self.assoc_type {
                            AssocType::Type1 | AssocType::Type2B => {
                                self.XA.powi(7) / (self.XA - 2.0).powi(5)
                                    * (self.XA.powi(2) - 6.0 * self.XA + 10.0)
                            }
                            AssocType::Type3B => {
                                (self.XA * (2.0 * self.XA - 1.0)).powi(4)
                                    / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0).powi(5)
                                    * (16.0 * self.XA.powi(6) - 96.0 * self.XA.powi(5)
                                        + (200.0 * self.XA.powi(4) - 160.0 * self.XA.powi(3))
                                        + (62.0 * self.XA.powi(2) - 12.0 * self.XA + 1.0))
                            }
                            AssocType::Type4C => {
                                8.0 * self.XA.powi(7) / (self.XA - 2.0).powi(5)
                                    * (self.XA.powi(2) - 6.0 * self.XA + 10.0)
                            }
                        },
                    )
                }
                self.xt3.1
            }
            fn xt4(&mut self) -> f64 {
                if self.xt4.0 != self.XA {
                    self.xt4 = (
                        self.XA,
                        24.0 * match self.assoc_type {
                            AssocType::Type1 | AssocType::Type2B => {
                                self.XA.powi(9) / (self.XA - 2.0).powi(7)
                                    * (self.XA.powi(3) - 9.0 * self.XA.powi(2) + 29.0 * self.XA
                                        - 35.0)
                            }
                            AssocType::Type3B => {
                                (self.XA * (2.0 * self.XA - 1.0)).powi(5)
                                    / (2.0 * self.XA.powi(2) - 4.0 * self.XA + 1.0).powi(7)
                                    * (64.0 * self.XA.powi(9) - 576.0 * self.XA.powi(8)
                                        + (2080.0 * self.XA.powi(7) - 3808.0 * self.XA.powi(6))
                                        + (3696.0 * self.XA.powi(5) - 2084.0 * self.XA.powi(4))
                                        + (716.0 * self.XA.powi(3) - 150.0 * self.XA.powi(2))
                                        + (18.0 * self.XA - 1.0))
                            }
                            AssocType::Type4C => {
                                16.0 * self.XA.powi(9) / (self.XA - 2.0).powi(7)
                                    * (self.XA.powi(3) - 9.0 * self.XA.powi(2) + 29.0 * self.XA
                                        - 35.0)
                            }
                        },
                    )
                }
                self.xt4.1
            }
        }
    };
}
/// macro_rules! fn_polar
macro_rules! fn_polar {
    ($name:ty) => {
        impl $name {
            pub fn t0d0(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                (self.temp, self.rho_num) = (temp, rho_num);
                self.a2t0d0(eta).powi(2) / (self.a2t0d0(eta) - self.a3t0d0(eta))
            }
            pub fn t0d1(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                (self.temp, self.rho_num) = (temp, rho_num);
                (self.a2t0d1(eta) * self.a2t0d0(eta) * (self.a2t0d0(eta) - 2.0 * self.a3t0d0(eta))
                    + self.a3t0d1(eta) * self.a2t0d0(eta).powi(2))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2)
            }
            pub fn t0d2(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                (self.temp, self.rho_num) = (temp, rho_num);
                (self.a2t0d2(eta)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(2) - 3.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                        + 2.0 * self.a3t0d0(eta).powi(2))
                    + self.a3t0d2(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 2.0 * self.a2t0d1(eta).powi(2) * self.a3t0d0(eta).powi(2)
                    + 2.0 * self.a3t0d1(eta).powi(2) * self.a2t0d0(eta).powi(2)
                    - 4.0
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(3)
            }
            pub fn t0d3(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                (self.temp, self.rho_num) = (temp, rho_num);
                (self.a2t0d3(eta)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(3)
                        - 4.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta)
                        + 5.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(2)
                        - 2.0 * self.a3t0d0(eta).powi(3))
                    + self.a3t0d3(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a2t0d2(eta)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 6.0
                        * self.a2t0d2(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 6.0
                        * self.a3t0d2(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 6.0
                        * self.a3t0d2(eta)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 6.0 * self.a2t0d1(eta).powi(3) * self.a3t0d0(eta).powi(2)
                    + 6.0
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta) + self.a3t0d0(eta))
                    + 6.0 * self.a3t0d1(eta).powi(3) * self.a2t0d0(eta).powi(2)
                    - 6.0
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta) + self.a2t0d0(eta)))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(4)
            }
            pub fn t0d4(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
                (self.temp, self.rho_num) = (temp, rho_num);
                (self.a2t0d4(eta)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(4)
                        - 5.0 * self.a2t0d0(eta).powi(3) * self.a3t0d0(eta)
                        + 9.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta).powi(2)
                        - 7.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(3)
                        + 2.0 * self.a3t0d0(eta).powi(4))
                    + self.a3t0d4(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(3)
                            - 3.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta)
                            + 3.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(2)
                            - self.a3t0d0(eta).powi(3))
                    + 8.0
                        * self.a2t0d3(eta)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 8.0
                        * self.a2t0d3(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 8.0
                        * self.a3t0d3(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 8.0
                        * self.a3t0d3(eta)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a2t0d2(eta).powi(2)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t0d2(eta).powi(2)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 12.0
                        * self.a2t0d2(eta)
                        * self.a3t0d2(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 36.0
                        * self.a2t0d2(eta)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 12.0
                        * self.a2t0d2(eta)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) + self.a2t0d0(eta) * self.a3t0d0(eta)
                            - 2.0 * self.a3t0d0(eta).powi(2))
                    + 24.0
                        * self.a2t0d2(eta)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a3t0d0(eta).powi(2))
                    + 36.0
                        * self.a3t0d2(eta)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 12.0
                        * self.a3t0d2(eta)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta)
                        * (self.a3t0d0(eta).powi(2) + self.a2t0d0(eta) * self.a3t0d0(eta)
                            - 2.0 * self.a2t0d0(eta).powi(2))
                    + 24.0
                        * self.a3t0d2(eta)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a2t0d0(eta).powi(2))
                    + 24.0 * self.a2t0d1(eta).powi(4) * self.a3t0d0(eta).powi(2)
                    + 24.0 * self.a3t0d1(eta).powi(4) * self.a2t0d0(eta).powi(2)
                    + 24.0
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d1(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2)
                            + 4.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 48.0
                        * self.a2t0d1(eta).powi(3)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta))
                    - 48.0
                        * self.a3t0d1(eta).powi(3)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta)))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(5)
            }
            pub fn t1d0(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                (self.temp, self.rho_num) = (temp, rho_num);
                (self.a2t1d0(eta, eta1)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta) - 2.0 * self.a3t0d0(eta))
                    + self.a3t1d0(eta, eta1) * self.a2t0d0(eta).powi(2))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(2)
            }
            pub fn t1d1(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                (self.temp, self.rho_num) = (temp, rho_num);
                (self.a2t1d1(eta, eta1)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(2) - 3.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                        + 2.0 * self.a3t0d0(eta).powi(2))
                    + self.a3t1d1(eta, eta1)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 2.0 * self.a2t1d0(eta, eta1) * self.a2t0d1(eta) * self.a3t0d0(eta).powi(2)
                    - 2.0
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                    + 2.0 * self.a3t1d0(eta, eta1) * self.a3t0d1(eta) * self.a2t0d0(eta).powi(2)
                    - 2.0
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(3)
            }
            pub fn t1d2(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                (self.temp, self.rho_num) = (temp, rho_num);
                (self.a2t1d2(eta, eta1)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(3)
                        - 4.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta)
                        + 5.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(2)
                        - 2.0 * self.a3t0d0(eta).powi(3))
                    + self.a3t1d2(eta, eta1)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 4.0
                        * self.a2t1d1(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 4.0
                        * self.a2t1d1(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 4.0
                        * self.a3t1d1(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 4.0
                        * self.a3t1d1(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 2.0
                        * self.a2t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 2.0
                        * self.a2t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 2.0
                        * self.a3t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 2.0
                        * self.a3t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 6.0
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta).powi(2)
                    - 2.0
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta) + 2.0 * self.a3t0d0(eta))
                    + 4.0
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta) + self.a3t0d0(eta))
                    + 6.0
                        * self.a3t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta).powi(2)
                    + 2.0
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta) + self.a3t0d0(eta))
                    - 4.0
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta) + self.a2t0d0(eta)))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(4)
            }
            pub fn t1d3(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64) -> f64 {
                (self.temp, self.rho_num) = (temp, rho_num);
                (self.a2t1d3(eta, eta1)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(4)
                        - 5.0 * self.a2t0d0(eta).powi(3) * self.a3t0d0(eta)
                        + 9.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta).powi(2)
                        - 7.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(3)
                        + 2.0 * self.a3t0d0(eta).powi(4))
                    + self.a3t1d3(eta, eta1)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(3)
                            - 3.0 * self.a2t0d0(eta).powi(2) * self.a3t0d0(eta)
                            + 3.0 * self.a2t0d0(eta) * self.a3t0d0(eta).powi(2)
                            - self.a3t0d0(eta).powi(3))
                    + 6.0
                        * self.a2t1d2(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t1d2(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 6.0
                        * self.a2t1d2(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 6.0
                        * self.a3t1d2(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 2.0
                        * self.a2t0d3(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 2.0
                        * self.a3t0d3(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 2.0
                        * self.a2t0d3(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 2.0
                        * self.a3t0d3(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a2t1d1(eta, eta1)
                        * self.a2t0d2(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 6.0
                        * self.a2t1d1(eta, eta1)
                        * self.a3t0d2(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t1d1(eta, eta1)
                        * self.a3t0d2(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 6.0
                        * self.a3t1d1(eta, eta1)
                        * self.a2t0d2(eta)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) - 2.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 18.0
                        * self.a2t1d1(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 18.0
                        * self.a3t1d1(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    - 6.0
                        * self.a2t1d1(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta).powi(2) + self.a2t0d0(eta) * self.a3t0d0(eta)
                            - 2.0 * self.a3t0d0(eta).powi(2))
                    - 6.0
                        * self.a3t1d1(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d0(eta)
                        * (self.a3t0d0(eta).powi(2) + self.a2t0d0(eta) * self.a3t0d0(eta)
                            - 2.0 * self.a2t0d0(eta).powi(2))
                    + 12.0
                        * self.a2t1d1(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a3t0d0(eta).powi(2))
                    + 12.0
                        * self.a3t1d1(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a2t0d0(eta).powi(2))
                    - 18.0
                        * self.a2t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 18.0
                        * self.a3t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 6.0
                        * self.a2t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a2t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a3t0d0(eta).powi(2))
                    + 6.0
                        * self.a2t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a2t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t0d2(eta)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a2t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta)
                        * self.a2t0d0(eta)
                        * (2.0 * self.a3t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a2t0d0(eta).powi(2))
                    + 6.0
                        * self.a3t0d2(eta)
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta)
                        * self.a3t0d0(eta)
                        * (2.0 * self.a2t0d0(eta).powi(2)
                            - self.a2t0d0(eta) * self.a3t0d0(eta)
                            - self.a3t0d0(eta).powi(2))
                    + 24.0
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(3)
                        * self.a3t0d0(eta).powi(2)
                    + 24.0
                        * self.a3t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(3)
                        * self.a2t0d0(eta).powi(2)
                    - 36.0
                        * self.a2t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d1(eta)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta))
                    - 36.0
                        * self.a3t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d1(eta)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta))
                    + 12.0
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(2)
                        * self.a2t0d1(eta)
                        * (self.a2t0d0(eta).powi(2)
                            + 4.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    + 12.0
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(2)
                        * self.a3t0d1(eta)
                        * (self.a2t0d0(eta).powi(2)
                            + 4.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                            + self.a3t0d0(eta).powi(2))
                    - 12.0
                        * self.a2t1d0(eta, eta1)
                        * self.a3t0d1(eta).powi(3)
                        * self.a2t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta))
                    - 12.0
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d1(eta).powi(3)
                        * self.a3t0d0(eta)
                        * (self.a2t0d0(eta) + self.a3t0d0(eta)))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(5)
            }
            pub fn t2d0(&mut self, temp: f64, rho_num: f64, eta: f64, eta1: f64, eta2: f64) -> f64 {
                (self.temp, self.rho_num) = (temp, rho_num);
                (self.a2t2d0(eta, eta1, eta2)
                    * self.a2t0d0(eta)
                    * (self.a2t0d0(eta).powi(2) - 3.0 * self.a2t0d0(eta) * self.a3t0d0(eta)
                        + 2.0 * self.a3t0d0(eta).powi(2))
                    + self.a3t2d0(eta, eta1, eta2)
                        * self.a2t0d0(eta).powi(2)
                        * (self.a2t0d0(eta) - self.a3t0d0(eta))
                    + 2.0 * self.a2t1d0(eta, eta1).powi(2) * self.a3t0d0(eta).powi(2)
                    + 2.0 * self.a3t1d0(eta, eta1).powi(2) * self.a2t0d0(eta).powi(2)
                    - 4.0
                        * self.a2t1d0(eta, eta1)
                        * self.a3t1d0(eta, eta1)
                        * self.a2t0d0(eta)
                        * self.a3t0d0(eta))
                    / (self.a2t0d0(eta) - self.a3t0d0(eta)).powi(3)
            }
        }
        impl $name {
            fn a2t0d0(&mut self, eta: f64) -> f64 {
                self.a2_coef * self.rho_num / self.temp.powi(2)
                    * (self.j2a.t0d0(eta) + self.epsilon / self.temp * self.j2b.t0d0(eta))
            }
            fn a2t0d1(&mut self, eta: f64) -> f64 {
                self.a2_coef * self.rho_num / self.temp.powi(2)
                    * (self.j2a.t0d1(eta) + self.epsilon / self.temp * self.j2b.t0d1(eta))
            }
            fn a2t0d2(&mut self, eta: f64) -> f64 {
                self.a2_coef * self.rho_num / self.temp.powi(2)
                    * (self.j2a.t0d2(eta) + self.epsilon / self.temp * self.j2b.t0d2(eta))
            }
            fn a2t0d3(&mut self, eta: f64) -> f64 {
                self.a2_coef * self.rho_num / self.temp.powi(2)
                    * (self.j2a.t0d3(eta) + self.epsilon / self.temp * self.j2b.t0d3(eta))
            }
            fn a2t0d4(&mut self, eta: f64) -> f64 {
                self.a2_coef * self.rho_num / self.temp.powi(2)
                    * (self.j2a.t0d4(eta) + self.epsilon / self.temp * self.j2b.t0d4(eta))
            }
            fn a2t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
                self.a2_coef * self.rho_num / self.temp.powi(2)
                    * (self.j2a.t1d0(eta, eta1)
                        + self.epsilon / self.temp * self.j2b.t1d0(eta, eta1))
            }
            fn a2t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
                self.a2_coef * self.rho_num / self.temp.powi(2)
                    * (self.j2a.t1d1(eta, eta1)
                        + self.epsilon / self.temp * self.j2b.t1d1(eta, eta1))
            }
            fn a2t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
                self.a2_coef * self.rho_num / self.temp.powi(2)
                    * (self.j2a.t1d2(eta, eta1)
                        + self.epsilon / self.temp * self.j2b.t1d2(eta, eta1))
            }
            fn a2t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
                self.a2_coef * self.rho_num / self.temp.powi(2)
                    * (self.j2a.t1d3(eta, eta1)
                        + self.epsilon / self.temp * self.j2b.t1d3(eta, eta1))
            }
            fn a2t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
                self.a2_coef * self.rho_num / self.temp.powi(2)
                    * (self.j2a.t2d0(eta, eta1, eta2)
                        + self.epsilon / self.temp * self.j2b.t2d0(eta, eta1, eta2))
            }
        }
        impl $name {
            fn a3t0d0(&mut self, eta: f64) -> f64 {
                self.a3_coef * self.rho_num.powi(2) / self.temp.powi(3) * self.j3c.t0d0(eta)
            }
            fn a3t0d1(&mut self, eta: f64) -> f64 {
                self.a3_coef * self.rho_num.powi(2) / self.temp.powi(3) * self.j3c.t0d1(eta)
            }
            fn a3t0d2(&mut self, eta: f64) -> f64 {
                self.a3_coef * self.rho_num.powi(2) / self.temp.powi(3) * self.j3c.t0d2(eta)
            }
            fn a3t0d3(&mut self, eta: f64) -> f64 {
                self.a3_coef * self.rho_num.powi(2) / self.temp.powi(3) * self.j3c.t0d3(eta)
            }
            fn a3t0d4(&mut self, eta: f64) -> f64 {
                self.a3_coef * self.rho_num.powi(2) / self.temp.powi(3) * self.j3c.t0d4(eta)
            }
            fn a3t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
                self.a3_coef * self.rho_num.powi(2) / self.temp.powi(3) * self.j3c.t1d0(eta, eta1)
            }
            fn a3t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
                self.a3_coef * self.rho_num.powi(2) / self.temp.powi(3) * self.j3c.t1d1(eta, eta1)
            }
            fn a3t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
                self.a3_coef * self.rho_num.powi(2) / self.temp.powi(3) * self.j3c.t1d2(eta, eta1)
            }
            fn a3t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
                self.a3_coef * self.rho_num.powi(2) / self.temp.powi(3) * self.j3c.t1d3(eta, eta1)
            }
            fn a3t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
                self.a3_coef * self.rho_num.powi(2) / self.temp.powi(3)
                    * self.j3c.t2d0(eta, eta1, eta2)
            }
        }
        /// J2aTerm
        struct J2aTerm {
            a0: f64,
            a1: f64,
            a2: f64,
            a3: f64,
            a4: f64,
            // cached variables
            t0d0: (f64, f64),
            t0d1: (f64, f64),
            t0d2: (f64, f64),
            t0d3: (f64, f64),
        }
        impl J2aTerm {
            fn new(m: f64) -> Self {
                let m1 = (m - 1.0) / m; // (m-1)/m
                let m12 = (m - 1.0) * (m - 2.0) / m.powi(2); // (m-1)/m * (m-2)/m
                Self {
                    a0: A00 + m1 * A10 + m12 * A20,
                    a1: A01 + m1 * A11 + m12 * A21,
                    a2: A02 + m1 * A12 + m12 * A22,
                    a3: A03 + m1 * A13 + m12 * A23,
                    a4: A04 + m1 * A14 + m12 * A24,
                    // cached variables
                    t0d0: (0.0, 0.0),
                    t0d1: (0.0, 0.0),
                    t0d2: (0.0, 0.0),
                    t0d3: (0.0, 0.0),
                }
            }
            /// equal to = [rho/T^2 *J2]_t0d0 / {rho/T^2}
            /// equal to = J2t0d0
            fn t0d0(&mut self, eta: f64) -> f64 {
                if eta != self.t0d0.0 {
                    self.t0d0 = (
                        eta,
                        self.a0
                            + self.a1 * eta
                            + self.a2 * eta.powi(2)
                            + self.a3 * eta.powi(3)
                            + self.a4 * eta.powi(4),
                    )
                }
                self.t0d0.1
            }
            /// equal to = [rho/T^2 *J2]_t0d1 / {rho/T^2}
            /// equal to = J2t0d0 + rho * J2t0d1
            fn t0d1(&mut self, eta: f64) -> f64 {
                // 1 + n
                if eta != self.t0d1.0 {
                    self.t0d1 = (
                        eta,
                        self.a0
                            + self.a1 * 2.0 * eta
                            + self.a2 * 3.0 * eta.powi(2)
                            + self.a3 * 4.0 * eta.powi(3)
                            + self.a4 * 5.0 * eta.powi(4),
                    )
                }
                self.t0d1.1
            }
            /// equal to = [rho/T^2 *J2]_t0d2 / {rho/T^2}
            /// equal to = 2 * rho * J2t0d1 + rho^2 * J2t0d2
            fn t0d2(&mut self, eta: f64) -> f64 {
                // 2 * n + n * ( n - 1 ) = n * ( n + 1 )
                if eta != self.t0d2.0 {
                    self.t0d2 = (
                        eta,
                        self.a1 * 2.0 * eta
                            + self.a2 * 6.0 * eta.powi(2)
                            + self.a3 * 12.0 * eta.powi(3)
                            + self.a4 * 20.0 * eta.powi(4),
                    )
                }
                self.t0d2.1
            }
            /// equal to = [rho/T^2 *J2]_t0d3 / {rho/T^2}
            /// equal to = 3 * rho^2 * J2t0d2 + rho^3 * J2t0d3
            fn t0d3(&mut self, eta: f64) -> f64 {
                // 3 * n * ( n - 1 ) + n * ( n - 1 ) * ( n - 2 ) = n * ( n - 1 ) * ( n + 1 )
                if eta != self.t0d3.0 {
                    self.t0d3 = (
                        eta,
                        self.a2 * 6.0 * eta.powi(2)
                            + self.a3 * 24.0 * eta.powi(3)
                            + self.a4 * 60.0 * eta.powi(4),
                    )
                }
                self.t0d3.1
            }
            /// equal to = [rho/T^2 *J2]_t0d4 / {rho/T^2}
            /// equal to = 4 * rho^3 * J3t0d3 + rho^4 * J2t0d4
            fn t0d4(&mut self, eta: f64) -> f64 {
                // 4 * n * ( n - 1 ) * ( n - 2 ) + n * ( n - 1 ) * ( n - 2 ) * ( n - 3 )
                // = n * ( n - 1 ) * ( n - 2 ) * ( n + 1 )
                self.a3 * 24.0 * eta.powi(3) + self.a4 * 120.0 * eta.powi(4)
            }
            /// equal to = [rho/T^2 *J2]_t1d0 / {rho/T^2}
            /// equal to = -2 * J2t0d0 + T * J2t1d0
            fn t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
                -2.0 * self.t0d0(eta)
                    + eta1
                        * (self.a1
                            + self.a2 * 2.0 * eta
                            + self.a3 * 3.0 * eta.powi(2)
                            + self.a4 * 4.0 * eta.powi(3))
            }
            /// equal to = [rho/T^2 *J2]_t1d1 / {rho/T^2}
            /// equal to = -2 * ( J2t0d0 + rho * J2t0d1 )
            ///            +T * ( J2t1d0 + rho * J2t1d1 )
            fn t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
                -2.0 * self.t0d1(eta)
                    + eta1
                        * (self.a1 * 2.0
                            + self.a2 * 6.0 * eta
                            + self.a3 * 12.0 * eta.powi(2)
                            + self.a4 * 20.0 * eta.powi(3))
            }
            /// equal to = [rho/T^2 *J2]_t1d2 / {rho/T^2}
            /// equal to = -2 * ( 2 * rho * J2t0d1 + rho^2 * J2t0d2 )
            ///            +T * ( 2 * rho * J2t1d1 + rho^2 * J2t1d2 )
            fn t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
                -2.0 * self.t0d2(eta)
                    + eta1
                        * (self.a1 * 2.0
                            + self.a2 * 12.0 * eta
                            + self.a3 * 36.0 * eta.powi(2)
                            + self.a4 * 80.0 * eta.powi(3))
            }
            /// equal to = [rho/T^2 *J2]_t1d3 / {rho/T^2}
            /// equal to = -2 * ( 6 * rho * J2t0d1 + 6 * rho^2 * J2t0d1 + rho^3 * J2t0d3 )
            ///            +T * ( 6 * rho * J2t1d1 + 6 * rho^2 * J2t1d2 + rho^3 * J2t1d3 )
            fn t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
                -2.0 * self.t0d3(eta)
                    + eta1
                        * (self.a2 * 12.0 * eta
                            + self.a3 * 72.0 * eta.powi(2)
                            + self.a4 * 240.0 * eta.powi(3))
            }
            /// equal to = [rho/T^2 *J2]_t2d0 / {rho/T^2}
            /// equal to = 6 * J2t0d0 - 4 * T * J2t1d0 + T^2 * J2t2d0
            fn t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
                6.0 * self.t0d0(eta)
                    + (eta2 - 4.0 * eta1)
                        * (self.a1
                            + self.a2 * 2.0 * eta
                            + self.a3 * 3.0 * eta.powi(2)
                            + self.a4 * 4.0 * eta.powi(3))
                    + eta1.powi(2)
                        * (self.a2 * 2.0 + self.a3 * 6.0 * eta + self.a4 * 12.0 * eta.powi(2))
            }
        }
        /// J2bTerm
        struct J2bTerm {
            b0: f64,
            b1: f64,
            b2: f64,
            // cached variables
            t0d0: (f64, f64),
            t0d1: (f64, f64),
            t0d2: (f64, f64),
            t0d3: (f64, f64),
        }
        impl J2bTerm {
            fn new(m: f64) -> Self {
                let m1 = (m - 1.0) / m; // (m-1)/m
                let m12 = (m - 1.0) * (m - 2.0) / m.powi(2); // (m-1)/m * (m-2)/m
                Self {
                    b0: B00 + m1 * B10 + m12 * B20,
                    b1: B01 + m1 * B11 + m12 * B21,
                    b2: B02 + m1 * B12 + m12 * B22,
                    // cached variables
                    t0d0: (0.0, 0.0),
                    t0d1: (0.0, 0.0),
                    t0d2: (0.0, 0.0),
                    t0d3: (0.0, 0.0),
                }
            }
            /// equal to = [rho/T^3 *J2]_t0d0 / {rho/T^3}
            /// equal to = J2t0d0
            fn t0d0(&mut self, eta: f64) -> f64 {
                if eta != self.t0d0.0 {
                    self.t0d0 = (eta, self.b0 + self.b1 * eta + self.b2 * eta.powi(2))
                }
                self.t0d0.1
            }
            /// equal to = [rho/T^3 *J2]_t0d1 / {rho/T^3}
            /// equal to = J2t0d0 + rho * J2t0d1
            fn t0d1(&mut self, eta: f64) -> f64 {
                // 1 + n
                if eta != self.t0d1.0 {
                    self.t0d1 = (
                        eta,
                        self.b0 + self.b1 * 2.0 * eta + self.b2 * 3.0 * eta.powi(2),
                    )
                }
                self.t0d1.1
            }
            /// equal to = [rho/T^3 *J2]_t0d2 / {rho/T^3}
            /// equal to = 2 * rho * J2t0d1 + rho^2 * J2t0d2
            fn t0d2(&mut self, eta: f64) -> f64 {
                // 2 * n + n * ( n - 1 ) = n * ( n + 1 )
                if eta != self.t0d2.0 {
                    self.t0d2 = (eta, self.b1 * 2.0 * eta + self.b2 * 6.0 * eta.powi(2))
                }
                self.t0d2.1
            }
            /// equal to = [rho/T^3 *J2]_t0d3 / {rho/T^3}
            /// equal to = 3 * rho^2 * J2t0d2 + rho^3 * J2t0d3
            fn t0d3(&mut self, eta: f64) -> f64 {
                // 3 * n * ( n - 1 ) + n * ( n - 1 ) * ( n - 2 ) = n * ( n - 1 ) * ( n + 1 )
                if eta != self.t0d3.0 {
                    self.t0d3 = (eta, self.b2 * 6.0 * eta.powi(2))
                }
                self.t0d3.1
            }
            /// equal to = [rho/T^3 *J2]_t0d4 / {rho/T^3}
            /// equal to = 4 * rho^3 * J3t0d3 + rho^4 * J2t0d4
            #[inline]
            fn t0d4(&mut self, _eta: f64) -> f64 {
                // 4 * n * ( n - 1 ) * ( n - 2 ) + n * ( n - 1 ) * ( n - 2 ) * ( n - 3 )
                0.0
            }
            /// equal to = [rho/T^3 *J2]_t1d0 / {rho/T^3}
            /// equal to = -3 * J2t0d0 + T * J2t1d0
            fn t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
                -3.0 * self.t0d0(eta) + eta1 * (self.b1 + self.b2 * 2.0 * eta)
            }
            /// equal to = [rho/T^3 *J2]_t1d1 / {rho/T^3}
            /// equal to = -3 * ( J2t0d0 + rho * J2t0d1 )
            ///            +T * ( J2t1d0 + rho * J2t1d1 )
            fn t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
                -3.0 * self.t0d1(eta) + eta1 * (self.b1 * 2.0 + self.b2 * 6.0 * eta)
            }
            /// equal to = [rho/T^3 *J2]_t1d2 / {rho/T^3}
            /// equal to = -3 * ( 2 * rho * J2t0d1 + rho^2 * J2t0d2 )
            ///            +T * ( 2 * rho * J2t1d1 + rho^2 * J2t1d2 )
            fn t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
                -3.0 * self.t0d2(eta) + eta1 * (self.b1 * 2.0 + self.b2 * 12.0 * eta)
            }
            /// equal to = [rho/T^3 *J2]_t1d3 / {rho/T^3}
            /// equal to = -3 * ( 6 * rho * J2t0d1 + 6 * rho^2 * J2t0d1 + rho^3 * J2t0d3 )
            ///            +T * ( 6 * rho * J2t1d1 + 6 * rho^2 * J2t1d2 + rho^3 * J2t1d3 )
            fn t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
                -3.0 * self.t0d3(eta) + eta1 * self.b2 * 12.0 * eta
            }
            /// equal to = [rho/T^3 *J2]_t2d0 / {rho/T^3}
            /// equal to = 12 * J2t0d0 - 6 * T * J2t1d0 + T^2 * J2t2d0
            fn t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
                12.0 * self.t0d0(eta)
                    + (eta2 - 6.0 * eta1) * (self.b1 + self.b2 * 2.0 * eta)
                    + eta1.powi(2) * self.b2 * 2.0
            }
        }
        /// J3cTerm
        struct J3cTerm {
            c0: f64,
            c1: f64,
            c2: f64,
            c3: f64,
            // cached variables
            t0d0: (f64, f64),
            t0d1: (f64, f64),
            t0d2: (f64, f64),
            t0d3: (f64, f64),
        }
        impl J3cTerm {
            fn new(m: f64) -> Self {
                let m1 = (m - 1.0) / m; // (m-1)/m
                let m12 = (m - 1.0) * (m - 2.0) / m.powi(2); // (m-1)/m * (m-2)/m
                Self {
                    c0: C00 + m1 * C10 + m12 * C20,
                    c1: C01 + m1 * C11 + m12 * C21,
                    c2: C02 + m1 * C12 + m12 * C22,
                    c3: C03 + m1 * C13 + m12 * C23,
                    // cached variables
                    t0d0: (0.0, 0.0),
                    t0d1: (0.0, 0.0),
                    t0d2: (0.0, 0.0),
                    t0d3: (0.0, 0.0),
                }
            }
            /// equal to = [rho^2/T^3 *J3]_t0d0 / {rho^2/T^3}
            /// equal to = J3t0d0
            fn t0d0(&mut self, eta: f64) -> f64 {
                if eta != self.t0d0.0 {
                    self.t0d0 = (
                        eta,
                        self.c0 + self.c1 * eta + self.c2 * eta.powi(2) + self.c3 * eta.powi(3),
                    )
                }
                self.t0d0.1
            }
            /// equal to = [rho^2/T^3 *J3]_t0d1 / {rho^2/T^3}
            /// equal to = 2 * J3t0d0 + rho * J3t0d1
            fn t0d1(&mut self, eta: f64) -> f64 {
                // 2 + n
                if eta != self.t0d1.0 {
                    self.t0d1 = (
                        eta,
                        2.0 * self.c0
                            + self.c1 * 3.0 * eta
                            + self.c2 * 4.0 * eta.powi(2)
                            + self.c3 * 5.0 * eta.powi(3),
                    )
                }
                self.t0d1.1
            }
            /// equal to = [rho^2/T^3 *J3]_t0d2 / {rho^2/T^3}
            /// equal to = 2 * J3t0d0 + 4 * rho * J3t0d1 + rho^2 * J3t0d2
            fn t0d2(&mut self, eta: f64) -> f64 {
                // 2 + 4 * n + n * ( n - 1 )
                // = 2 + 3 * n + n^2
                if eta != self.t0d2.0 {
                    self.t0d2 = (
                        eta,
                        2.0 * self.c0
                            + self.c1 * 6.0 * eta
                            + self.c2 * 12.0 * eta.powi(2)
                            + self.c3 * 20.0 * eta.powi(3),
                    )
                }
                self.t0d2.1
            }
            /// equal to = [rho^2/T^3 *J3]_t0d3 / {rho^2/T^3}
            /// equal to = rho * ( 6 * J3t0d1 + 6 * rho * J3t0d2 + rho^2 * J3t0d3)
            fn t0d3(&mut self, eta: f64) -> f64 {
                // 6 * n + 6 * n * ( n - 1 ) + n * ( n - 1 ) * ( n - 2 )
                // = n * ( n^2 + 3 * n + 2 )
                if eta != self.t0d3.0 {
                    self.t0d3 = (
                        eta,
                        self.c1 * 6.0 * eta
                            + self.c2 * 24.0 * eta.powi(2)
                            + self.c3 * 60.0 * eta.powi(3),
                    )
                }
                self.t0d3.1
            }
            /// equal to = [rho^2/T^3 *J3]_t0d4 / {rho^2/T^3}
            /// equal to = rho^2 * ( 12 * J3t0d2 + 8 * rho * J3t0d3 + rho^2 * J3t0d4)
            fn t0d4(&mut self, eta: f64) -> f64 {
                // 12 * n * ( n - 1 ) + 8 * n * ( n - 1 ) * ( n - 2 )
                // + n * ( n - 1 ) * ( n - 2 ) * ( n - 3 )
                // = n * ( n - 1 ) * ( n^2 + 3 * n + 2 )
                self.c2 * 24.0 * eta.powi(2) + self.c3 * 120.0 * eta.powi(3)
            }
            /// equal to = [rho^2/T^3 *J3]_t1d0 / {rho^2/T^3}
            /// equal to = -3 * J3t0d0 + T * J3t1d0
            fn t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
                -3.0 * self.t0d0(eta)
                    + eta1 * (self.c1 + self.c2 * 2.0 * eta + self.c3 * 3.0 * eta.powi(2))
            }
            /// equal to = [rho^2/T^3 *J3]_t1d1 / {rho^2/T^3}
            /// equal to = -3 * ( 2 * J3t0d0 + rho * J3t0d1 )
            ///            +T * ( 2 * J3t1d0 + rho * J3t1d1 )
            fn t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
                -3.0 * self.t0d1(eta)
                    + eta1 * (self.c1 * 3.0 + self.c2 * 8.0 * eta + self.c3 * 15.0 * eta.powi(2))
            }
            /// equal to = [rho^2/T^3 *J3]_t1d2 / {rho^2/T^3}
            /// equal to = -3 * ( 2 * J3t0d0 + 4 * rho * J3t0d1 + rho^2 * J3t0d2 )
            ///            +T * ( 2 * J3t1d0 + 4 * rho * J3t1d1 + rho^2 * J3t1d2 )
            fn t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
                -3.0 * self.t0d2(eta)
                    + eta1 * (self.c1 * 6.0 + self.c2 * 24.0 * eta + self.c3 * 60.0 * eta.powi(2))
            }
            /// equal to = [rho^2/T^3 *J3]_t1d3 / {rho^2/T^3}
            /// equal to = -3 * ( 6 * rho * J3t0d1 + 6 * rho^2 * J3t0d2 + rho^3 * J3t0d3 )
            ///            +T * ( 6 * rho * J3t1d1 + 6 * rho^2 * J3t1d2 + rho^3 * J3t1d3 )
            fn t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
                -3.0 * self.t0d3(eta)
                    + eta1 * (self.c1 * 6.0 + self.c2 * 48.0 * eta + self.c3 * 180.0 * eta.powi(2))
            }
            /// equal to = [rho^2/T^3 *J3]_t2d0 / {rho^2/T^3}
            /// equal to = 12 * J3t0d0 - 6 * T * J3t1d0 + T^2 * J3t2d0
            fn t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
                12.0 * self.t0d0(eta)
                    + (eta2 - 6.0 * eta1)
                        * (self.c1 + self.c2 * 2.0 * eta + self.c3 * 3.0 * eta.powi(2))
                    + eta1.powi(2) * (self.c2 * 2.0 + self.c3 * 6.0 * eta)
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
