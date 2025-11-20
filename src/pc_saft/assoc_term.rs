use super::GiiTerm;
#[derive(Clone)]
pub enum AssocType {
    Type1,
    Type2B,
    Type3B,
    Type4C,
}
/// AssocTerm
#[allow(non_snake_case)]
pub struct AssocTerm {
    gii: GiiTerm,
    assoc_type: AssocType,
    kappa_AB_sigma3_dens: f64,
    epsilon_AB_temp: f64,
    kappa_AB_sigma3: f64,
    epsilon_AB: f64,
    dens: f64,
    temp: f64,
    XA: f64,
    x: f64,
    // cached variables
    xt1: (f64, f64),
    xt2: (f64, f64),
    xt3: (f64, f64),
    x_t0d1: (f64, f64),
    x_t0d2: (f64, f64),
    x_t0d3: (f64, f64),
    x_t1d0: (f64, f64),
    x_t1d1: (f64, f64),
    x_t2d0: (f64, f64),
    t_t0d0: (f64, f64),
    t_t0d1: (f64, f64),
    t_t0d2: (f64, f64),
    t_t0d3: (f64, f64),
    t_t1d0: (f64, f64),
    t_t1d1: (f64, f64),
    t_t2d0: (f64, f64),
}
#[allow(non_snake_case)]
impl AssocTerm {
    pub fn new_1_term(x: f64, kappa_AB_sigma3: f64, epsilon_AB: f64) -> Self {
        Self {
            gii: GiiTerm::new(&[0.0]),
            assoc_type: AssocType::Type1,
            kappa_AB_sigma3_dens: 0.0,
            epsilon_AB_temp: 0.0,
            kappa_AB_sigma3,
            epsilon_AB,
            dens: 0.0,
            temp: 0.0,
            XA: 1.0,
            x,
            // cached variables
            xt1: (0.0, 0.0),
            xt2: (0.0, 0.0),
            xt3: (0.0, 0.0),
            x_t0d1: (0.0, 0.0),
            x_t0d2: (0.0, 0.0),
            x_t0d3: (0.0, 0.0),
            x_t1d0: (0.0, 0.0),
            x_t1d1: (0.0, 0.0),
            x_t2d0: (0.0, 0.0),
            t_t0d0: (0.0, 0.0),
            t_t0d1: (0.0, 0.0),
            t_t0d2: (0.0, 0.0),
            t_t0d3: (0.0, 0.0),
            t_t1d0: (0.0, 0.0),
            t_t1d1: (0.0, 0.0),
            t_t2d0: (0.0, 0.0),
        }
    }
    pub fn new_x_frac(&self, x: f64) -> Self {
        Self {
            x,
            gii: self.gii.clone(),
            assoc_type: self.assoc_type.clone(),
            ..*self
        }
    }
    pub fn new_2B_term(x: f64, kappa_AB_sigma3: f64, epsilon_AB: f64) -> Self {
        Self {
            assoc_type: AssocType::Type2B,
            ..AssocTerm::new_1_term(x, kappa_AB_sigma3, epsilon_AB)
        }
    }
    pub fn new_3B_term(x: f64, kappa_AB_sigma3: f64, epsilon_AB: f64) -> Self {
        Self {
            assoc_type: AssocType::Type3B,
            ..AssocTerm::new_1_term(x, kappa_AB_sigma3, epsilon_AB)
        }
    }
    pub fn new_4C_term(x: f64, kappa_AB_sigma3: f64, epsilon_AB: f64) -> Self {
        Self {
            assoc_type: AssocType::Type4C,
            ..AssocTerm::new_1_term(x, kappa_AB_sigma3, epsilon_AB)
        }
    }
}
/// AssocGlyTerm
#[allow(non_snake_case)]
pub struct AssocGlyTerm {
    assoc_type: AssocType,
    kappa_AB_sigma3_dens: f64,
    epsilon_AB_temp: f64,
    kappa_AB_sigma3: f64,
    epsilon_AB: f64,
    dens: f64,
    temp: f64,
    XA: f64,
    x: f64,
    // parameters for gly
    c0: f64,
    c1: f64,
    c2: f64,
    // cached variables
    xt1: (f64, f64),
    xt2: (f64, f64),
    xt3: (f64, f64),
    x_t0d1: (f64, f64),
    x_t0d2: (f64, f64),
    x_t0d3: (f64, f64),
    x_t1d0: (f64, f64),
    x_t1d1: (f64, f64),
    x_t2d0: (f64, f64),
    t_t0d0: (f64, f64),
    t_t0d1: (f64, f64),
    t_t0d2: (f64, f64),
    t_t0d3: (f64, f64),
    t_t1d0: (f64, f64),
    t_t1d1: (f64, f64),
    t_t2d0: (f64, f64),
    g_t0d0: (f64, f64),
    g_t0d1: (f64, f64),
    g_t0d2: (f64, f64),
    g_t0d3: (f64, f64),
    g_t1d0: (f64, f64),
    g_t1d1: (f64, f64),
    g_t2d0: (f64, f64),
}
#[allow(non_snake_case)]
impl AssocGlyTerm {
    pub fn new_1_term(
        x: f64,
        kappa_AB_sigma3: f64,
        epsilon_AB: f64,
        c0: f64,
        c1: f64,
        c2: f64,
    ) -> Self {
        Self {
            assoc_type: AssocType::Type1,
            kappa_AB_sigma3_dens: 0.0,
            epsilon_AB_temp: 0.0,
            kappa_AB_sigma3,
            epsilon_AB,
            dens: 0.0,
            temp: 0.0,
            XA: 1.0,
            x,
            // parameters for gly
            c0,
            c1,
            c2,
            // cached variables
            xt1: (0.0, 0.0),
            xt2: (0.0, 0.0),
            xt3: (0.0, 0.0),
            x_t0d1: (0.0, 0.0),
            x_t0d2: (0.0, 0.0),
            x_t0d3: (0.0, 0.0),
            x_t1d0: (0.0, 0.0),
            x_t1d1: (0.0, 0.0),
            x_t2d0: (0.0, 0.0),
            t_t0d0: (0.0, 0.0),
            t_t0d1: (0.0, 0.0),
            t_t0d2: (0.0, 0.0),
            t_t0d3: (0.0, 0.0),
            t_t1d0: (0.0, 0.0),
            t_t1d1: (0.0, 0.0),
            t_t2d0: (0.0, 0.0),
            g_t0d0: (0.0, 0.0),
            g_t0d1: (0.0, 0.0),
            g_t0d2: (0.0, 0.0),
            g_t0d3: (0.0, 0.0),
            g_t1d0: (0.0, 0.0),
            g_t1d1: (0.0, 0.0),
            g_t2d0: (0.0, 0.0),
        }
    }
    pub fn new_x_frac(&self, x: f64) -> Self {
        Self {
            x,
            assoc_type: self.assoc_type.clone(),
            ..*self
        }
    }
    pub fn new_2B_term(
        x: f64,
        kappa_AB_sigma3: f64,
        epsilon_AB: f64,
        c0: f64,
        c1: f64,
        c2: f64,
    ) -> Self {
        Self {
            assoc_type: AssocType::Type2B,
            ..AssocGlyTerm::new_1_term(x, kappa_AB_sigma3, epsilon_AB, c0, c1, c2)
        }
    }
    pub fn new_3B_term(
        x: f64,
        kappa_AB_sigma3: f64,
        epsilon_AB: f64,
        c0: f64,
        c1: f64,
        c2: f64,
    ) -> Self {
        Self {
            assoc_type: AssocType::Type3B,
            ..AssocGlyTerm::new_1_term(x, kappa_AB_sigma3, epsilon_AB, c0, c1, c2)
        }
    }
    pub fn new_4C_term(
        x: f64,
        kappa_AB_sigma3: f64,
        epsilon_AB: f64,
        c0: f64,
        c1: f64,
        c2: f64,
    ) -> Self {
        Self {
            assoc_type: AssocType::Type4C,
            ..AssocGlyTerm::new_1_term(x, kappa_AB_sigma3, epsilon_AB, c0, c1, c2)
        }
    }
}
/// macro_rules! fn_assoc
macro_rules! fn_assoc {
    ($name:ty) => {
        impl $name {
            #[allow(non_snake_case)]
            fn XA_flash(&mut self, temp: f64, rho_num: f64, zeta2t0: f64, zeta3t0: f64, dit0: f64) {
                if temp != self.temp || rho_num != self.dens {
                    (self.temp, self.epsilon_AB_temp) = (temp, self.epsilon_AB / temp);
                    (self.dens, self.kappa_AB_sigma3_dens) =
                        (rho_num, self.kappa_AB_sigma3 * rho_num);
                    self.XA = match self.assoc_type {
                        AssocType::Type1 | AssocType::Type2B => {
                            (-1.0 + (1.0 + 4.0 * self.t_t0d0(zeta2t0, zeta3t0, dit0)).sqrt())
                                / (2.0 * self.t_t0d0(zeta2t0, zeta3t0, dit0))
                        }
                        AssocType::Type3B => {
                            (-(1.0 - self.t_t0d0(zeta2t0, zeta3t0, dit0))
                                + ((1.0 - self.t_t0d0(zeta2t0, zeta3t0, dit0)).powi(2)
                                    + 8.0 * self.t_t0d0(zeta2t0, zeta3t0, dit0))
                                .sqrt())
                                / (4.0 * self.t_t0d0(zeta2t0, zeta3t0, dit0))
                        }
                        AssocType::Type4C => {
                            (-1.0 + (1.0 + 8.0 * self.t_t0d0(zeta2t0, zeta3t0, dit0)).sqrt())
                                / (4.0 * self.t_t0d0(zeta2t0, zeta3t0, dit0))
                        }
                    };
                }
            }
            #[allow(clippy::too_many_arguments)]
            pub fn mu_k<'a>(
                &mut self,
                temp: f64,
                rho_num: f64,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
                k: usize,
                zeta2_k: &'a [f64],
                zeta3_k: &'a [f64],
            ) -> impl Iterator<Item = f64> + use<'a> {
                let val = self.t0d0(temp, rho_num, zeta2t0, zeta3t0, dit0) / self.x
                    + rho_num
                        * self.x
                        * match self.assoc_type {
                            AssocType::Type1 => self.site_mu_k::<1>(
                                self.XA, zeta2t0, zeta3t0, dit0, zeta2_k[k], zeta3_k[k],
                            ),
                            AssocType::Type2B => {
                                2.0 * self.site_mu_k::<1>(
                                    self.XA, zeta2t0, zeta3t0, dit0, zeta2_k[k], zeta3_k[k],
                                )
                            }
                            AssocType::Type3B => {
                                2.0 * self.site_mu_k::<1>(
                                    self.XA, zeta2t0, zeta3t0, dit0, zeta2_k[k], zeta3_k[k],
                                ) + self.site_mu_k::<2>(
                                    2.0 * self.XA - 1.0,
                                    zeta2t0,
                                    zeta3t0,
                                    dit0,
                                    zeta2_k[k],
                                    zeta3_k[k],
                                )
                            }
                            AssocType::Type4C => {
                                4.0 * self.site_mu_k::<1>(
                                    self.XA, zeta2t0, zeta3t0, dit0, zeta2_k[k], zeta3_k[k],
                                )
                            }
                        };
                let vec: Vec<f64> = zeta2_k
                    .iter()
                    .zip(zeta3_k)
                    .enumerate()
                    .map(|(index, (&zeta2_k, &zeta3_k))| {
                        if k == index {
                            0.0
                        } else {
                            rho_num
                                * self.x
                                * match self.assoc_type {
                                    AssocType::Type1 => self.site_mu_k_g::<1>(
                                        self.XA, zeta2t0, zeta3t0, dit0, zeta2_k, zeta3_k,
                                    ),
                                    AssocType::Type2B => {
                                        2.0 * self.site_mu_k_g::<1>(
                                            self.XA, zeta2t0, zeta3t0, dit0, zeta2_k, zeta3_k,
                                        )
                                    }
                                    AssocType::Type3B => {
                                        2.0 * self.site_mu_k_g::<1>(
                                            self.XA, zeta2t0, zeta3t0, dit0, zeta2_k, zeta3_k,
                                        ) + self.site_mu_k_g::<2>(
                                            2.0 * self.XA - 1.0,
                                            zeta2t0,
                                            zeta3t0,
                                            dit0,
                                            zeta2_k,
                                            zeta3_k,
                                        )
                                    }
                                    AssocType::Type4C => {
                                        4.0 * self.site_mu_k_g::<1>(
                                            self.XA, zeta2t0, zeta3t0, dit0, zeta2_k, zeta3_k,
                                        )
                                    }
                                }
                        }
                    })
                    .collect();
                zeta3_k
                    .iter()
                    .zip(vec)
                    .enumerate()
                    .map(move |(index, (_, vec))| if k == index { val } else { vec })
            }
            pub fn t0d0(
                &mut self,
                temp: f64,
                rho_num: f64,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
            ) -> f64 {
                self.XA_flash(temp, rho_num, zeta2t0, zeta3t0, dit0);
                self.x
                    * match self.assoc_type {
                        AssocType::Type1 => self.site_t0d0::<1>(self.XA),
                        AssocType::Type2B => 2.0 * self.site_t0d0::<1>(self.XA),
                        AssocType::Type3B => {
                            2.0 * self.site_t0d0::<1>(self.XA)
                                + self.site_t0d0::<2>(2.0 * self.XA - 1.0)
                        }
                        AssocType::Type4C => 4.0 * self.site_t0d0::<1>(self.XA),
                    }
            }
            pub fn t0d1(
                &mut self,
                temp: f64,
                rho_num: f64,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
            ) -> f64 {
                self.XA_flash(temp, rho_num, zeta2t0, zeta3t0, dit0);
                self.x
                    * match self.assoc_type {
                        AssocType::Type1 => self.site_t0d1::<1>(self.XA, zeta2t0, zeta3t0, dit0),
                        AssocType::Type2B => {
                            2.0 * self.site_t0d1::<1>(self.XA, zeta2t0, zeta3t0, dit0)
                        }
                        AssocType::Type3B => {
                            2.0 * self.site_t0d1::<1>(self.XA, zeta2t0, zeta3t0, dit0)
                                + self.site_t0d1::<2>(2.0 * self.XA - 1.0, zeta2t0, zeta3t0, dit0)
                        }
                        AssocType::Type4C => {
                            4.0 * self.site_t0d1::<1>(self.XA, zeta2t0, zeta3t0, dit0)
                        }
                    }
            }
            pub fn t0d2(
                &mut self,
                temp: f64,
                rho_num: f64,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
            ) -> f64 {
                self.XA_flash(temp, rho_num, zeta2t0, zeta3t0, dit0);
                self.x
                    * match self.assoc_type {
                        AssocType::Type1 => self.site_t0d2::<1>(self.XA, zeta2t0, zeta3t0, dit0),
                        AssocType::Type2B => {
                            2.0 * self.site_t0d2::<1>(self.XA, zeta2t0, zeta3t0, dit0)
                        }
                        AssocType::Type3B => {
                            2.0 * self.site_t0d2::<1>(self.XA, zeta2t0, zeta3t0, dit0)
                                + self.site_t0d2::<2>(2.0 * self.XA - 1.0, zeta2t0, zeta3t0, dit0)
                        }
                        AssocType::Type4C => {
                            4.0 * self.site_t0d2::<1>(self.XA, zeta2t0, zeta3t0, dit0)
                        }
                    }
            }
            pub fn t0d3(
                &mut self,
                temp: f64,
                rho_num: f64,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
            ) -> f64 {
                self.XA_flash(temp, rho_num, zeta2t0, zeta3t0, dit0);
                self.x
                    * match self.assoc_type {
                        AssocType::Type1 => self.site_t0d3::<1>(self.XA, zeta2t0, zeta3t0, dit0),
                        AssocType::Type2B => {
                            2.0 * self.site_t0d3::<1>(self.XA, zeta2t0, zeta3t0, dit0)
                        }
                        AssocType::Type3B => {
                            2.0 * self.site_t0d3::<1>(self.XA, zeta2t0, zeta3t0, dit0)
                                + self.site_t0d3::<2>(2.0 * self.XA - 1.0, zeta2t0, zeta3t0, dit0)
                        }
                        AssocType::Type4C => {
                            4.0 * self.site_t0d3::<1>(self.XA, zeta2t0, zeta3t0, dit0)
                        }
                    }
            }
            pub fn t1d0(
                &mut self,
                temp: f64,
                rho_num: f64,
                (zeta2t0, zeta2t1): (f64, f64),
                (zeta3t0, zeta3t1): (f64, f64),
                (dit0, dit1): (f64, f64),
            ) -> f64 {
                self.XA_flash(temp, rho_num, zeta2t0, zeta3t0, dit0);
                self.x
                    * match self.assoc_type {
                        AssocType::Type1 => self.site_t1d0::<1>(
                            self.XA,
                            (zeta2t0, zeta2t1),
                            (zeta3t0, zeta3t1),
                            (dit0, dit1),
                        ),
                        AssocType::Type2B => {
                            2.0 * self.site_t1d0::<1>(
                                self.XA,
                                (zeta2t0, zeta2t1),
                                (zeta3t0, zeta3t1),
                                (dit0, dit1),
                            )
                        }
                        AssocType::Type3B => {
                            2.0 * self.site_t1d0::<1>(
                                self.XA,
                                (zeta2t0, zeta2t1),
                                (zeta3t0, zeta3t1),
                                (dit0, dit1),
                            ) + self.site_t1d0::<2>(
                                2.0 * self.XA - 1.0,
                                (zeta2t0, zeta2t1),
                                (zeta3t0, zeta3t1),
                                (dit0, dit1),
                            )
                        }
                        AssocType::Type4C => {
                            4.0 * self.site_t1d0::<1>(
                                self.XA,
                                (zeta2t0, zeta2t1),
                                (zeta3t0, zeta3t1),
                                (dit0, dit1),
                            )
                        }
                    }
            }
            pub fn t1d1(
                &mut self,
                temp: f64,
                rho_num: f64,
                (zeta2t0, zeta2t1): (f64, f64),
                (zeta3t0, zeta3t1): (f64, f64),
                (dit0, dit1): (f64, f64),
            ) -> f64 {
                self.XA_flash(temp, rho_num, zeta2t0, zeta3t0, dit0);
                self.x
                    * match self.assoc_type {
                        AssocType::Type1 => self.site_t1d1::<1>(
                            self.XA,
                            (zeta2t0, zeta2t1),
                            (zeta3t0, zeta3t1),
                            (dit0, dit1),
                        ),
                        AssocType::Type2B => {
                            2.0 * self.site_t1d1::<1>(
                                self.XA,
                                (zeta2t0, zeta2t1),
                                (zeta3t0, zeta3t1),
                                (dit0, dit1),
                            )
                        }
                        AssocType::Type3B => {
                            2.0 * self.site_t1d1::<1>(
                                self.XA,
                                (zeta2t0, zeta2t1),
                                (zeta3t0, zeta3t1),
                                (dit0, dit1),
                            ) + self.site_t1d1::<2>(
                                2.0 * self.XA - 1.0,
                                (zeta2t0, zeta2t1),
                                (zeta3t0, zeta3t1),
                                (dit0, dit1),
                            )
                        }
                        AssocType::Type4C => {
                            4.0 * self.site_t1d1::<1>(
                                self.XA,
                                (zeta2t0, zeta2t1),
                                (zeta3t0, zeta3t1),
                                (dit0, dit1),
                            )
                        }
                    }
            }
            pub fn t2d0(
                &mut self,
                temp: f64,
                rho_num: f64,
                (zeta2t0, zeta2t1, zeta2t2): (f64, f64, f64),
                (zeta3t0, zeta3t1, zeta3t2): (f64, f64, f64),
                (dit0, dit1, dit2): (f64, f64, f64),
            ) -> f64 {
                self.XA_flash(temp, rho_num, zeta2t0, zeta3t0, dit0);
                self.x
                    * match self.assoc_type {
                        AssocType::Type1 => self.site_t2d0::<1>(
                            self.XA,
                            (zeta2t0, zeta2t1, zeta2t2),
                            (zeta3t0, zeta3t1, zeta3t2),
                            (dit0, dit1, dit2),
                        ),
                        AssocType::Type2B => {
                            2.0 * self.site_t2d0::<1>(
                                self.XA,
                                (zeta2t0, zeta2t1, zeta2t2),
                                (zeta3t0, zeta3t1, zeta3t2),
                                (dit0, dit1, dit2),
                            )
                        }
                        AssocType::Type3B => {
                            2.0 * self.site_t2d0::<1>(
                                self.XA,
                                (zeta2t0, zeta2t1, zeta2t2),
                                (zeta3t0, zeta3t1, zeta3t2),
                                (dit0, dit1, dit2),
                            ) + self.site_t2d0::<2>(
                                2.0 * self.XA - 1.0,
                                (zeta2t0, zeta2t1, zeta2t2),
                                (zeta3t0, zeta3t1, zeta3t2),
                                (dit0, dit1, dit2),
                            )
                        }
                        AssocType::Type4C => {
                            4.0 * self.site_t2d0::<1>(
                                self.XA,
                                (zeta2t0, zeta2t1, zeta2t2),
                                (zeta3t0, zeta3t1, zeta3t2),
                                (dit0, dit1, dit2),
                            )
                        }
                    }
            }
        }
        impl $name {
            fn site_mu_k<const C: i32>(
                &mut self,
                x: f64,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
                zeta2_k: f64,
                zeta3_k: f64,
            ) -> f64 {
                self.site_x1::<C>(x) * self.x_mu_k(zeta2t0, zeta3t0, dit0, zeta2_k, zeta3_k)
            }
            fn site_mu_k_g<const C: i32>(
                &mut self,
                x: f64,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
                zeta2_k: f64,
                zeta3_k: f64,
            ) -> f64 {
                self.site_x1::<C>(x) * self.x_mu_k_g(zeta2t0, zeta3t0, dit0, zeta2_k, zeta3_k)
            }
            fn site_t0d0<const C: i32>(&mut self, x: f64) -> f64 {
                x.ln() - x / 2.0 + 0.5
            }
            fn site_t0d1<const C: i32>(
                &mut self,
                x: f64,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
            ) -> f64 {
                self.site_x1::<C>(x) * self.x_t0d1(zeta2t0, zeta3t0, dit0)
            }
            fn site_t0d2<const C: i32>(
                &mut self,
                x: f64,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
            ) -> f64 {
                self.site_x2::<C>(x) * self.x_t0d1(zeta2t0, zeta3t0, dit0).powi(2)
                    + self.site_x1::<C>(x) * self.x_t0d2(zeta2t0, zeta3t0, dit0)
            }
            fn site_t0d3<const C: i32>(
                &mut self,
                x: f64,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
            ) -> f64 {
                self.site_x3::<C>(x) * self.x_t0d1(zeta2t0, zeta3t0, dit0).powi(3)
                    + 3.0
                        * self.site_x2::<C>(x)
                        * self.x_t0d1(zeta2t0, zeta3t0, dit0)
                        * self.x_t0d2(zeta2t0, zeta3t0, dit0)
                    + self.site_x1::<C>(x) * self.x_t0d3(zeta2t0, zeta3t0, dit0)
            }
            fn site_t1d0<const C: i32>(
                &mut self,
                x: f64,
                (zeta2t0, zeta2t1): (f64, f64),
                (zeta3t0, zeta3t1): (f64, f64),
                (dit0, dit1): (f64, f64),
            ) -> f64 {
                self.site_x1::<C>(x)
                    * self.x_t1d0((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1))
            }
            fn site_t1d1<const C: i32>(
                &mut self,
                x: f64,
                (zeta2t0, zeta2t1): (f64, f64),
                (zeta3t0, zeta3t1): (f64, f64),
                (dit0, dit1): (f64, f64),
            ) -> f64 {
                self.site_x2::<C>(x)
                    * self.x_t1d0((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1))
                    * self.x_t0d1(zeta2t0, zeta3t0, dit0)
                    + self.site_x1::<C>(x)
                        * self.x_t1d1((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1))
            }
            fn site_t2d0<const C: i32>(
                &mut self,
                x: f64,
                (zeta2t0, zeta2t1, zeta2t2): (f64, f64, f64),
                (zeta3t0, zeta3t1, zeta3t2): (f64, f64, f64),
                (dit0, dit1, dit2): (f64, f64, f64),
            ) -> f64 {
                self.site_x2::<C>(x)
                    * self
                        .x_t1d0((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1))
                        .powi(2)
                    + self.site_x1::<C>(x)
                        * self.x_t2d0(
                            (zeta2t0, zeta2t1, zeta2t2),
                            (zeta3t0, zeta3t1, zeta3t2),
                            (dit0, dit1, dit2),
                        )
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
        }
        impl $name {
            fn x_mu_k(
                &mut self,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
                zeta2_k: f64,
                zeta3_k: f64,
            ) -> f64 {
                self.xt1() * self.t_mu_k(zeta2t0, zeta3t0, dit0, zeta2_k, zeta3_k)
            }
            fn x_mu_k_g(
                &mut self,
                zeta2t0: f64,
                zeta3t0: f64,
                dit0: f64,
                zeta2_k: f64,
                zeta3_k: f64,
            ) -> f64 {
                self.xt1() * self.t_mu_k_g(zeta2t0, zeta3t0, dit0, zeta2_k, zeta3_k)
            }
            fn x_t0d1(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
                if self.x_t0d1.0 != self.XA {
                    self.x_t0d1 = (self.XA, self.xt1() * self.t_t0d1(zeta2t0, zeta3t0, dit0))
                }
                self.x_t0d1.1
            }
            fn x_t0d2(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
                if self.x_t0d2.0 != self.XA {
                    self.x_t0d2 = (
                        self.XA,
                        self.xt2() * self.t_t0d1(zeta2t0, zeta3t0, dit0).powi(2)
                            + self.xt1() * self.t_t0d2(zeta2t0, zeta3t0, dit0),
                    )
                }
                self.x_t0d2.1
            }
            fn x_t0d3(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
                if self.x_t0d3.0 != self.XA {
                    self.x_t0d3 = (
                        self.XA,
                        self.xt3() * self.t_t0d1(zeta2t0, zeta3t0, dit0).powi(3)
                            + 3.0
                                * self.xt2()
                                * self.t_t0d1(zeta2t0, zeta3t0, dit0)
                                * self.t_t0d2(zeta2t0, zeta3t0, dit0)
                            + self.xt1() * self.t_t0d3(zeta2t0, zeta3t0, dit0),
                    )
                }
                self.x_t0d3.1
            }
            fn x_t1d0(
                &mut self,
                (zeta2t0, zeta2t1): (f64, f64),
                (zeta3t0, zeta3t1): (f64, f64),
                (dit0, dit1): (f64, f64),
            ) -> f64 {
                if self.x_t1d0.0 != self.XA {
                    self.x_t1d0 = (
                        self.XA,
                        self.xt1()
                            * self.t_t1d0((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1)),
                    )
                }
                self.x_t1d0.1
            }
            fn x_t1d1(
                &mut self,
                (zeta2t0, zeta2t1): (f64, f64),
                (zeta3t0, zeta3t1): (f64, f64),
                (dit0, dit1): (f64, f64),
            ) -> f64 {
                if self.x_t1d1.0 != self.XA {
                    self.x_t1d1 = (
                        self.XA,
                        self.xt2()
                            * self.t_t1d0((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1))
                            * self.t_t0d1(zeta2t0, zeta3t0, dit0)
                            + self.xt1()
                                * self.t_t1d1((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1)),
                    )
                }
                self.x_t1d1.1
            }
            fn x_t2d0(
                &mut self,
                (zeta2t0, zeta2t1, zeta2t2): (f64, f64, f64),
                (zeta3t0, zeta3t1, zeta3t2): (f64, f64, f64),
                (dit0, dit1, dit2): (f64, f64, f64),
            ) -> f64 {
                if self.x_t2d0.0 != self.XA {
                    self.x_t2d0 = (
                        self.XA,
                        self.xt2()
                            * self
                                .t_t1d0((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1))
                                .powi(2)
                            + self.xt1()
                                * self.t_t2d0(
                                    (zeta2t0, zeta2t1, zeta2t2),
                                    (zeta3t0, zeta3t1, zeta3t2),
                                    (dit0, dit1, dit2),
                                ),
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
        }
    };
}
fn_assoc!(AssocTerm); // fn_assoc!(AssocTerm);
impl AssocTerm {
    fn t_mu_k(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64, zeta2_k: f64, zeta3_k: f64) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0)
            * (self.kappa_AB_sigma3 * self.gii.t0d0(zeta2t0, zeta3t0, &[dit0])[0]
                + (self.x * self.kappa_AB_sigma3_dens)
                    * (zeta3_k / (1.0 - zeta3t0).powi(2)
                        + dit0
                            * (1.5 * zeta2_k / (1.0 - zeta3t0).powi(2)
                                + 3.0 * zeta3_k * zeta2t0 / (1.0 - zeta3t0).powi(3))
                        + dit0.powi(2)
                            * (zeta2_k * zeta2t0 / (1.0 - zeta3t0).powi(3)
                                + 1.5 * zeta3_k * zeta2t0.powi(2) / (1.0 - zeta3t0).powi(4))))
    }
    fn t_mu_k_g(
        &mut self,
        zeta2t0: f64,
        zeta3t0: f64,
        dit0: f64,
        zeta2_k: f64,
        zeta3_k: f64,
    ) -> f64 {
        self.x
            * self.kappa_AB_sigma3_dens
            * (self.epsilon_AB_temp.exp() - 1.0)
            * (zeta3_k / (1.0 - zeta3t0).powi(2)
                + dit0
                    * (1.5 * zeta2_k / (1.0 - zeta3t0).powi(2)
                        + 3.0 * zeta3_k * zeta2t0 / (1.0 - zeta3t0).powi(3))
                + dit0.powi(2)
                    * (zeta2_k * zeta2t0 / (1.0 - zeta3t0).powi(3)
                        + 1.5 * zeta3_k * zeta2t0.powi(2) / (1.0 - zeta3t0).powi(4)))
    }
    fn t_t0d0(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.t_t0d0.0 {
            self.t_t0d0 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * self.gii.t0d0(zeta2t0, zeta3t0, &[dit0])[0],
            )
        }
        self.t_t0d0.1
    }
    fn t_t0d1(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.t_t0d1.0 {
            self.t_t0d1 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii.t0d1(zeta2t0, zeta3t0, &[dit0])[0]
                        + self.gii.t0d0(zeta2t0, zeta3t0, &[dit0])[0]),
            )
        }
        self.t_t0d1.1
    }
    fn t_t0d2(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.t_t0d2.0 {
            self.t_t0d2 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii.t0d2(zeta2t0, zeta3t0, &[dit0])[0]
                        + 2.0 * self.gii.t0d1(zeta2t0, zeta3t0, &[dit0])[0]),
            )
        }
        self.t_t0d2.1
    }
    fn t_t0d3(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.t_t0d3.0 {
            self.t_t0d3 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii.t0d3(zeta2t0, zeta3t0, &[dit0])[0]
                        + 3.0 * self.gii.t0d2(zeta2t0, zeta3t0, &[dit0])[0]),
            )
        }
        self.t_t0d3.1
    }
    fn t_t1d0(
        &mut self,
        (zeta2t0, zeta2t1): (f64, f64),
        (zeta3t0, zeta3t1): (f64, f64),
        (dit0, dit1): (f64, f64),
    ) -> f64 {
        if zeta3t0 != self.t_t1d0.0 {
            self.t_t1d0 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * self.gii.t1d0(
                            (zeta2t0, zeta2t1),
                            (zeta3t0, zeta3t1),
                            (&[dit0], &[dit1]),
                        )[0]
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * self.gii.t0d0(zeta2t0, zeta3t0, &[dit0])[0]),
            )
        }
        self.t_t1d0.1
    }
    fn t_t1d1(
        &mut self,
        (zeta2t0, zeta2t1): (f64, f64),
        (zeta3t0, zeta3t1): (f64, f64),
        (dit0, dit1): (f64, f64),
    ) -> f64 {
        if zeta3t0 != self.t_t1d1.0 {
            self.t_t1d1 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii.t1d1(
                            (zeta2t0, zeta2t1),
                            (zeta3t0, zeta3t1),
                            (&[dit0], &[dit1]),
                        )[0] + self.gii.t1d0(
                            (zeta2t0, zeta2t1),
                            (zeta3t0, zeta3t1),
                            (&[dit0], &[dit1]),
                        )[0])
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii.t0d1(zeta2t0, zeta3t0, &[dit0])[0]
                                + self.gii.t0d0(zeta2t0, zeta3t0, &[dit0])[0])),
            )
        }
        self.t_t1d1.1
    }
    fn t_t2d0(
        &mut self,
        (zeta2t0, zeta2t1, zeta2t2): (f64, f64, f64),
        (zeta3t0, zeta3t1, zeta3t2): (f64, f64, f64),
        (dit0, dit1, dit2): (f64, f64, f64),
    ) -> f64 {
        if zeta3t0 != self.t_t2d0.0 {
            self.t_t2d0 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * self.gii.t2d0(
                            (zeta2t0, zeta2t1, zeta2t2),
                            (zeta3t0, zeta3t1, zeta3t2),
                            (&[dit0], &[dit1], &[dit2]),
                        )[0]
                        - 2.0
                            * self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                            * self.gii.t1d0(
                                (zeta2t0, zeta2t1),
                                (zeta3t0, zeta3t1),
                                (&[dit0], &[dit1]),
                            )[0]
                        + self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                            * (self.epsilon_AB_temp + 2.0)
                            * self.gii.t0d0(zeta2t0, zeta3t0, &[dit0])[0]),
            )
        }
        self.t_t2d0.1
    }
}
fn_assoc!(AssocGlyTerm); // fn_assoc!(AssocGlyTerm);
impl AssocGlyTerm {
    fn t_mu_k(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64, zeta2_k: f64, zeta3_k: f64) -> f64 {
        (self.epsilon_AB_temp.exp() - 1.0)
            * (self.kappa_AB_sigma3 * self.gii_gly_t0d0(zeta2t0, zeta3t0, dit0)
                + (self.x * self.kappa_AB_sigma3_dens)
                    * self.gii_gly_mu_k(zeta2t0, zeta3t0, dit0, zeta2_k, zeta3_k))
    }
    fn t_mu_k_g(
        &mut self,
        zeta2t0: f64,
        zeta3t0: f64,
        dit0: f64,
        zeta2_k: f64,
        zeta3_k: f64,
    ) -> f64 {
        self.x
            * self.kappa_AB_sigma3_dens
            * (self.epsilon_AB_temp.exp() - 1.0)
            * self.gii_gly_mu_k(zeta2t0, zeta3t0, dit0, zeta2_k, zeta3_k)
    }
    fn t_t0d0(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.t_t0d0.0 {
            self.t_t0d0 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * self.gii_gly_t0d0(zeta2t0, zeta3t0, dit0),
            )
        }
        self.t_t0d0.1
    }
    fn t_t0d1(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.t_t0d1.0 {
            self.t_t0d1 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii_gly_t0d1(zeta2t0, zeta3t0, dit0)
                        + self.gii_gly_t0d0(zeta2t0, zeta3t0, dit0)),
            )
        }
        self.t_t0d1.1
    }
    fn t_t0d2(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.t_t0d2.0 {
            self.t_t0d2 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii_gly_t0d2(zeta2t0, zeta3t0, dit0)
                        + 2.0 * self.gii_gly_t0d1(zeta2t0, zeta3t0, dit0)),
            )
        }
        self.t_t0d2.1
    }
    fn t_t0d3(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.t_t0d3.0 {
            self.t_t0d3 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii_gly_t0d3(zeta2t0, zeta3t0, dit0)
                        + 3.0 * self.gii_gly_t0d2(zeta2t0, zeta3t0, dit0)),
            )
        }
        self.t_t0d3.1
    }
    fn t_t1d0(
        &mut self,
        (zeta2t0, zeta2t1): (f64, f64),
        (zeta3t0, zeta3t1): (f64, f64),
        (dit0, dit1): (f64, f64),
    ) -> f64 {
        if zeta3t0 != self.t_t1d0.0 {
            self.t_t1d0 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * self.gii_gly_t1d0((zeta2t0, zeta2t1), (zeta3t0, zeta3t1), (dit0, dit1))
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * self.gii_gly_t0d0(zeta2t0, zeta3t0, dit0)),
            )
        }
        self.t_t1d0.1
    }
    fn t_t1d1(
        &mut self,
        (zeta2t0, zeta2t1): (f64, f64),
        (zeta3t0, zeta3t1): (f64, f64),
        (dit0, dit1): (f64, f64),
    ) -> f64 {
        if zeta3t0 != self.t_t1d1.0 {
            self.t_t1d1 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii_gly_t1d1(
                            (zeta2t0, zeta2t1),
                            (zeta3t0, zeta3t1),
                            (dit0, dit1),
                        ) + self.gii_gly_t1d0(
                            (zeta2t0, zeta2t1),
                            (zeta3t0, zeta3t1),
                            (dit0, dit1),
                        ))
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii_gly_t0d1(zeta2t0, zeta3t0, dit0)
                                + self.gii_gly_t0d0(zeta2t0, zeta3t0, dit0))),
            )
        }
        self.t_t1d1.1
    }
    fn t_t2d0(
        &mut self,
        (zeta2t0, zeta2t1, zeta2t2): (f64, f64, f64),
        (zeta3t0, zeta3t1, zeta3t2): (f64, f64, f64),
        (dit0, dit1, dit2): (f64, f64, f64),
    ) -> f64 {
        if zeta3t0 != self.t_t2d0.0 {
            self.t_t2d0 = (
                zeta3t0,
                self.x
                    * self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * self.gii_gly_t2d0(
                            (zeta2t0, zeta2t1, zeta2t2),
                            (zeta3t0, zeta3t1, zeta3t2),
                            (dit0, dit1, dit2),
                        )
                        - 2.0
                            * self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                            * self.gii_gly_t1d0(
                                (zeta2t0, zeta2t1),
                                (zeta3t0, zeta3t1),
                                (dit0, dit1),
                            )
                        + self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                            * (self.epsilon_AB_temp + 2.0)
                            * self.gii_gly_t0d0(zeta2t0, zeta3t0, dit0)),
            )
        }
        self.t_t2d0.1
    }
}
impl AssocGlyTerm {
    fn gii_gly_mu_k(
        &mut self,
        zeta2t0: f64,
        zeta3t0: f64,
        dit0: f64,
        zeta2_k: f64,
        zeta3_k: f64,
    ) -> f64 {
        self.c0 * zeta3_k / (1.0 - zeta3t0).powi(2)
            + (self.c1 * dit0)
                * (1.5 * zeta2_k / (1.0 - zeta3t0).powi(2)
                    + 3.0 * zeta3_k * zeta2t0 / (1.0 - zeta3t0).powi(3))
            + (self.c2 * dit0.powi(2))
                * (zeta2_k * zeta2t0 / (1.0 - zeta3t0).powi(3)
                    + 1.5 * zeta3_k * zeta2t0.powi(2) / (1.0 - zeta3t0).powi(4))
    }
    fn gii_gly_t0d0(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.g_t0d0.0 {
            self.g_t0d0 = (
                zeta3t0,
                self.c0 / (1.0 - zeta3t0)
                    + self.c1 * dit0 * 1.5 * zeta2t0 / (1.0 - zeta3t0).powi(2)
                    + self.c2 * dit0.powi(2) * 0.5 * zeta2t0.powi(2) / (1.0 - zeta3t0).powi(3),
            )
        }
        self.g_t0d0.1
    }
    fn gii_gly_t0d1(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.g_t0d1.0 {
            self.g_t0d1 = (
                zeta3t0,
                self.c0 * zeta3t0 / (1.0 - zeta3t0).powi(2)
                    + self.c1 * dit0 * 1.5 * zeta2t0 / (1.0 - zeta3t0).powi(3) * (1.0 + zeta3t0)
                    + self.c2 * dit0.powi(2) * 0.5 * zeta2t0.powi(2) / (1.0 - zeta3t0).powi(4)
                        * (2.0 + zeta3t0),
            )
        }
        self.g_t0d1.1
    }
    fn gii_gly_t0d2(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.g_t0d2.0 {
            self.g_t0d2 = (
                zeta3t0,
                self.c0 * 2.0 * zeta3t0.powi(2) / (1.0 - zeta3t0).powi(3)
                    + self.c1 * dit0 * 3.0 * zeta2t0 * zeta3t0 / (1.0 - zeta3t0).powi(4)
                        * (2.0 + zeta3t0)
                    + self.c2 * dit0.powi(2) * zeta2t0.powi(2) / (1.0 - zeta3t0).powi(5)
                        * (1.0 + 4.0 * zeta3t0 + zeta3t0.powi(2)),
            )
        }
        self.g_t0d2.1
    }
    fn gii_gly_t0d3(&mut self, zeta2t0: f64, zeta3t0: f64, dit0: f64) -> f64 {
        if zeta3t0 != self.g_t0d3.0 {
            self.g_t0d3 = (
                zeta3t0,
                self.c0 * 6.0 * zeta3t0.powi(3) / (1.0 - zeta3t0).powi(4)
                    + self.c1 * dit0 * 9.0 * zeta2t0 * zeta3t0.powi(2) / (1.0 - zeta3t0).powi(5)
                        * (3.0 + zeta3t0)
                    + self.c2 * dit0.powi(2) * 3.0 * zeta2t0.powi(2) * zeta3t0
                        / (1.0 - zeta3t0).powi(6)
                        * (3.0 + 6.0 * zeta3t0 + zeta3t0.powi(2)),
            )
        }
        self.g_t0d3.1
    }
    fn gii_gly_t1d0(
        &mut self,
        (zeta2t0, zeta2t1): (f64, f64),
        (zeta3t0, zeta3t1): (f64, f64),
        (dit0, dit1): (f64, f64),
    ) -> f64 {
        if zeta3t0 != self.g_t1d0.0 {
            self.g_t1d0 = (
                zeta3t0,
                self.c0 * zeta3t1 / (1.0 - zeta3t0).powi(2)
                    + self.c1 * dit1 * 1.5 * zeta2t0 / (1.0 - zeta3t0).powi(2)
                    + (self.c1 * dit0)
                        * (1.5 * zeta2t1 / (1.0 - zeta3t0).powi(2)
                            + 3.0 * zeta2t0 * zeta3t1 / (1.0 - zeta3t0).powi(3))
                    + self.c2 * dit0 * dit1 * zeta2t0.powi(2) / (1.0 - zeta3t0).powi(3)
                    + (self.c2 * dit0.powi(2))
                        * (zeta2t0 * zeta2t1 / (1.0 - zeta3t0).powi(3)
                            + 1.5 * zeta2t0.powi(2) * zeta3t1 / (1.0 - zeta3t0).powi(4)),
            )
        }
        self.g_t1d0.1
    }
    fn gii_gly_t1d1(
        &mut self,
        (zeta2t0, zeta2t1): (f64, f64),
        (zeta3t0, zeta3t1): (f64, f64),
        (dit0, dit1): (f64, f64),
    ) -> f64 {
        if zeta3t0 != self.g_t1d1.0 {
            self.g_t1d1 = (
                zeta3t0,
                self.c0 * zeta3t1 / (1.0 - zeta3t0).powi(3) * (1.0 + zeta3t0)
                    + self.c1 * dit1 * 1.5 * zeta2t0 / (1.0 - zeta3t0).powi(3) * (1.0 + zeta3t0)
                    + (self.c1 * dit0)
                        * (1.5 * zeta2t1 / (1.0 - zeta3t0).powi(3) * (1.0 + zeta3t0)
                            + 3.0 * zeta2t0 * zeta3t1 / (1.0 - zeta3t0).powi(4) * (2.0 + zeta3t0))
                    + self.c2 * dit0 * dit1 * zeta2t0.powi(2) / (1.0 - zeta3t0).powi(4)
                        * (2.0 + zeta3t0)
                    + (self.c2 * dit0.powi(2))
                        * (zeta2t0 * zeta2t1 / (1.0 - zeta3t0).powi(4) * (2.0 + zeta3t0)
                            + 1.5 * zeta2t0.powi(2) * zeta3t1 / (1.0 - zeta3t0).powi(5)
                                * (3.0 + zeta3t0)),
            )
        }
        self.g_t1d1.1
    }
    fn gii_gly_t2d0(
        &mut self,
        (zeta2t0, zeta2t1, zeta2t2): (f64, f64, f64),
        (zeta3t0, zeta3t1, zeta3t2): (f64, f64, f64),
        (dit0, dit1, dit2): (f64, f64, f64),
    ) -> f64 {
        if zeta3t0 != self.g_t2d0.0 {
            self.g_t2d0 = (
                zeta3t0,
                self.c0
                    * (zeta3t2 / (1.0 - zeta3t0).powi(2)
                        + 2.0 * zeta3t1.powi(2) / (1.0 - zeta3t0).powi(3))
                    + self.c1 * dit2 * 1.5 * zeta2t0 / (1.0 - zeta3t0).powi(2)
                    + (self.c1 * dit1)
                        * (3.0 * zeta2t1 / (1.0 - zeta3t0).powi(2)
                            + 6.0 * zeta2t0 * zeta3t1 / (1.0 - zeta3t0).powi(3))
                    + (self.c1 * dit0)
                        * (1.5 * zeta2t2 / (1.0 - zeta3t0).powi(2)
                            + 6.0 * zeta2t1 * zeta3t1 / (1.0 - zeta3t0).powi(3)
                            + 3.0 * zeta2t0 * zeta3t2 / (1.0 - zeta3t0).powi(3)
                            + 9.0 * zeta2t0 * zeta3t1.powi(2) / (1.0 - zeta3t0).powi(4))
                    + self.c2 * (dit1.powi(2) + dit0 * dit2) * zeta2t0.powi(2)
                        / (1.0 - zeta3t0).powi(3)
                    + (self.c2 * dit0 * dit1)
                        * (4.0 * zeta2t0 * zeta2t1 / (1.0 - zeta3t0).powi(3)
                            + 6.0 * zeta2t0.powi(2) * zeta3t1 / (1.0 - zeta3t0).powi(4))
                    + (self.c2 * dit0.powi(2))
                        * (zeta2t1.powi(2) / (1.0 - zeta3t0).powi(3)
                            + zeta2t0 * zeta2t2 / (1.0 - zeta3t0).powi(3)
                            + 6.0 * zeta2t0 * zeta2t1 * zeta3t1 / (1.0 - zeta3t0).powi(4)
                            + 1.5 * zeta2t0.powi(2) * zeta3t2 / (1.0 - zeta3t0).powi(4)
                            + 6.0 * zeta2t0.powi(2) * zeta3t1.powi(2) / (1.0 - zeta3t0).powi(5)),
            )
        }
        self.g_t2d0.1
    }
}
