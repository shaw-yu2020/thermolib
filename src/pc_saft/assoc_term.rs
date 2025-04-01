use super::GiiTerm;
enum AssocType {
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
    // cached variables
    xt1: (f64, f64),
    xt2: (f64, f64),
    xt3: (f64, f64),
    xt4: (f64, f64),
    x_t0d1: (f64, f64),
    x_t0d2: (f64, f64),
    x_t0d3: (f64, f64),
    x_t0d4: (f64, f64),
    x_t1d0: (f64, f64),
    x_t1d1: (f64, f64),
    x_t1d2: (f64, f64),
    x_t1d3: (f64, f64),
    x_t2d0: (f64, f64),
    t_t0d0: (f64, f64),
    t_t0d1: (f64, f64),
    t_t0d2: (f64, f64),
    t_t0d3: (f64, f64),
    t_t0d4: (f64, f64),
    t_t1d0: (f64, f64),
    t_t1d1: (f64, f64),
    t_t1d2: (f64, f64),
    t_t1d3: (f64, f64),
    t_t2d0: (f64, f64),
}
#[allow(non_snake_case)]
impl AssocTerm {
    pub fn new_1_term(kappa_AB_sigma3: f64, epsilon_AB: f64) -> Self {
        Self {
            gii: GiiTerm::new(0.0),
            assoc_type: AssocType::Type1,
            kappa_AB_sigma3_dens: 0.0,
            epsilon_AB_temp: 0.0,
            kappa_AB_sigma3,
            epsilon_AB,
            dens: 0.0,
            temp: 0.0,
            XA: 1.0,
            // cached variables
            xt1: (0.0, 0.0),
            xt2: (0.0, 0.0),
            xt3: (0.0, 0.0),
            xt4: (0.0, 0.0),
            x_t0d1: (0.0, 0.0),
            x_t0d2: (0.0, 0.0),
            x_t0d3: (0.0, 0.0),
            x_t0d4: (0.0, 0.0),
            x_t1d0: (0.0, 0.0),
            x_t1d1: (0.0, 0.0),
            x_t1d2: (0.0, 0.0),
            x_t1d3: (0.0, 0.0),
            x_t2d0: (0.0, 0.0),
            t_t0d0: (0.0, 0.0),
            t_t0d1: (0.0, 0.0),
            t_t0d2: (0.0, 0.0),
            t_t0d3: (0.0, 0.0),
            t_t0d4: (0.0, 0.0),
            t_t1d0: (0.0, 0.0),
            t_t1d1: (0.0, 0.0),
            t_t1d2: (0.0, 0.0),
            t_t1d3: (0.0, 0.0),
            t_t2d0: (0.0, 0.0),
        }
    }
    pub fn new_2B_term(kappa_AB_sigma3: f64, epsilon_AB: f64) -> Self {
        Self {
            assoc_type: AssocType::Type2B,
            ..AssocTerm::new_1_term(kappa_AB_sigma3, epsilon_AB)
        }
    }
    pub fn new_3B_term(kappa_AB_sigma3: f64, epsilon_AB: f64) -> Self {
        Self {
            assoc_type: AssocType::Type3B,
            ..AssocTerm::new_1_term(kappa_AB_sigma3, epsilon_AB)
        }
    }
    pub fn new_4C_term(kappa_AB_sigma3: f64, epsilon_AB: f64) -> Self {
        Self {
            assoc_type: AssocType::Type4C,
            ..AssocTerm::new_1_term(kappa_AB_sigma3, epsilon_AB)
        }
    }
    fn XA_flash(&mut self, temp: f64, rho_num: f64, eta: f64) {
        if temp != self.temp || rho_num != self.dens {
            (self.temp, self.epsilon_AB_temp) = (temp, self.epsilon_AB / temp);
            (self.dens, self.kappa_AB_sigma3_dens) = (rho_num, self.kappa_AB_sigma3 * rho_num);
            self.XA = match self.assoc_type {
                AssocType::Type1 | AssocType::Type2B => {
                    (-1.0 + (1.0 + 4.0 * self.t_t0d0(eta)).sqrt()) / (2.0 * self.t_t0d0(eta))
                }
                AssocType::Type3B => {
                    (-(1.0 - self.t_t0d0(eta))
                        + ((1.0 - self.t_t0d0(eta)).powi(2) + 8.0 * self.t_t0d0(eta)).sqrt())
                        / (4.0 * self.t_t0d0(eta))
                }
                AssocType::Type4C => {
                    (-1.0 + (1.0 + 8.0 * self.t_t0d0(eta)).sqrt()) / (4.0 * self.t_t0d0(eta))
                }
            };
        }
    }
}
fn_assoc!(AssocTerm);
impl AssocTerm {
    fn t_t0d0(&mut self, eta: f64) -> f64 {
        if eta != self.t_t0d0.0 {
            self.t_t0d0 = (
                eta,
                self.kappa_AB_sigma3_dens * (self.epsilon_AB_temp.exp() - 1.0) * self.gii.t0d0(eta),
            )
        }
        self.t_t0d0.1
    }
    fn t_t0d1(&mut self, eta: f64) -> f64 {
        if eta != self.t_t0d1.0 {
            self.t_t0d1 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii.t0d1(eta) + self.gii.t0d0(eta)),
            )
        }
        self.t_t0d1.1
    }
    fn t_t0d2(&mut self, eta: f64) -> f64 {
        if eta != self.t_t0d2.0 {
            self.t_t0d2 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii.t0d2(eta) + 2.0 * self.gii.t0d1(eta)),
            )
        }
        self.t_t0d2.1
    }
    fn t_t0d3(&mut self, eta: f64) -> f64 {
        if eta != self.t_t0d3.0 {
            self.t_t0d3 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii.t0d3(eta) + 3.0 * self.gii.t0d2(eta)),
            )
        }
        self.t_t0d3.1
    }
    fn t_t0d4(&mut self, eta: f64) -> f64 {
        if eta != self.t_t0d4.0 {
            self.t_t0d4 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii.t0d4(eta) + 4.0 * self.gii.t0d3(eta)),
            )
        }
        self.t_t0d4.1
    }
    fn t_t1d0(&mut self, eta: f64, eta1: f64) -> f64 {
        if eta != self.t_t1d0.0 {
            self.t_t1d0 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0) * self.gii.t1d0(eta, eta1)
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp) * self.gii.t0d0(eta)),
            )
        }
        self.t_t1d0.1
    }
    fn t_t1d1(&mut self, eta: f64, eta1: f64) -> f64 {
        if eta != self.t_t1d1.0 {
            self.t_t1d1 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii.t1d1(eta, eta1) + self.gii.t1d0(eta, eta1))
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii.t0d1(eta) + self.gii.t0d0(eta))),
            )
        }
        self.t_t1d1.1
    }
    fn t_t1d2(&mut self, eta: f64, eta1: f64) -> f64 {
        if eta != self.t_t1d2.0 {
            self.t_t1d2 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii.t1d2(eta, eta1) + 2.0 * self.gii.t1d1(eta, eta1))
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii.t0d2(eta) + 2.0 * self.gii.t0d1(eta))),
            )
        }
        self.t_t1d2.1
    }
    fn t_t1d3(&mut self, eta: f64, eta1: f64) -> f64 {
        if eta != self.t_t1d3.0 {
            self.t_t1d3 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii.t1d3(eta, eta1) + 3.0 * self.gii.t1d2(eta, eta1))
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii.t0d3(eta) + 3.0 * self.gii.t0d2(eta))),
            )
        }
        self.t_t1d3.1
    }
    fn t_t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        if eta != self.t_t2d0.0 {
            self.t_t2d0 = (
                eta,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0) * self.gii.t2d0(eta, eta1, eta2)
                        - 2.0
                            * self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                            * self.gii.t1d0(eta, eta1)
                        + self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                            * (self.epsilon_AB_temp + 2.0)
                            * self.gii.t0d0(eta)),
            )
        }
        self.t_t2d0.1
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
    // parameters for gly
    c0: f64,
    c1: f64,
    c2: f64,
    // cached variables
    xt1: (f64, f64),
    xt2: (f64, f64),
    xt3: (f64, f64),
    xt4: (f64, f64),
    x_t0d1: (f64, f64),
    x_t0d2: (f64, f64),
    x_t0d3: (f64, f64),
    x_t0d4: (f64, f64),
    x_t1d0: (f64, f64),
    x_t1d1: (f64, f64),
    x_t1d2: (f64, f64),
    x_t1d3: (f64, f64),
    x_t2d0: (f64, f64),
    t_t0d0: (f64, f64),
    t_t0d1: (f64, f64),
    t_t0d2: (f64, f64),
    t_t0d3: (f64, f64),
    t_t0d4: (f64, f64),
    t_t1d0: (f64, f64),
    t_t1d1: (f64, f64),
    t_t1d2: (f64, f64),
    t_t1d3: (f64, f64),
    t_t2d0: (f64, f64),
    g_t0d0: (f64, f64),
    g_t0d1: (f64, f64),
    g_t0d2: (f64, f64),
    g_t0d3: (f64, f64),
    g_t0d4: (f64, f64),
    g_t1d0: (f64, f64),
    g_t1d1: (f64, f64),
    g_t1d2: (f64, f64),
    g_t1d3: (f64, f64),
    g_t2d0: (f64, f64),
}
#[allow(non_snake_case)]
impl AssocGlyTerm {
    pub fn new_1_term(kappa_AB_sigma3: f64, epsilon_AB: f64, c0: f64, c1: f64, c2: f64) -> Self {
        Self {
            assoc_type: AssocType::Type1,
            kappa_AB_sigma3_dens: 0.0,
            epsilon_AB_temp: 0.0,
            kappa_AB_sigma3,
            epsilon_AB,
            dens: 0.0,
            temp: 0.0,
            XA: 1.0,
            // parameters for gly
            c0,
            c1: -2.0 * c0 + 1.5 * c1,
            c2: c0 - 1.5 * c1 + 0.5 * c2,
            // cached variables
            xt1: (0.0, 0.0),
            xt2: (0.0, 0.0),
            xt3: (0.0, 0.0),
            xt4: (0.0, 0.0),
            x_t0d1: (0.0, 0.0),
            x_t0d2: (0.0, 0.0),
            x_t0d3: (0.0, 0.0),
            x_t0d4: (0.0, 0.0),
            x_t1d0: (0.0, 0.0),
            x_t1d1: (0.0, 0.0),
            x_t1d2: (0.0, 0.0),
            x_t1d3: (0.0, 0.0),
            x_t2d0: (0.0, 0.0),
            t_t0d0: (0.0, 0.0),
            t_t0d1: (0.0, 0.0),
            t_t0d2: (0.0, 0.0),
            t_t0d3: (0.0, 0.0),
            t_t0d4: (0.0, 0.0),
            t_t1d0: (0.0, 0.0),
            t_t1d1: (0.0, 0.0),
            t_t1d2: (0.0, 0.0),
            t_t1d3: (0.0, 0.0),
            t_t2d0: (0.0, 0.0),
            g_t0d0: (0.0, 0.0),
            g_t0d1: (0.0, 0.0),
            g_t0d2: (0.0, 0.0),
            g_t0d3: (0.0, 0.0),
            g_t0d4: (0.0, 0.0),
            g_t1d0: (0.0, 0.0),
            g_t1d1: (0.0, 0.0),
            g_t1d2: (0.0, 0.0),
            g_t1d3: (0.0, 0.0),
            g_t2d0: (0.0, 0.0),
        }
    }
    pub fn new_2B_term(kappa_AB_sigma3: f64, epsilon_AB: f64, c0: f64, c1: f64, c2: f64) -> Self {
        Self {
            assoc_type: AssocType::Type2B,
            ..AssocGlyTerm::new_1_term(kappa_AB_sigma3, epsilon_AB, c0, c1, c2)
        }
    }
    pub fn new_3B_term(kappa_AB_sigma3: f64, epsilon_AB: f64, c0: f64, c1: f64, c2: f64) -> Self {
        Self {
            assoc_type: AssocType::Type3B,
            ..AssocGlyTerm::new_1_term(kappa_AB_sigma3, epsilon_AB, c0, c1, c2)
        }
    }
    pub fn new_4C_term(kappa_AB_sigma3: f64, epsilon_AB: f64, c0: f64, c1: f64, c2: f64) -> Self {
        Self {
            assoc_type: AssocType::Type4C,
            ..AssocGlyTerm::new_1_term(kappa_AB_sigma3, epsilon_AB, c0, c1, c2)
        }
    }
    fn XA_flash(&mut self, temp: f64, rho_num: f64, eta: f64) {
        if temp != self.temp || rho_num != self.dens {
            (self.temp, self.epsilon_AB_temp) = (temp, self.epsilon_AB / temp);
            (self.dens, self.kappa_AB_sigma3_dens) = (rho_num, self.kappa_AB_sigma3 * rho_num);
            self.XA = match self.assoc_type {
                AssocType::Type1 | AssocType::Type2B => {
                    (-1.0 + (1.0 + 4.0 * self.t_t0d0(eta)).sqrt()) / (2.0 * self.t_t0d0(eta))
                }
                AssocType::Type3B => {
                    (-(1.0 - self.t_t0d0(eta))
                        + ((1.0 - self.t_t0d0(eta)).powi(2) + 8.0 * self.t_t0d0(eta)).sqrt())
                        / (4.0 * self.t_t0d0(eta))
                }
                AssocType::Type4C => {
                    (-1.0 + (1.0 + 8.0 * self.t_t0d0(eta)).sqrt()) / (4.0 * self.t_t0d0(eta))
                }
            };
        }
    }
}
fn_assoc!(AssocGlyTerm);
impl AssocGlyTerm {
    fn t_t0d0(&mut self, eta0: f64) -> f64 {
        if eta0 != self.t_t0d0.0 {
            self.t_t0d0 = (
                eta0,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * self.gii_gly_t0d0(eta0),
            )
        }
        self.t_t0d0.1
    }
    fn t_t0d1(&mut self, eta0: f64) -> f64 {
        if eta0 != self.t_t0d1.0 {
            self.t_t0d1 = (
                eta0,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii_gly_t0d1(eta0) + self.gii_gly_t0d0(eta0)),
            )
        }
        self.t_t0d1.1
    }
    fn t_t0d2(&mut self, eta0: f64) -> f64 {
        if eta0 != self.t_t0d2.0 {
            self.t_t0d2 = (
                eta0,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii_gly_t0d2(eta0) + 2.0 * self.gii_gly_t0d1(eta0)),
            )
        }
        self.t_t0d2.1
    }
    fn t_t0d3(&mut self, eta0: f64) -> f64 {
        if eta0 != self.t_t0d3.0 {
            self.t_t0d3 = (
                eta0,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii_gly_t0d3(eta0) + 3.0 * self.gii_gly_t0d2(eta0)),
            )
        }
        self.t_t0d3.1
    }
    fn t_t0d4(&mut self, eta0: f64) -> f64 {
        if eta0 != self.t_t0d4.0 {
            self.t_t0d4 = (
                eta0,
                self.kappa_AB_sigma3_dens
                    * (self.epsilon_AB_temp.exp() - 1.0)
                    * (self.gii_gly_t0d4(eta0) + 4.0 * self.gii_gly_t0d3(eta0)),
            )
        }
        self.t_t0d4.1
    }
    fn t_t1d0(&mut self, eta0: f64, eta1: f64) -> f64 {
        if eta0 != self.t_t1d0.0 {
            self.t_t1d0 = (
                eta0,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0) * self.gii_gly_t1d0(eta0, eta1)
                        - self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                            * self.gii_gly_t0d0(eta0)),
            )
        }
        self.t_t1d0.1
    }
    fn t_t1d1(&mut self, eta0: f64, eta1: f64) -> f64 {
        if eta0 != self.t_t1d1.0 {
            self.t_t1d1 = (
                eta0,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii_gly_t1d1(eta0, eta1) + self.gii_gly_t1d0(eta0, eta1))
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii_gly_t0d1(eta0) + self.gii_gly_t0d0(eta0))),
            )
        }
        self.t_t1d1.1
    }
    fn t_t1d2(&mut self, eta0: f64, eta1: f64) -> f64 {
        if eta0 != self.t_t1d2.0 {
            self.t_t1d2 = (
                eta0,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii_gly_t1d2(eta0, eta1) + 2.0 * self.gii_gly_t1d1(eta0, eta1))
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii_gly_t0d2(eta0) + 2.0 * self.gii_gly_t0d1(eta0))),
            )
        }
        self.t_t1d2.1
    }
    fn t_t1d3(&mut self, eta0: f64, eta1: f64) -> f64 {
        if eta0 != self.t_t1d3.0 {
            self.t_t1d3 = (
                eta0,
                self.kappa_AB_sigma3_dens
                    * ((self.epsilon_AB_temp.exp() - 1.0)
                        * (self.gii_gly_t1d3(eta0, eta1) + 3.0 * self.gii_gly_t1d2(eta0, eta1))
                        - (self.epsilon_AB_temp.exp() * self.epsilon_AB_temp)
                            * (self.gii_gly_t0d3(eta0) + 3.0 * self.gii_gly_t0d2(eta0))),
            )
        }
        self.t_t1d3.1
    }
    fn t_t2d0(&mut self, eta0: f64, eta1: f64, eta2: f64) -> f64 {
        if eta0 != self.t_t2d0.0 {
            self.t_t2d0 = (
                eta0,
                self.kappa_AB_sigma3_dens
                    * (self.gii_gly_t2d0(eta0, eta1, eta2) * (self.epsilon_AB_temp.exp() - 1.0)
                        - 2.0
                            * self.gii_gly_t1d0(eta0, eta1)
                            * self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                        + self.gii_gly_t0d0(eta0)
                            * self.epsilon_AB_temp.exp()
                            * self.epsilon_AB_temp
                            * (self.epsilon_AB_temp + 2.0)),
            )
        }
        self.t_t2d0.1
    }
}
impl AssocGlyTerm {
    fn gii_gly_t0d0(&mut self, eta0: f64) -> f64 {
        if eta0 != self.g_t0d0.0 {
            self.g_t0d0 = (
                eta0,
                (1.0 - eta0).powi(3).recip() * (self.c0 + self.c1 * eta0 + self.c2 * eta0.powi(2)),
            )
        }
        self.g_t0d0.1
    }
    fn gii_gly_t0d1(&mut self, eta0: f64) -> f64 {
        if eta0 != self.g_t0d1.0 {
            self.g_t0d1 = (
                eta0,
                eta0 / (1.0 - eta0).powi(4)
                    * ((3.0 * self.c0 + self.c1)
                        + (2.0 * self.c1 + 2.0 * self.c2) * eta0
                        + self.c2 * eta0.powi(2)),
            )
        }
        self.g_t0d1.1
    }
    fn gii_gly_t0d2(&mut self, eta0: f64) -> f64 {
        if eta0 != self.g_t0d2.0 {
            self.g_t0d2 = (
                eta0,
                2.0 * eta0.powi(2) / (1.0 - eta0).powi(5)
                    * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                        + (3.0 * self.c1 + 4.0 * self.c2) * eta0
                        + self.c2 * eta0.powi(2)),
            )
        }
        self.g_t0d2.1
    }
    fn gii_gly_t0d3(&mut self, eta0: f64) -> f64 {
        if eta0 != self.g_t0d3.0 {
            self.g_t0d3 = (
                eta0,
                6.0 * eta0.powi(3) / (1.0 - eta0).powi(6)
                    * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                        + (4.0 * self.c1 + 6.0 * self.c2) * eta0
                        + self.c2 * eta0.powi(2)),
            )
        }
        self.g_t0d3.1
    }
    fn gii_gly_t0d4(&mut self, eta0: f64) -> f64 {
        if eta0 != self.g_t0d4.0 {
            self.g_t0d4 = (
                eta0,
                24.0 * eta0.powi(4) / (1.0 - eta0).powi(7)
                    * ((15.0 * self.c0 + 10.0 * self.c1 + 6.0 * self.c2)
                        + (5.0 * self.c1 + 8.0 * self.c2) * eta0
                        + self.c2 * eta0.powi(2)),
            )
        }
        self.g_t0d4.1
    }
    fn gii_gly_t1d0(&mut self, eta0: f64, eta1: f64) -> f64 {
        if eta0 != self.g_t1d0.0 {
            self.g_t1d0 = (
                eta0,
                eta1 / (1.0 - eta0).powi(4)
                    * ((3.0 * self.c0 + self.c1)
                        + (2.0 * self.c1 + 2.0 * self.c2) * eta0
                        + self.c2 * eta0.powi(2)),
            )
        }
        self.g_t1d0.1
    }
    fn gii_gly_t1d1(&mut self, eta0: f64, eta1: f64) -> f64 {
        if eta0 != self.g_t1d1.0 {
            self.g_t1d1 = (
                eta0,
                eta1 / (1.0 - eta0).powi(4)
                    * ((3.0 * self.c0 + self.c1)
                        + (2.0 * self.c1 + 2.0 * self.c2) * eta0
                        + self.c2 * eta0.powi(2))
                    + 2.0 * eta1 * eta0 / (1.0 - eta0).powi(5)
                        * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                            + (3.0 * self.c1 + 4.0 * self.c2) * eta0
                            + self.c2 * eta0.powi(2)),
            )
        }
        self.g_t1d1.1
    }
    fn gii_gly_t1d2(&mut self, eta0: f64, eta1: f64) -> f64 {
        if eta0 != self.g_t1d2.0 {
            self.g_t1d2 = (
                eta0,
                4.0 * eta1 * eta0 / (1.0 - eta0).powi(5)
                    * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                        + (3.0 * self.c1 + 4.0 * self.c2) * eta0
                        + self.c2 * eta0.powi(2))
                    + 6.0 * eta1 * eta0.powi(2) / (1.0 - eta0).powi(6)
                        * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                            + (4.0 * self.c1 + 6.0 * self.c2) * eta0
                            + self.c2 * eta0.powi(2)),
            )
        }
        self.g_t1d2.1
    }
    fn gii_gly_t1d3(&mut self, eta0: f64, eta1: f64) -> f64 {
        if eta0 != self.g_t1d3.0 {
            self.g_t1d3 = (
                eta0,
                18.0 * eta1 * eta0.powi(2) / (1.0 - eta0).powi(6)
                    * ((10.0 * self.c0 + 6.0 * self.c1 + 3.0 * self.c2)
                        + (4.0 * self.c1 + 6.0 * self.c2) * eta0
                        + self.c2 * eta0.powi(2))
                    + 24.0 * eta1 * eta0.powi(3) / (1.0 - eta0).powi(7)
                        * ((15.0 * self.c0 + 10.0 * self.c1 + 6.0 * self.c2)
                            + (5.0 * self.c1 + 8.0 * self.c2) * eta0
                            + self.c2 * eta0.powi(2)),
            )
        }
        self.g_t1d3.1
    }
    fn gii_gly_t2d0(&mut self, eta0: f64, eta1: f64, eta2: f64) -> f64 {
        if eta0 != self.g_t2d0.0 {
            self.g_t2d0 = (
                eta0,
                eta2 / (1.0 - eta0).powi(4)
                    * ((3.0 * self.c0 + self.c1)
                        + (2.0 * self.c1 + 2.0 * self.c2) * eta0
                        + self.c2 * eta0.powi(2))
                    + 2.0 * eta1.powi(2) / (1.0 - eta0).powi(5)
                        * ((6.0 * self.c0 + 3.0 * self.c1 + self.c2)
                            + (3.0 * self.c1 + 4.0 * self.c2) * eta0
                            + self.c2 * eta0.powi(2)),
            )
        }
        self.g_t2d0.1
    }
}
