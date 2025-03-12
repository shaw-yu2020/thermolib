use super::GiiTerm;
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
            gii: GiiTerm::new(),
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
enum AssocType {
    Type1,
    Type2B,
    Type3B,
    Type4C,
}
impl AssocTerm {
    pub fn t0d0(&mut self, temp: f64, rho_num: f64, eta: f64) -> f64 {
        self.XA_flash(temp, rho_num, eta);
        match self.assoc_type {
            AssocType::Type1 => self.site_t0d0::<1>(self.XA),
            AssocType::Type2B => 2.0 * self.site_t0d0::<1>(self.XA),
            AssocType::Type3B => {
                2.0 * self.site_t0d0::<1>(self.XA) + self.site_t0d0::<2>(2.0 * self.XA - 1.0)
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
impl AssocTerm {
    fn site_t0d0<const C: i32>(&mut self, x: f64) -> f64 {
        x.ln() - x / 2.0 + 0.5
    }
    fn site_t0d1<const C: i32>(&mut self, x: f64, eta: f64) -> f64 {
        self.site_x1::<C>(x) * self.x_t0d1(eta)
    }
    fn site_t0d2<const C: i32>(&mut self, x: f64, eta: f64) -> f64 {
        self.site_x2::<C>(x) * self.x_t0d1(eta).powi(2) + self.site_x1::<C>(x) * self.x_t0d2(eta)
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
impl AssocTerm {
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
impl AssocTerm {
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
impl AssocTerm {
    fn xt1(&mut self) -> f64 {
        if self.xt1.0 != self.XA {
            self.xt1 = (
                self.XA,
                match self.assoc_type {
                    AssocType::Type1 | AssocType::Type2B => self.XA.powi(3) / (self.XA - 2.0),
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
                            * (4.0 * self.XA.powi(3) - 12.0 * self.XA.powi(2) + 6.0 * self.XA - 1.0)
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
                            * (self.XA.powi(3) - 9.0 * self.XA.powi(2) + 29.0 * self.XA - 35.0)
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
                            * (self.XA.powi(3) - 9.0 * self.XA.powi(2) + 29.0 * self.XA - 35.0)
                    }
                },
            )
        }
        self.xt4.1
    }
}
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
