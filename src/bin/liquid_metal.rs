use eframe::egui;
use thermolib::LiquidMetal;
fn main() {
    let _ = eframe::run_native(
        "LiquidMetalApp",
        eframe::NativeOptions::default(),
        Box::new(|_| Box::<LiquidMetalApp>::default()),
    );
}
struct LiquidMetalApp {
    liquid_metal: LiquidMetal,
    metal: String,
    temperature: String,
    rho: String,
    eta: String,
    lambda: String,
}
impl Default for LiquidMetalApp {
    fn default() -> Self {
        LiquidMetalApp {
            liquid_metal: LiquidMetal::new_metal("Si").unwrap(),
            metal: String::from("Si"),
            temperature: String::from("1800.0"),
            rho: String::new(),
            eta: String::new(),
            lambda: String::new(),
        }
    }
}
impl eframe::App for LiquidMetalApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let mut no_err: bool = true;
        let temperature: f64 = match self.temperature.parse::<f64>() {
            Ok(t) => t,
            Err(_) => {
                self.rho = String::from("temperature is not a number.");
                self.eta = String::from("temperature is not a number.");
                self.lambda = String::from("temperature is not a number.");
                no_err = false;
                0.0
            }
        };
        if no_err {
            match thermolib::LiquidMetal::new_metal(&self.metal) {
                Ok(metal) => self.liquid_metal = metal,
                Err(_) => {
                    self.rho = String::from("there is no metal in database.");
                    self.eta = String::from("there is no metal in database.");
                    self.lambda = String::from("there is no metal in database.");
                    no_err = false;
                }
            }
        }
        if no_err {
            self.rho = match self.liquid_metal.calc_rho(temperature) {
                Ok(rho) => rho.to_string(),
                Err(_) => String::from("out of range"),
            };
            self.eta = match self.liquid_metal.calc_eta(temperature) {
                Ok(eta) => eta.to_string(),
                Err(_) => String::from("out of range"),
            };
            self.lambda = match self.liquid_metal.calc_lambda(temperature) {
                Ok(lambda) => lambda.to_string(),
                Err(_) => String::from("out of range"),
            };
        }
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                let _ = ui.button(" App Visuals: ");
                if ui.button(" Light ").clicked() {
                    ctx.set_visuals(egui::Visuals::light());
                }
                if ui.button(" Dark ").clicked() {
                    ctx.set_visuals(egui::Visuals::dark());
                }
            });
            ui.horizontal(|ui| {
                ui.label(" metal: ");
                ui.text_edit_singleline(&mut self.metal);
            });
            ui.horizontal(|ui| {
                ui.label(" temperature [K]: ");
                ui.text_edit_singleline(&mut self.temperature);
            });
            ui.horizontal(|ui| {
                ui.label(" density [kg/m3]: ");
                ui.label(&self.rho);
            });
            ui.horizontal(|ui| {
                ui.label(" viscosity [mPa*s]: ");
                ui.label(&self.eta);
            });
            ui.horizontal(|ui| {
                ui.label(" thermal conductivity [W/m/K]: ");
                ui.label(&self.lambda);
            })
        });
    }
}
