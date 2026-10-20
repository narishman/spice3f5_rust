use crate::circuit::Circuit;

pub mod bjt;
pub mod capacitor;
pub mod devsup;
pub mod diode;
pub mod inductor;
pub mod isource;
pub mod mos1;
pub mod resistor;
pub mod vsource;

pub trait Device {
    fn setup(&mut self, circuit: &mut Circuit) -> Result<(), String>;
    fn load(&mut self, circuit: &mut Circuit) -> Result<(), String>;
}

pub enum DeviceEnum {
    Resistor(resistor::Resistor),
    VoltageSource(vsource::VoltageSource),
    CurrentSource(isource::CurrentSource),
    Capacitor(capacitor::Capacitor),
    Inductor(inductor::Inductor),
    Diode(diode::Diode),
    Bjt(bjt::Bjt),
    Mos1(mos1::Mos1),
}

impl Device for DeviceEnum {
    fn setup(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        match self {
            DeviceEnum::Resistor(d) => d.setup(circuit),
            DeviceEnum::VoltageSource(d) => d.setup(circuit),
            DeviceEnum::CurrentSource(d) => d.setup(circuit),
            DeviceEnum::Capacitor(d) => d.setup(circuit),
            DeviceEnum::Inductor(d) => d.setup(circuit),
            DeviceEnum::Diode(d) => d.setup(circuit),
            DeviceEnum::Bjt(d) => d.setup(circuit),
            DeviceEnum::Mos1(d) => d.setup(circuit),
        }
    }

    fn load(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        match self {
            DeviceEnum::Resistor(d) => d.load(circuit),
            DeviceEnum::VoltageSource(d) => d.load(circuit),
            DeviceEnum::CurrentSource(d) => d.load(circuit),
            DeviceEnum::Capacitor(d) => d.load(circuit),
            DeviceEnum::Inductor(d) => d.load(circuit),
            DeviceEnum::Diode(d) => d.load(circuit),
            DeviceEnum::Bjt(d) => d.load(circuit),
            DeviceEnum::Mos1(d) => d.load(circuit),
        }
    }
}
