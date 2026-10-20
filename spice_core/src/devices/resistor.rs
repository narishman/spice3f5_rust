use super::Device;

use crate::circuit::Circuit;

pub struct Resistor {
    pub name: String,
    pub pos_node: usize,
    pub neg_node: usize,
    pub resist: f64,
    pub conduct: f64,

    pub pos_pos_ptr: Option<usize>,
    pub neg_neg_ptr: Option<usize>,
    pub pos_neg_ptr: Option<usize>,
    pub neg_pos_ptr: Option<usize>,
}

impl Resistor {
    pub fn new(name: impl Into<String>, pos_node: usize, neg_node: usize, resist: f64) -> Self {
        let conduct = if resist != 0.0 { 1.0 / resist } else { 1e12 };
        Self {
            name: name.into(),
            pos_node,
            neg_node,
            resist,
            conduct,
            pos_pos_ptr: None,
            neg_neg_ptr: None,
            pos_neg_ptr: None,
            neg_pos_ptr: None,
        }
    }
}

impl Device for Resistor {
    fn setup(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        self.pos_pos_ptr = circuit.matrix.get_element(self.pos_node, self.pos_node);
        self.neg_neg_ptr = circuit.matrix.get_element(self.neg_node, self.neg_node);
        self.pos_neg_ptr = circuit.matrix.get_element(self.pos_node, self.neg_node);
        self.neg_pos_ptr = circuit.matrix.get_element(self.neg_node, self.pos_node);
        Ok(())
    }

    fn load(&mut self, circuit: &mut Circuit) -> Result<(), String> {
        if let Some(idx) = self.pos_pos_ptr {
            circuit.matrix.elements_arena[idx].real += self.conduct;
        }
        if let Some(idx) = self.neg_neg_ptr {
            circuit.matrix.elements_arena[idx].real += self.conduct;
        }
        if let Some(idx) = self.pos_neg_ptr {
            circuit.matrix.elements_arena[idx].real -= self.conduct;
        }
        if let Some(idx) = self.neg_pos_ptr {
            circuit.matrix.elements_arena[idx].real -= self.conduct;
        }
        Ok(())
    }
}
