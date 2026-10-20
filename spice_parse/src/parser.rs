use nom::{
    IResult,
    branch::alt,
    bytes::complete::tag_no_case,
    character::complete::{not_line_ending, space1},
    combinator::opt,
    sequence::tuple,
};

use crate::ast::{ControlStmt, DeviceStmt, ModelStmt, Stmt};
use crate::lexer::{Lexer, identifier, param_pair, sp, spice_val};

/// High-level SPICE netlist parser consuming logical line tokens from `Lexer`
#[derive(Debug, Clone)]
pub struct Parser {
    lines: Vec<String>,
}

impl Parser {
    pub fn new(lines: Vec<String>) -> Self {
        Self { lines }
    }

    pub fn parse(&mut self) -> Result<Vec<Stmt>, String> {
        let mut stmts = Vec::new();
        for (idx, line) in self.lines.iter().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            match parse_line(trimmed) {
                Ok((_, stmt)) => stmts.push(stmt),
                Err(e) => {
                    return Err(format!(
                        "Parse error on line {}: '{}' ({})",
                        idx + 1,
                        line,
                        e
                    ));
                }
            }
        }
        Ok(stmts)
    }
}

/// Convenience helper to lex and parse a full multi-line SPICE netlist string
pub fn parse_netlist(input: &str) -> Result<Vec<Stmt>, String> {
    let mut lexer = Lexer::new(input);
    let tokens = lexer.tokenize();
    let mut parser = Parser::new(tokens);
    parser.parse()
}

/// Parse trailing `KEY=VAL` parameters from the remainder of a line
fn parse_trailing_params(mut rest: &str) -> Vec<(String, f64)> {
    let mut params = Vec::new();
    loop {
        rest = rest
            .trim_start_matches(|c: char| c.is_whitespace() || c == '(' || c == ')' || c == ',');
        if rest.is_empty() || rest.starts_with('*') || rest.starts_with(';') {
            break;
        }
        if let Ok((next_rest, kv)) = param_pair(rest) {
            params.push(kv);
            rest = next_rest;
        } else {
            // Skip one token if unrecognized
            if let Some(idx) = rest.find(|c: char| c.is_whitespace() || c == ')' || c == ',') {
                rest = &rest[idx..];
            } else {
                break;
            }
        }
    }
    params
}

/// Parse optional model identifier and trailing `KEY=VAL` parameters
fn parse_optional_model_and_params(rest: &str) -> (Option<String>, Vec<(String, f64)>) {
    let trimmed = rest.trim();
    if trimmed.is_empty() || trimmed.starts_with('*') || trimmed.starts_with(';') {
        return (None, Vec::new());
    }
    // Check if the first token is a KEY=VAL pair or a model identifier
    if param_pair(trimmed).is_ok() {
        (None, parse_trailing_params(trimmed))
    } else if let Ok((after_model, model_name)) = identifier(trimmed) {
        (Some(model_name), parse_trailing_params(after_model))
    } else {
        (None, parse_trailing_params(trimmed))
    }
}

/// Parse a Resistor: RNAME N1 N2 VALUE [PARAMS...]
fn parse_resistor(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, name, _, n1, _, n2, _, val, _, rest)) = tuple((
        tag_no_case("R"),
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        space1,
        spice_val,
        sp,
        not_line_ending,
    ))(input)?;

    Ok((
        input,
        Stmt::Device(DeviceStmt {
            ident: format!("R{}", name),
            nodes: vec![n1, n2],
            model: None,
            value: val,
            params: parse_trailing_params(rest),
        }),
    ))
}

/// Parse a Capacitor: CNAME N1 N2 VALUE [PARAMS...]
fn parse_capacitor(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, name, _, n1, _, n2, _, val, _, rest)) = tuple((
        tag_no_case("C"),
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        space1,
        spice_val,
        sp,
        not_line_ending,
    ))(input)?;

    Ok((
        input,
        Stmt::Device(DeviceStmt {
            ident: format!("C{}", name),
            nodes: vec![n1, n2],
            model: None,
            value: val,
            params: parse_trailing_params(rest),
        }),
    ))
}

/// Parse an Inductor: LNAME N1 N2 VALUE [PARAMS...]
fn parse_inductor(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, name, _, n1, _, n2, _, val, _, rest)) = tuple((
        tag_no_case("L"),
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        space1,
        spice_val,
        sp,
        not_line_ending,
    ))(input)?;

    Ok((
        input,
        Stmt::Device(DeviceStmt {
            ident: format!("L{}", name),
            nodes: vec![n1, n2],
            model: None,
            value: val,
            params: parse_trailing_params(rest),
        }),
    ))
}

/// Parse a Voltage Source: VNAME N+ N- [DC] VALUE
fn parse_vsource(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, name, _, n1, _, n2, _, _dc_tag, _, val, _, rest)) = tuple((
        tag_no_case("V"),
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        sp,
        opt(tag_no_case("DC")),
        sp,
        spice_val,
        sp,
        not_line_ending,
    ))(input)?;

    Ok((
        input,
        Stmt::Device(DeviceStmt {
            ident: format!("V{}", name),
            nodes: vec![n1, n2],
            model: None,
            value: val,
            params: parse_trailing_params(rest),
        }),
    ))
}

/// Parse a Current Source: INAME N+ N- [DC] VALUE
fn parse_isource(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, name, _, n1, _, n2, _, _dc_tag, _, val, _, rest)) = tuple((
        tag_no_case("I"),
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        sp,
        opt(tag_no_case("DC")),
        sp,
        spice_val,
        sp,
        not_line_ending,
    ))(input)?;

    Ok((
        input,
        Stmt::Device(DeviceStmt {
            ident: format!("I{}", name),
            nodes: vec![n1, n2],
            model: None,
            value: val,
            params: parse_trailing_params(rest),
        }),
    ))
}

/// Parse a Diode: DNAME N+ N- [MODEL] [PARAMS...]
fn parse_diode(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, name, _, n1, _, n2, _, rest)) = tuple((
        tag_no_case("D"),
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        sp,
        not_line_ending,
    ))(input)?;

    let (model, params) = parse_optional_model_and_params(rest);
    Ok((
        input,
        Stmt::Device(DeviceStmt {
            ident: format!("D{}", name),
            nodes: vec![n1, n2],
            model,
            value: 0.0,
            params,
        }),
    ))
}

/// Parse a BJT: QNAME NC NB NE [NS] [MODEL] [PARAMS...]
fn parse_bjt(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, name, _, nc, _, nb, _, ne, _, rest)) = tuple((
        tag_no_case("Q"),
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        sp,
        not_line_ending,
    ))(input)?;

    let (model, params) = parse_optional_model_and_params(rest);
    Ok((
        input,
        Stmt::Device(DeviceStmt {
            ident: format!("Q{}", name),
            nodes: vec![nc, nb, ne],
            model,
            value: 0.0,
            params,
        }),
    ))
}

/// Parse a MOSFET: MNAME ND NG NS NB [MODEL] [PARAMS...]
fn parse_mosfet(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, name, _, nd, _, ng, _, ns, _, nb, _, rest)) = tuple((
        tag_no_case("M"),
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        sp,
        not_line_ending,
    ))(input)?;

    let (model, params) = parse_optional_model_and_params(rest);
    Ok((
        input,
        Stmt::Device(DeviceStmt {
            ident: format!("M{}", name),
            nodes: vec![nd, ng, ns, nb],
            model,
            value: 0.0,
            params,
        }),
    ))
}

/// Parse a `.MODEL` card: .MODEL MNAME TYPE [(] [KEY=VAL ...] [)]
fn parse_model(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, _, name, _, model_type, _, rest)) = tuple((
        tag_no_case(".MODEL"),
        space1,
        identifier,
        space1,
        identifier,
        sp,
        not_line_ending,
    ))(input)?;

    let clean_type = model_type
        .trim_matches(|c: char| c == '(' || c == ')')
        .to_uppercase();
    let params = parse_trailing_params(rest);

    Ok((
        input,
        Stmt::Model(ModelStmt {
            name,
            model_type: clean_type,
            params,
        }),
    ))
}

/// Parse `.DC`, `.TRAN`, `.OP`, `.END` control statements
fn parse_control(input: &str) -> IResult<&str, Stmt> {
    if let Ok((rem, (_, _, src, _, start, _, stop, _, step, _, _))) = tuple((
        tag_no_case(".DC"),
        space1,
        identifier,
        space1,
        spice_val,
        space1,
        spice_val,
        space1,
        spice_val,
        sp,
        not_line_ending,
    ))(input)
    {
        return Ok((
            rem,
            Stmt::Control(ControlStmt::Dc {
                source: src,
                start,
                stop,
                step,
            }),
        ));
    }

    if let Ok((rem, (_, _, step, _, stop, _, _))) = tuple((
        tag_no_case(".TRAN"),
        space1,
        spice_val,
        space1,
        spice_val,
        sp,
        not_line_ending,
    ))(input)
    {
        return Ok((rem, Stmt::Control(ControlStmt::Tran { step, stop })));
    }

    if let Ok((rem, (_, _, _))) = tuple((tag_no_case(".OP"), sp, not_line_ending))(input) {
        return Ok((rem, Stmt::Control(ControlStmt::Op)));
    }

    let (rem, (_, _, _)) = tuple((tag_no_case(".END"), sp, not_line_ending))(input)?;
    Ok((rem, Stmt::Control(ControlStmt::End)))
}

pub fn parse_line(input: &str) -> IResult<&str, Stmt> {
    let (input, _) = sp(input)?;

    // Check for comment
    if input.starts_with('*') {
        let (input, comment) = not_line_ending(input)?;
        return Ok((input, Stmt::Comment(comment.to_string())));
    }

    alt((
        parse_model,
        parse_control,
        parse_resistor,
        parse_capacitor,
        parse_inductor,
        parse_vsource,
        parse_isource,
        parse_diode,
        parse_bjt,
        parse_mosfet,
    ))(input)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_resistor() {
        let (rem, stmt) = parse_line("R1 in out 10k").unwrap();
        assert_eq!(rem, "");
        if let Stmt::Device(d) = stmt {
            assert_eq!(d.ident, "R1");
            assert_eq!(d.nodes, vec!["in", "out"]);
            assert_eq!(d.value, 10000.0);
        } else {
            panic!("Expected device");
        }
    }

    #[test]
    fn test_parse_sources() {
        let (_, stmt) = parse_line("V_in 1 0 DC 5.0").unwrap();
        if let Stmt::Device(d) = stmt {
            assert_eq!(d.ident, "V_in");
            assert_eq!(d.nodes, vec!["1", "0"]);
            assert_eq!(d.value, 5.0);
        }

        let (_, stmt) = parse_line("I1 2 0 1m").unwrap();
        if let Stmt::Device(d) = stmt {
            assert_eq!(d.ident, "I1");
            assert_eq!(d.nodes, vec!["2", "0"]);
            assert_eq!(d.value, 0.001);
        }
    }

    #[test]
    fn test_parse_active_and_passive_devices() {
        let (_, stmt) = parse_line("C1 out 0 1u").unwrap();
        if let Stmt::Device(d) = stmt {
            assert_eq!(d.ident, "C1");
            assert_eq!(d.nodes, vec!["out", "0"]);
            assert!((d.value - 1e-6).abs() < 1e-15);
        } else {
            panic!("Expected Capacitor");
        }

        let (_, stmt) = parse_line("L1 1 2 10m").unwrap();
        if let Stmt::Device(d) = stmt {
            assert_eq!(d.ident, "L1");
            assert_eq!(d.nodes, vec!["1", "2"]);
            assert!((d.value - 1e-2).abs() < 1e-12);
        } else {
            panic!("Expected Inductor");
        }

        let (_, stmt) = parse_line("D1 anode cathode DMOD").unwrap();
        if let Stmt::Device(d) = stmt {
            assert_eq!(d.ident, "D1");
            assert_eq!(d.nodes, vec!["anode", "cathode"]);
            assert_eq!(d.model.as_deref(), Some("DMOD"));
        } else {
            panic!("Expected Diode");
        }

        let (_, stmt) = parse_line("Q1 c b e Q2N2222").unwrap();
        if let Stmt::Device(d) = stmt {
            assert_eq!(d.ident, "Q1");
            assert_eq!(d.nodes, vec!["c", "b", "e"]);
            assert_eq!(d.model.as_deref(), Some("Q2N2222"));
        } else {
            panic!("Expected BJT");
        }

        let (_, stmt) = parse_line("M1 d g s b NMOS_MOD W=20u L=2u").unwrap();
        if let Stmt::Device(d) = stmt {
            assert_eq!(d.ident, "M1");
            assert_eq!(d.nodes, vec!["d", "g", "s", "b"]);
            assert_eq!(d.model.as_deref(), Some("NMOS_MOD"));
            assert_eq!(d.params.len(), 2);
        } else {
            panic!("Expected MOSFET");
        }
    }

    #[test]
    fn test_parse_model_and_multiline() {
        let netlist = "
        * Sample Netlist
        V1 1 0 5.0
        R1 1 2 1k
        Q1 2 1 0 QMOD
        .MODEL QMOD NPN (BF=150 IS=1e-15)
        .END
        ";
        let mut lexer = Lexer::new(netlist);
        let tokens = lexer.tokenize();
        let mut parser = Parser::new(tokens);
        let stmts = parser.parse().unwrap();
        assert_eq!(stmts.len(), 6);
    }
}
