use nom::{
    IResult,
    branch::alt,
    bytes::complete::tag_no_case,
    character::complete::{not_line_ending, space1},
    combinator::opt,
    sequence::tuple,
};

use crate::ast::{DeviceStmt, Stmt};
use crate::lexer::{identifier, sp, spice_val};

/// Parse a Resistor: RNAME N1 N2 VALUE [PARAMS...]
fn parse_resistor(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, name, _, n1, _, n2, _, val, _, _rest)) = tuple((
        tag_no_case("R"),
        identifier,
        space1,
        identifier,
        space1,
        identifier,
        space1,
        spice_val,
        sp,
        not_line_ending, // Handle extra params later when supported
    ))(input)?;

    Ok((
        input,
        Stmt::Device(DeviceStmt {
            ident: format!("R{}", name),
            nodes: vec![n1, n2],
            value: val,
            params: vec![],
        }),
    ))
}

/// Parse a Voltage Source: VNAME N+ N- [DC] VALUE
fn parse_vsource(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, name, _, n1, _, n2, _, _dc_tag, _, val, _, _rest)) = tuple((
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
            value: val,
            params: vec![],
        }),
    ))
}

/// Parse a Current Source: INAME N+ N- [DC] VALUE
fn parse_isource(input: &str) -> IResult<&str, Stmt> {
    let (input, (_, name, _, n1, _, n2, _, _dc_tag, _, val, _, _rest)) = tuple((
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
            value: val,
            params: vec![],
        }),
    ))
}

pub fn parse_line(input: &str) -> IResult<&str, Stmt> {
    let (input, _) = sp(input)?;

    // Check for comment
    if input.starts_with('*') {
        let (input, comment) = not_line_ending(input)?;
        return Ok((input, Stmt::Comment(comment.to_string())));
    }

    alt((parse_resistor, parse_vsource, parse_isource))(input)
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
}
