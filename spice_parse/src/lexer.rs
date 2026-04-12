use nom::{
    IResult,
    branch::alt,
    bytes::complete::{tag_no_case, take_while1},
    character::complete::{alpha1, alphanumeric1, char, multispace0, space0, space1},
    combinator::{map, opt, recognize},
    number::complete::double,
    sequence::{delimited, pair, tuple},
};

/// Consume leading whitespace except newlines
pub fn sp(i: &str) -> IResult<&str, &str> {
    space0(i)
}

/// Identifiers for nodes and device names
pub fn identifier(input: &str) -> IResult<&str, String> {
    let res: IResult<&str, &str> = recognize(pair(
        take_while1(|c: char| c.is_alphabetic() || c == '_'),
        take_while1(|c: char| c.is_alphanumeric() || c == '_' || c == '.' || c == '-'),
    ))(input);

    let (input, out) = match res {
        Ok(v) => v,
        Err(_) => take_while1(|c: char| c.is_alphanumeric() || c == '_')(input)?,
    };
    Ok((input, out.to_string()))
}

/// Parse SPICE floats with SI Suffixes exactly like INPevaluate
/// Support T, G, MEG, K, MIL, M, U, N, P, F
pub fn spice_val(input: &str) -> IResult<&str, f64> {
    let (input, val) = double(input)?;

    // Look for scale factors
    let ext = alt((
        tag_no_case("meg"),
        tag_no_case("mil"),
        tag_no_case("t"),
        tag_no_case("g"),
        tag_no_case("k"),
        tag_no_case("u"),
        tag_no_case("n"),
        tag_no_case("p"),
        tag_no_case("f"),
        tag_no_case("m"),
    ));

    let (input, suffix) = opt(ext)(input)?;

    let scale = match suffix.map(|s: &str| s.to_lowercase()).as_deref() {
        Some("t") => 1e12,
        Some("g") => 1e9,
        Some("meg") => 1e6,
        Some("k") => 1e3,
        Some("mil") => 25.4e-6,
        Some("m") => 1e-3,
        Some("u") => 1e-6,
        Some("n") => 1e-9,
        Some("p") => 1e-12,
        Some("f") => 1e-15,
        _ => 1.0,
    };

    // SPICE ignores anything trailing the suffix until the next boundary space
    let res: IResult<&str, &str> = take_while1(|c: char| c.is_alphanumeric())(input);
    let (input, _) = match res {
        Ok(v) => v,
        Err(_) => (input, ""),
    };

    Ok((input, val * scale))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spice_val() {
        assert!((spice_val("10k").unwrap().1 - 10000.0).abs() < f64::EPSILON);
        assert!((spice_val("1.5meg").unwrap().1 - 1500000.0).abs() < f64::EPSILON);
        assert!((spice_val("2.2MegOhms").unwrap().1 - 2200000.0).abs() < f64::EPSILON); // Checks trailing gobble
        assert!((spice_val("10u").unwrap().1 - 0.00001).abs() < 1e-15);
        assert!((spice_val("4mil").unwrap().1 - 4.0 * 25.4e-6).abs() < f64::EPSILON);
        assert!((spice_val("-5.5").unwrap().1 - -5.5).abs() < f64::EPSILON);
    }
}
