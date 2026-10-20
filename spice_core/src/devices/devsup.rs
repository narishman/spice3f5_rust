/// `pnjlim` dynamically limits the per-iteration change of PN junction voltages
/// to avoid exponential overflow inside Newton-Raphson solvers.
pub fn dev_pnjlim(vnew: f64, vold: f64, vt: f64, vcrit: f64, icheck: &mut bool) -> f64 {
    let mut vnew_lim = vnew;

    if (vnew_lim > vcrit) && ((vnew_lim - vold).abs() > (vt + vt)) {
        if vold > 0.0 {
            let arg = 1.0 + (vnew_lim - vold) / vt;
            if arg > 0.0 {
                vnew_lim = vold + vt * arg.ln();
            } else {
                vnew_lim = vcrit;
            }
        } else {
            vnew_lim = vt * (vnew_lim / vt).ln();
        }
        *icheck = true;
    } else {
        *icheck = false;
    }

    vnew_lim
}

/// Limit the per-iteration change of VDS for FETs
pub fn dev_limvds(mut vnew: f64, vold: f64) -> f64 {
    if vold >= 3.5 {
        if vnew > vold {
            vnew = vnew.min(3.0 * vold + 2.0);
        } else {
            if vnew < 3.5 {
                vnew = vnew.max(2.0);
            }
        }
    } else {
        if vnew > vold {
            vnew = vnew.min(4.0);
        } else {
            vnew = vnew.max(-0.5);
        }
    }
    vnew
}

/// Limit the per-iteration change of FET voltages
pub fn dev_fetlim(mut vnew: f64, vold: f64, vto: f64) -> f64 {
    let vtsthi = (2.0 * (vold - vto)).abs() + 2.0;
    let vtstlo = vtsthi / 2.0 + 2.0;
    let vtox = vto + 3.5;
    let delv = vnew - vold;

    if vold >= vto {
        if vold >= vtox {
            if delv <= 0.0 {
                // going off
                if vnew >= vtox {
                    if -delv > vtstlo {
                        vnew = vold - vtstlo;
                    }
                } else {
                    vnew = vnew.max(vto + 2.0);
                }
            } else {
                // staying on
                if delv >= vtsthi {
                    vnew = vold + vtsthi;
                }
            }
        } else {
            // middle region
            if delv <= 0.0 {
                vnew = vnew.max(vto - 0.5);
            } else {
                vnew = vnew.min(vto + 4.0);
            }
        }
    } else {
        // off
        if delv <= 0.0 {
            if -delv > vtsthi {
                vnew = vold - vtsthi;
            }
        } else {
            let vtemp = vto + 0.5;
            if vnew <= vtemp {
                if delv > vtstlo {
                    vnew = vold + vtstlo;
                }
            } else {
                vnew = vtemp;
            }
        }
    }
    vnew
}

/// Compute the MOS overlap capacitances (1/2 of non-constant portion) as functions
/// of the device terminal voltages (1:1 port of `DEVqmeyer` in `devsup.c`).
pub fn dev_qmeyer(
    vgs: f64,
    vgd: f64,
    _vgb: f64,
    von: f64,
    vdsat: f64,
    capgs: &mut f64,
    capgd: &mut f64,
    capgb: &mut f64,
    phi: f64,
    cox: f64,
) {
    let vgst = vgs - von;
    if vgst <= -phi {
        *capgb = cox / 2.0;
        *capgs = 0.0;
        *capgd = 0.0;
    } else if vgst <= -phi / 2.0 {
        *capgb = -vgst * cox / (2.0 * phi);
        *capgs = 0.0;
        *capgd = 0.0;
    } else if vgst <= 0.0 {
        *capgb = -vgst * cox / (2.0 * phi);
        *capgs = vgst * cox / (1.5 * phi) + cox / 3.0;
        *capgd = 0.0;
    } else {
        let vds = vgs - vgd;
        if vdsat <= vds {
            *capgs = cox / 3.0;
            *capgd = 0.0;
            *capgb = 0.0;
        } else {
            let vddif = 2.0 * vdsat - vds;
            let vddif1 = vdsat - vds;
            let vddif2 = vddif * vddif;
            *capgd = cox * (1.0 - vdsat * vdsat / vddif2) / 3.0;
            *capgs = cox * (1.0 - vddif1 * vddif1 / vddif2) / 3.0;
            *capgb = 0.0;
        }
    }
}
