//! semantic validation of a parsed case: termination criteria and physical ranges.

use super::schema::{BoundaryKind, CaseConfig, InitialCondition};

impl CaseConfig {
    /// return every violated rule; empty means the case is valid.
    pub fn validate(&self) -> Vec<String> {
        let mut problems = Vec::new();

        if self.metadata.name.trim().is_empty() {
            problems.push("metadata.name must not be empty".to_owned());
        }
        if self.metadata.case_revision < 1 {
            problems.push(format!(
                "metadata.case_revision must be >= 1 (got {})",
                self.metadata.case_revision
            ));
        }

        if self.physics.gamma <= 1.0 {
            problems.push(format!(
                "physics.gamma must be > 1 (got {})",
                self.physics.gamma
            ));
        }
        if self.physics.gas_constant <= 0.0 {
            problems.push(format!(
                "physics.gas_constant must be > 0 (got {})",
                self.physics.gas_constant
            ));
        }
        if let Some(reference) = &self.physics.reference {
            check_positive(reference.rho, "physics.reference.rho", &mut problems);
            check_positive(reference.p, "physics.reference.p", &mut problems);
            check_positive(reference.l, "physics.reference.l", &mut problems);
        }
        if let Some(mu) = self.physics.mu {
            check_positive(mu, "physics.mu", &mut problems);
        }
        check_positive(self.physics.pr, "physics.pr", &mut problems);
        if self.physics.equations == super::schema::Equations::Rans2dSa {
            match &self.physics.turbulence {
                None => problems.push(
                    "physics.turbulence is required when equations = \"rans_2d_sa\"".to_owned(),
                ),
                Some(t) => {
                    check_positive(
                        t.nu_tilde_inf,
                        "physics.turbulence.nu_tilde_inf",
                        &mut problems,
                    );
                    check_positive(t.pr_t, "physics.turbulence.pr_t", &mut problems);
                }
            }
            if self.physics.mu.is_none() {
                problems.push("physics.mu is required when equations = \"rans_2d_sa\"".to_owned());
            }
        }

        if self.mesh.nx < 2 {
            problems.push(format!("mesh.nx must be >= 2 (got {})", self.mesh.nx));
        }
        if self.mesh.x1 <= self.mesh.x0 {
            problems.push(format!(
                "mesh.x1 must be > mesh.x0 (got [{}, {}])",
                self.mesh.x0, self.mesh.x1
            ));
        }
        if self.mesh.source == "clustered" {
            match self.mesh.first_cell {
                Some(h) if h <= 0.0 => {
                    problems.push("mesh.first_cell must be > 0".to_owned());
                }
                None => {
                    problems
                        .push("mesh.first_cell is required when source = \"clustered\"".to_owned());
                }
                _ => {}
            }
            if let Some(r) = self.mesh.growth {
                if r <= 1.0 {
                    problems.push(
                        "mesh.growth must be > 1 (use source = \"uniform\" otherwise)".to_owned(),
                    );
                }
            }
            if let Some(mode) = &self.mesh.cluster {
                if !matches!(mode.as_str(), "wall" | "channel" | "top") {
                    problems.push(format!(
                        "mesh.cluster must be \"wall\", \"channel\", or \"top\" (got \"{mode}\")"
                    ));
                }
            }
        }

        // a profile-inflow side needs its profile spec, and a spec without a
        // profile-inflow side would silently do nothing.
        let is_profile = |k: &BoundaryKind| matches!(k, BoundaryKind::ProfileInflow);
        let has_profile_side = is_profile(&self.boundaries.left)
            || is_profile(&self.boundaries.right)
            || self.boundaries.top.as_ref().is_some_and(is_profile)
            || self.boundaries.bottom.as_ref().is_some_and(is_profile);
        if has_profile_side && self.boundaries.inflow_profile.is_none() {
            problems.push(
                "boundaries.inflow_profile is required when a side is \"profile_inflow\""
                    .to_owned(),
            );
        }
        if !has_profile_side && self.boundaries.inflow_profile.is_some() {
            problems.push(
                "boundaries.inflow_profile is set but no side is \"profile_inflow\"".to_owned(),
            );
        }

        match &self.initial_condition {
            InitialCondition::Uniform { rho, p, .. } => {
                check_positive(*rho, "initial_condition.rho", &mut problems);
                check_positive(*p, "initial_condition.p", &mut problems);
            }
            InitialCondition::TwoState { left, right } => {
                check_positive(left.rho, "initial_condition.left.rho", &mut problems);
                check_positive(left.p, "initial_condition.left.p", &mut problems);
                check_positive(right.rho, "initial_condition.right.rho", &mut problems);
                check_positive(right.p, "initial_condition.right.p", &mut problems);
            }
            InitialCondition::Blast {
                radius,
                ambient,
                fireball,
            } => {
                if *radius <= 0.0 {
                    problems.push(format!(
                        "initial_condition.radius must be > 0 (got {radius})"
                    ));
                }
                check_positive(ambient.rho, "initial_condition.ambient.rho", &mut problems);
                check_positive(ambient.p, "initial_condition.ambient.p", &mut problems);
                check_positive(
                    fireball.rho,
                    "initial_condition.fireball.rho",
                    &mut problems,
                );
                check_positive(fireball.p, "initial_condition.fireball.p", &mut problems);
            }
        }

        let cfl = self.numerics.cfl;
        if cfl <= 0.0 || cfl > 1.0 {
            problems.push(format!("numerics.cfl must be in (0, 1] (got {cfl})"));
        }

        if self.time.final_time < 0.0 {
            problems.push(format!(
                "time.final_time must be >= 0 (got {})",
                self.time.final_time
            ));
        }
        if let Some(target) = self.time.residual_target {
            if target <= 0.0 {
                problems.push(format!("time.residual_target must be > 0 (got {target})"));
            }
        }
        let has_termination = self.time.final_time > 0.0
            || self.time.max_steps > 0
            || self.time.residual_target.is_some();
        if !has_termination {
            problems.push(
                "time needs a termination criterion: final_time > 0, max_steps > 0, or residual_target"
                    .to_owned(),
            );
        }

        if self.output.interval_steps < 1 {
            problems.push(format!(
                "output.interval_steps must be >= 1 (got {})",
                self.output.interval_steps
            ));
        }
        if self.output.formats.is_empty() {
            problems.push("output.formats must list at least one format".to_owned());
        }
        if self.output.fields.is_empty() {
            problems.push("output.fields must list at least one field".to_owned());
        }
        for field in &self.output.fields {
            if field.trim().is_empty() {
                problems.push("output.fields entries must not be empty".to_owned());
            }
        }

        problems
    }
}

fn check_positive(value: f64, context: &str, problems: &mut Vec<String>) {
    if value <= 0.0 {
        problems.push(format!("{context} must be > 0 (got {value})"));
    }
}

#[cfg(test)]
mod tests {
    use super::super::parse_case_toml;
    use super::CaseConfig;

    fn valid_case() -> &'static str {
        r#"
schema_version = 1

[metadata]
name = "sod-shock-tube"
case_revision = 1

[physics]
equations = "euler_1d"
gamma = 1.4
gas_constant = 287.0

[mesh]
nx = 100
x0 = 0.0
x1 = 1.0

[initial_condition]
type = "two_state"
left = { rho = 1.0, u = 0.0, p = 1.0 }
right = { rho = 0.125, u = 0.0, p = 0.1 }

[boundaries]
left = "wall"
right = "wall"

[numerics]
flux = "hll"
reconstruction = "first_order"
cfl = 0.5

[time]
final_time = 0.2
max_steps = 10000

[output]
interval_steps = 100
formats = ["csv", "vtk"]
fields = ["rho", "u", "p"]
"#
    }

    fn problems_for(toml: &str) -> Vec<String> {
        parse_case_toml(toml)
            .err()
            .map(|err| match err {
                super::super::ConfigError::Invalid { problems } => problems,
                other => vec![other.to_string()],
            })
            .unwrap_or_default()
    }

    #[test]
    fn valid_case_passes() {
        let config = parse_case_toml(valid_case()).expect("valid case should parse");
        assert!(CaseConfig::validate(&config).is_empty());
    }

    #[test]
    fn termination_by_residual_alone_is_enough() {
        let toml = valid_case().replace("[time]\nfinal_time = 0.2\n", "[time]\n");
        assert!(problems_for(&toml).is_empty(), "max_steps across the check");
    }

    #[test]
    fn rejects_missing_termination() {
        let toml =
            valid_case().replace("[time]\nfinal_time = 0.2\nmax_steps = 10000\n", "[time]\n");
        let problems = problems_for(&toml);
        assert!(
            problems.iter().any(|p| p.contains("termination criterion")),
            "expected a termination error, got {problems:?}"
        );
    }

    #[test]
    fn rejects_bad_density() {
        let toml = valid_case().replace("right = { rho = 0.125", "right = { rho = -0.125");
        let problems = problems_for(&toml);
        assert!(problems.iter().any(|p| p.contains("right.rho")));
    }

    #[test]
    fn rejects_zero_pressure() {
        let toml = valid_case().replace(
            "left = { rho = 1.0, u = 0.0, p = 1.0 }",
            "left = { rho = 1.0, u = 0.0, p = 0.0 }",
        );
        let problems = problems_for(&toml);
        assert!(problems.iter().any(|p| p.contains("left.p")));
    }

    #[test]
    fn rejects_gamma_not_above_one() {
        for gamma in ["1.0", "0.9"] {
            let toml = valid_case().replace("gamma = 1.4", &format!("gamma = {gamma}"));
            let problems = problems_for(&toml);
            assert!(
                problems.iter().any(|p| p.contains("gamma")),
                "gamma = {gamma} should fail, got {problems:?}"
            );
        }
    }

    #[test]
    fn rejects_gas_constant_nonpositive() {
        let toml = valid_case().replace("gas_constant = 287.0", "gas_constant = 0.0");
        let problems = problems_for(&toml);
        assert!(problems.iter().any(|p| p.contains("gas_constant")));
    }

    #[test]
    fn rejects_mesh_with_fewer_than_two_cells() {
        let toml = valid_case().replace("nx = 100", "nx = 1");
        let problems = problems_for(&toml);
        assert!(problems.iter().any(|p| p.contains("nx")));
    }

    #[test]
    fn rejects_inverted_or_flat_domain() {
        for (x0, x1) in [("1.0", "1.0"), ("1.0", "0.5")] {
            let toml = valid_case()
                .replace("x0 = 0.0", &format!("x0 = {x0}"))
                .replace("x1 = 1.0", &format!("x1 = {x1}"));
            let problems = problems_for(&toml);
            assert!(
                problems.iter().any(|p| p.contains("mesh.x1")),
                "domain [{x0}, {x1}] should fail, got {problems:?}"
            );
        }
    }

    #[test]
    fn rejects_cfl_out_of_range() {
        for cfl in ["0.0", "-0.1", "1.5"] {
            let toml = valid_case().replace("cfl = 0.5", &format!("cfl = {cfl}"));
            let problems = problems_for(&toml);
            assert!(
                problems.iter().any(|p| p.contains("cfl")),
                "cfl = {cfl} should fail, got {problems:?}"
            );
        }
    }

    #[test]
    fn rejects_negative_final_time() {
        let toml = valid_case().replace("final_time = 0.2", "final_time = -1.0");
        let problems = problems_for(&toml);
        assert!(problems.iter().any(|p| p.contains("final_time")));
    }

    #[test]
    fn rejects_negative_residual_target() {
        let toml = valid_case().replace(
            "max_steps = 10000\n",
            "max_steps = 10000\nresidual_target = -1.0\n",
        );
        let problems = problems_for(&toml);
        assert!(problems.iter().any(|p| p.contains("residual_target")));
    }

    #[test]
    fn rejects_zero_output_interval() {
        let toml = valid_case().replace("interval_steps = 100", "interval_steps = 0");
        let problems = problems_for(&toml);
        assert!(problems.iter().any(|p| p.contains("interval_steps")));
    }

    #[test]
    fn rejects_empty_formats_or_fields() {
        let empty_formats = valid_case().replace("formats = [\"csv\", \"vtk\"]", "formats = []");
        assert!(problems_for(&empty_formats)
            .iter()
            .any(|p| p.contains("formats")));
        let empty_fields = valid_case().replace("fields = [\"rho\", \"u\", \"p\"]", "fields = []");
        assert!(problems_for(&empty_fields)
            .iter()
            .any(|p| p.contains("fields")));
    }

    #[test]
    fn rejects_blank_name() {
        let toml = valid_case().replace("name = \"sod-shock-tube\"", "name = \"  \"");
        let problems = problems_for(&toml);
        assert!(problems.iter().any(|p| p.contains("name")));
    }

    #[test]
    fn defaults_apply_when_omitted() {
        let toml = valid_case()
            .replace("unit_system = \"si\"\n", "")
            .replace("case_revision = 1\n", "");
        let config = parse_case_toml(&toml).expect("defaults should fill omitted fields");
        assert_eq!(config.metadata.case_revision, 1);
        assert_eq!(config.physics.unit_system, super::super::UnitSystem::Si);
    }

    fn sa_case() -> String {
        valid_case()
            .replace("equations = \"euler_1d\"", "equations = \"rans_2d_sa\"")
            .replace("gamma = 1.4\n", "gamma = 1.4\nmu = 0.05\npr = 0.72\n")
            .replace("[mesh]\nnx = 100\n", "[mesh]\nnx = 32\nny = 32\n")
            + "\n[physics.turbulence]\nnu_tilde_inf = 0.15\npr_t = 0.9\n"
    }

    #[test]
    fn rans_case_parses_and_validates() {
        let config = parse_case_toml(&sa_case()).expect("rans case should parse");
        assert_eq!(config.physics.equations, super::super::Equations::Rans2dSa);
        assert!((config.physics.mu.unwrap() - 0.05).abs() < 1e-12);
        assert!((config.physics.turbulence.unwrap().nu_tilde_inf - 0.15).abs() < 1e-12);
        assert!(config.validate().is_empty(), "{:?}", config.validate());
    }

    #[test]
    fn rans_case_requires_turbulence_block_and_mu() {
        let no_turb = sa_case().replace(
            "[physics.turbulence]\nnu_tilde_inf = 0.15\npr_t = 0.9\n",
            "",
        );
        let problems = problems_for(&no_turb);
        assert!(
            problems.iter().any(|p| p.contains("physics.turbulence")),
            "got {problems:?}"
        );
        let no_mu = sa_case().replace("mu = 0.05\n", "");
        let problems = problems_for(&no_mu);
        assert!(
            problems.iter().any(|p| p.contains("physics.mu")),
            "got {problems:?}"
        );
    }

    #[test]
    fn boundary_sides_parse_for_2d() {
        let with_sides = valid_case().replace(
            "[boundaries]\nleft = \"wall\"\nright = \"wall\"",
            "[boundaries]\nleft = \"wall\"\nright = \"wall\"\ntop = \"outflow\"\nbottom = \"outflow\"",
        );
        let config = parse_case_toml(&with_sides).expect("2d sides should parse");
        assert_eq!(
            config.boundaries.bottom,
            Some(super::super::BoundaryKind::Outflow)
        );
    }
}
