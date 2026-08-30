//! Method descriptions and the method library.
//!
//! Every integration method in this project is data, not code: a JSON file
//! carries the coefficients, the claimed properties and the reference to the
//! paper it was published in. The code paths are generic over that data.

pub mod coeff_serde;
pub mod lmm;
pub mod rk;
pub mod rosenbrock;

pub use coeff_serde::{CoeffValue, Slot};
pub use lmm::{LmmCoefficients, LmmFamily, LmmFile, Normalization};
pub use rk::{RkRuntime, RkTableau, RkTableauFile, Structure};
pub use rosenbrock::{RosenbrockFile, RosenbrockRuntime, RosenbrockTableau};

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Broad category a method belongs to. Decides which stepper runs it.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum MethodClass {
    RungeKutta,
    LinearMultistep,
    Rosenbrock,
}

/// Bibliographic reference including the DOI of the original publication.
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct Reference {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub authors: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub title: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub year: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub doi: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
}

impl Reference {
    /// Resolvable link, preferring the DOI.
    pub fn link(&self) -> Option<String> {
        if let Some(doi) = &self.doi {
            return Some(format!("https://doi.org/{}", doi.trim_start_matches("https://doi.org/")));
        }
        self.url.clone()
    }
}

/// Properties a method file claims. The analysis verifies them and reports
/// disagreements rather than trusting either side.
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct ClaimedProperties {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub a_stable: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub l_stable: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stiffly_accurate: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub symplectic: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub symmetric: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stage_order: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ssp_coefficient: Option<CoeffValue>,
}

/// Suggested defaults a method ships with.
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct MethodDefaults {
    /// Error controller preset name, see `control::ControllerPreset`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub controller: Option<String>,
    /// Nonlinear solver name, see `nonlinear::SolverKind`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nonlinear_solver: Option<String>,
}

/// A method file as it is stored on disk.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct MethodFile {
    pub id: String,
    pub name: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub aliases: Vec<String>,
    pub class: MethodClass,
    /// Sub family used for grouping in the UI, e.g. `esdirk`, `bdf`, `ssprk`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub family: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub tags: Vec<String>,
    /// Order claimed by the publication.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub order: Option<u32>,
    /// Order of the embedded solution used for error control.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub embedded_order: Option<u32>,
    #[serde(default)]
    pub properties: ClaimedProperties,
    #[serde(default)]
    pub defaults: MethodDefaults,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub references: Vec<Reference>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tableau: Option<RkTableauFile>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub multistep: Option<LmmFile>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rosenbrock: Option<RosenbrockFile>,
}

/// The coefficient carrying part of a method.
#[derive(Clone, Debug)]
pub enum MethodKind {
    RungeKutta(RkTableau),
    LinearMultistep(LmmFamily),
    Rosenbrock(RosenbrockTableau),
}

/// A validated method: metadata plus coefficients.
#[derive(Clone, Debug)]
pub struct Method {
    pub id: String,
    pub name: String,
    pub aliases: Vec<String>,
    pub family: String,
    pub description: Option<String>,
    pub tags: Vec<String>,
    pub declared_order: Option<u32>,
    pub declared_embedded_order: Option<u32>,
    pub properties: ClaimedProperties,
    pub defaults: MethodDefaults,
    pub references: Vec<Reference>,
    pub kind: MethodKind,
}

#[derive(Debug)]
pub enum MethodError {
    Parse(String),
    Invalid(String),
}

impl std::fmt::Display for MethodError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MethodError::Parse(m) => write!(f, "parse error: {m}"),
            MethodError::Invalid(m) => write!(f, "invalid method: {m}"),
        }
    }
}

impl std::error::Error for MethodError {}

impl Method {
    pub fn from_file(file: MethodFile) -> Result<Method, MethodError> {
        let id = file.id.clone();
        let invalid = move |e: String| MethodError::Invalid(format!("{id}: {e}"));
        let blocks = (
            file.tableau.is_some(),
            file.multistep.is_some(),
            file.rosenbrock.is_some(),
        );
        let kind = match (file.class, blocks) {
            (MethodClass::RungeKutta, (true, false, false)) => MethodKind::RungeKutta(
                RkTableau::from_file(file.tableau.as_ref().unwrap())
                    .map_err(|e| invalid(e.to_string()))?,
            ),
            (MethodClass::LinearMultistep, (false, true, false)) => MethodKind::LinearMultistep(
                LmmFamily::from_file(file.multistep.as_ref().unwrap())
                    .map_err(|e| invalid(e.to_string()))?,
            ),
            (MethodClass::Rosenbrock, (false, false, true)) => MethodKind::Rosenbrock(
                RosenbrockTableau::from_file(file.rosenbrock.as_ref().unwrap())
                    .map_err(|e| invalid(e.to_string()))?,
            ),
            _ => return Err(invalid("the class and the coefficient block do not match".into())),
        };

        Ok(Method {
            id: file.id,
            name: file.name,
            aliases: file.aliases,
            family: file.family.unwrap_or_else(|| "other".to_string()),
            description: file.description,
            tags: file.tags,
            declared_order: file.order,
            declared_embedded_order: file.embedded_order,
            properties: file.properties,
            defaults: file.defaults,
            references: file.references,
            kind,
        })
    }

    pub fn parse(json: &str) -> Result<Method, MethodError> {
        let file: MethodFile = serde_json::from_str(json).map_err(|e| MethodError::Parse(e.to_string()))?;
        Method::from_file(file)
    }

    pub fn class(&self) -> MethodClass {
        match self.kind {
            MethodKind::RungeKutta(_) => MethodClass::RungeKutta,
            MethodKind::LinearMultistep(_) => MethodClass::LinearMultistep,
            MethodKind::Rosenbrock(_) => MethodClass::Rosenbrock,
        }
    }

    pub fn tableau(&self) -> Option<&RkTableau> {
        match &self.kind {
            MethodKind::RungeKutta(t) => Some(t),
            _ => None,
        }
    }

    pub fn multistep(&self) -> Option<&LmmFamily> {
        match &self.kind {
            MethodKind::LinearMultistep(m) => Some(m),
            _ => None,
        }
    }

    pub fn rosenbrock(&self) -> Option<&RosenbrockTableau> {
        match &self.kind {
            MethodKind::Rosenbrock(r) => Some(r),
            _ => None,
        }
    }

    /// Number of stages for Runge-Kutta and Rosenbrock, number of steps for
    /// multistep.
    pub fn size(&self) -> usize {
        match &self.kind {
            MethodKind::RungeKutta(t) => t.stages,
            MethodKind::LinearMultistep(m) => m.steps,
            MethodKind::Rosenbrock(r) => r.stages,
        }
    }

    pub fn is_implicit(&self) -> bool {
        match &self.kind {
            MethodKind::RungeKutta(t) => !t.is_explicit(),
            MethodKind::LinearMultistep(m) => m.implicit,
            // A Rosenbrock method solves a linear system at every stage. It
            // never iterates, but it is not explicit either.
            MethodKind::Rosenbrock(_) => true,
        }
    }

    pub fn is_adaptive(&self) -> bool {
        match &self.kind {
            MethodKind::RungeKutta(t) => t.has_embedded(),
            // Multistep error estimates come from the order reduced formula,
            // which only exists while the family still has a shorter member.
            MethodKind::LinearMultistep(m) => m.steps > m.min_steps,
            MethodKind::Rosenbrock(r) => r.has_embedded(),
        }
    }
}

/// A collection of methods, addressable by id or alias.
#[derive(Clone, Debug, Default)]
pub struct MethodLibrary {
    methods: Vec<Method>,
    index: HashMap<String, usize>,
}

impl MethodLibrary {
    pub fn new() -> Self {
        MethodLibrary::default()
    }

    pub fn insert(&mut self, method: Method) -> Result<(), MethodError> {
        let position = self.methods.len();
        let mut keys = vec![method.id.to_ascii_lowercase()];
        keys.extend(method.aliases.iter().map(|a| a.to_ascii_lowercase()));
        for key in &keys {
            if let Some(&existing) = self.index.get(key) {
                return Err(MethodError::Invalid(format!(
                    "duplicate method key {key:?}, already used by {}",
                    self.methods[existing].id
                )));
            }
        }
        for key in keys {
            self.index.insert(key, position);
        }
        self.methods.push(method);
        Ok(())
    }

    pub fn get(&self, key: &str) -> Option<&Method> {
        self.index
            .get(&key.to_ascii_lowercase())
            .map(|&i| &self.methods[i])
    }

    pub fn len(&self) -> usize {
        self.methods.len()
    }

    pub fn is_empty(&self) -> bool {
        self.methods.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Method> {
        self.methods.iter()
    }

    pub fn ids(&self) -> Vec<&str> {
        self.methods.iter().map(|m| m.id.as_str()).collect()
    }

    /// Load every `*.json` below a directory tree.
    pub fn from_directory(root: impl AsRef<std::path::Path>) -> Result<MethodLibrary, MethodError> {
        let mut library = MethodLibrary::new();
        let mut files = Vec::new();
        collect_json(root.as_ref(), &mut files)
            .map_err(|e| MethodError::Parse(format!("cannot read method directory: {e}")))?;
        files.sort();
        for path in files {
            let text = std::fs::read_to_string(&path)
                .map_err(|e| MethodError::Parse(format!("{}: {e}", path.display())))?;
            let method = Method::parse(&text)
                .map_err(|e| MethodError::Parse(format!("{}: {e}", path.display())))?;
            library.insert(method)?;
        }
        Ok(library)
    }

    /// The method library baked into the binary at compile time.
    #[cfg(feature = "embedded-methods")]
    pub fn embedded() -> Result<MethodLibrary, MethodError> {
        let mut library = MethodLibrary::new();
        for (path, text) in crate::embedded::EMBEDDED_METHODS {
            let method = Method::parse(text).map_err(|e| MethodError::Parse(format!("{path}: {e}")))?;
            library.insert(method)?;
        }
        Ok(library)
    }
}

fn collect_json(dir: &std::path::Path, out: &mut Vec<std::path::PathBuf>) -> std::io::Result<()> {
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            collect_json(&path, out)?;
        } else if path.extension().map_or(false, |e| e == "json") {
            out.push(path);
        }
    }
    Ok(())
}
