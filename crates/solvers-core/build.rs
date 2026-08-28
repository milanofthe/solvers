//! Bakes the JSON method library into the crate.
//!
//! The files under `methods/` stay the single source of truth; this only
//! generates the `include_str!` table so a binary or a wasm bundle does not
//! need the directory at runtime.

use std::path::{Path, PathBuf};

fn main() {
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let methods = manifest.join("../../methods");
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("embedded_methods.rs");

    let mut files = Vec::new();
    if methods.is_dir() {
        collect(&methods, &mut files);
    }
    files.sort();

    println!("cargo:rerun-if-changed={}", methods.display());
    for file in &files {
        println!("cargo:rerun-if-changed={}", file.display());
    }

    let mut source = String::from(
        "/// Method files embedded at compile time as (relative path, contents).\n\
         pub static EMBEDDED_METHODS: &[(&str, &str)] = &[\n",
    );
    for file in &files {
        let rel = file
            .strip_prefix(&methods)
            .unwrap_or(file)
            .to_string_lossy()
            .replace('\\', "/");
        let abs = file.to_string_lossy().replace('\\', "/");
        source.push_str(&format!("    ({rel:?}, include_str!({abs:?})),\n"));
    }
    source.push_str("];\n");

    std::fs::write(&out, source).expect("cannot write embedded method table");
}

fn collect(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            collect(&path, out);
        } else if path.extension().map_or(false, |e| e == "json") {
            out.push(path);
        }
    }
}
