//! Command line access to the method library and its analyses.

use solvers_core::analysis;
use solvers_core::method::MethodLibrary;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let command = args.first().map(String::as_str).unwrap_or("help");

    let library = match MethodLibrary::embedded() {
        Ok(library) => library,
        Err(error) => {
            eprintln!("method library failed to load: {error}");
            std::process::exit(1);
        }
    };

    match command {
        "list" => {
            let mut rows: Vec<&solvers_core::Method> = library.iter().collect();
            rows.sort_by(|a, b| (a.family.as_str(), a.id.as_str()).cmp(&(&b.family, &b.id)));
            println!("{:<22} {:<14} {:>5} {:>7}  {}", "id", "family", "order", "stages", "name");
            for method in rows {
                println!(
                    "{:<22} {:<14} {:>5} {:>7}  {}",
                    method.id,
                    method.family,
                    method
                        .declared_order
                        .map(|o| o.to_string())
                        .unwrap_or_else(|| "-".into()),
                    method.size(),
                    method.name
                );
            }
            println!("\n{} methods", library.len());
        }
        "analyze" => {
            let Some(id) = args.get(1) else {
                eprintln!("usage: solvers analyze <method-id>");
                std::process::exit(2);
            };
            let Some(method) = library.get(id) else {
                eprintln!("unknown method: {id}");
                std::process::exit(2);
            };
            let report = analysis::analyze(method);
            println!("{}", serde_json::to_string_pretty(&report).unwrap());
        }
        "verify" => {
            let mut failures = 0;
            for method in library.iter() {
                let report = analysis::analyze(method);
                if report.discrepancies.is_empty() {
                    continue;
                }
                failures += 1;
                println!("{}:", method.id);
                for issue in &report.discrepancies {
                    println!("  {issue}");
                }
            }
            if failures == 0 {
                println!("all {} methods agree with their analysis", library.len());
            } else {
                println!("\n{failures} of {} methods disagree", library.len());
                std::process::exit(1);
            }
        }
        _ => {
            println!("usage: solvers <command>");
            println!("  list                 list every method in the library");
            println!("  analyze <method-id>  derive order and stability from the coefficients");
            println!("  verify               check every method file against its analysis");
        }
    }
}
