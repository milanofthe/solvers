use solvers_core::analysis::{stability_region, suggested_window};
use solvers_core::method::MethodLibrary;
use std::time::Instant;

fn main() {
    let library = MethodLibrary::embedded().unwrap();
    for (id, n) in [
        ("rkdp87", 140usize),
        ("radau_iia_9", 140),
        ("rodas5", 140),
        ("bdf6", 72),
        ("adams_bashforth_8", 72),
        ("adams_moulton_7", 72),
    ] {
        let method = library.get(id).unwrap();
        let start = Instant::now();
        let (re, im) = suggested_window(method);
        let window = start.elapsed();
        let start = Instant::now();
        let grid = stability_region(method, re, im, n, n).unwrap();
        let sample = start.elapsed();
        println!(
            "{id:<20} {n}x{n}  window {:>8.2?}  grid {:>8.2?}  per point {:>6.2?}",
            window,
            sample,
            sample / (grid.magnitude.len() as u32)
        );
    }
}
