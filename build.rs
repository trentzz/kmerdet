fn main() {
    println!("cargo::rustc-check-cfg=cfg(has_jellyfish)");
}
