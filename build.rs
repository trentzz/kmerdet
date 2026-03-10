fn main() {
    println!("cargo::rustc-check-cfg=cfg(has_jellyfish)");

    // Attempt to find jellyfish via pkg-config
    let jf_found = pkg_config::probe_library("jellyfish-2.0").is_ok();

    if jf_found {
        cc::Build::new()
            .cpp(true)
            .file("jf_wrapper/jf_wrapper.cpp")
            .flag("-std=c++11")
            .compile("jf_wrapper");
        println!("cargo:rustc-cfg=has_jellyfish");
    } else {
        println!("cargo:warning=Jellyfish not found; FFI bindings will not be available.");
        println!("cargo:warning=Install jellyfish 2.2+ and ensure pkg-config can find it.");
    }
}
