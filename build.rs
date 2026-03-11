fn main() {
    println!("cargo::rustc-check-cfg=cfg(has_jellyfish)");

    // Attempt to find jellyfish via pkg-config
    match pkg_config::probe_library("jellyfish-2.0") {
        Ok(lib) => {
            let mut build = cc::Build::new();
            build
                .cpp(true)
                .file("jf_wrapper/jf_wrapper.cpp")
                .flag("-std=c++11");

            // Add include paths from pkg-config
            for path in &lib.include_paths {
                build.include(path);
            }

            build.compile("jf_wrapper");

            // Link the jellyfish library
            for path in &lib.link_paths {
                println!("cargo:rustc-link-search=native={}", path.display());
            }
            for lib_name in &lib.libs {
                println!("cargo:rustc-link-lib={}", lib_name);
            }
            // Jellyfish is C++ — need to link libstdc++
            println!("cargo:rustc-link-lib=stdc++");

            println!("cargo:rustc-cfg=has_jellyfish");
        }
        Err(_) => {
            println!("cargo:warning=Jellyfish not found; FFI bindings will not be available.");
            println!("cargo:warning=Install jellyfish 2.2+ and ensure pkg-config can find it.");
        }
    }
}
