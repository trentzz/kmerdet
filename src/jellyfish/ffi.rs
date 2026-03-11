/// FFI bindings to the jellyfish C++ library via jf_wrapper.
///
/// Only available when compiled with jellyfish (cfg `has_jellyfish`).
/// See `jf_wrapper/jf_wrapper.h` for the C API.

#[cfg(has_jellyfish)]
mod bindings {
    use std::ffi::CString;
    use std::path::Path;

    extern "C" {
        fn jf_open(path: *const std::ffi::c_char) -> *mut std::ffi::c_void;
        fn jf_close(db: *mut std::ffi::c_void);
        fn jf_query(db: *const std::ffi::c_void, kmer: *const std::ffi::c_char) -> u64;
        fn jf_kmer_length(db: *const std::ffi::c_void) -> u32;
    }

    pub struct JellyfishDb {
        handle: *mut std::ffi::c_void,
        k: u8,
    }

    // Safety: JellyfishDb is read-only after construction.
    unsafe impl Send for JellyfishDb {}
    unsafe impl Sync for JellyfishDb {}

    impl JellyfishDb {
        pub fn open(path: &Path) -> anyhow::Result<Self> {
            let c_path =
                CString::new(path.to_str().ok_or_else(|| anyhow::anyhow!("invalid path"))?)
                    .map_err(|_| anyhow::anyhow!("path contains null byte"))?;

            let handle = unsafe { jf_open(c_path.as_ptr()) };
            if handle.is_null() {
                anyhow::bail!("failed to open jellyfish database: {}", path.display());
            }

            let k = unsafe { jf_kmer_length(handle) } as u8;
            Ok(Self { handle, k })
        }
    }

    impl crate::jellyfish::KmerDatabase for JellyfishDb {
        fn query(&self, kmer: &str) -> u64 {
            let Ok(c_kmer) = CString::new(kmer) else {
                return 0;
            };
            unsafe { jf_query(self.handle, c_kmer.as_ptr()) }
        }

        fn kmer_length(&self) -> u8 {
            self.k
        }
    }

    impl Drop for JellyfishDb {
        fn drop(&mut self) {
            if !self.handle.is_null() {
                unsafe { jf_close(self.handle) };
            }
        }
    }
}

#[cfg(has_jellyfish)]
pub use bindings::JellyfishDb;
