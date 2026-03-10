/// Thin C wrapper around jellyfish C++ API for FFI from Rust.
///
/// This file is only compiled when jellyfish is found via pkg-config.
/// See build.rs for the conditional compilation logic.

// TODO Phase 1: Implement when jellyfish headers are available
//
// Expected implementation:
// #include <jellyfish/mer_dna.hpp>
// #include <jellyfish/jellyfish.hpp>
// #include "jf_wrapper.h"
//
// struct JfDatabase {
//     std::unique_ptr<jellyfish::file_header> header;
//     std::unique_ptr<binary_reader> reader;
//     // hash table handle
// };
//
// extern "C" JfDatabase* jf_open(const char* path) { ... }
// extern "C" void jf_close(JfDatabase* db) { ... }
// extern "C" uint64_t jf_query(const JfDatabase* db, const char* kmer) { ... }
// extern "C" uint32_t jf_kmer_length(const JfDatabase* db) { ... }
// extern "C" bool jf_canonical(const char* kmer, uint32_t k, char* out) { ... }
