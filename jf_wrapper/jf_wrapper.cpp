/// Thin C wrapper around jellyfish C++ API for FFI from Rust.
///
/// This file is only compiled when jellyfish is found via pkg-config.
/// See build.rs for the conditional compilation logic.

#include <jellyfish/jellyfish.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/file_header.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/large_hash_array.hpp>
#include "jf_wrapper.h"

#include <fstream>
#include <memory>
#include <string>
#include <cstring>

// The raw array type for memory-mapped random-access queries.
typedef jellyfish::large_hash::array_raw<jellyfish::mer_dna> array_raw_type;

struct JfDatabase {
    std::unique_ptr<jellyfish::mapped_file> mapped;
    std::unique_ptr<array_raw_type>         hash;
    uint32_t kmer_length;
    bool     canonical;
};

extern "C" {

JfDatabase* jf_open(const char* path) {
    try {
        auto db = new JfDatabase();

        // Read the file header via an input stream
        std::ifstream is(path, std::ios::in | std::ios::binary);
        if (!is.good()) {
            delete db;
            return nullptr;
        }
        jellyfish::file_header header(is);
        is.close();

        // Set global k-mer length
        db->kmer_length = header.key_len() / 2;
        db->canonical   = header.canonical();
        jellyfish::mer_dna::k(db->kmer_length);

        // Memory-map the file
        db->mapped.reset(new jellyfish::mapped_file(path));

        // Create the raw array view over the mapped data
        char* data_start = db->mapped->base() + header.offset();
        size_t data_bytes = db->mapped->length() - header.offset();

        db->hash.reset(new array_raw_type(
            data_start,
            data_bytes,
            header.size(),
            header.key_len(),
            header.val_len(),
            header.max_reprobe(),
            header.matrix()
        ));

        return db;
    } catch (...) {
        return nullptr;
    }
}

void jf_close(JfDatabase* db) {
    delete db;
}

uint64_t jf_query(const JfDatabase* db, const char* kmer) {
    if (!db || !db->hash || !kmer) return 0;

    try {
        jellyfish::mer_dna::k(db->kmer_length);
        jellyfish::mer_dna mer(kmer);

        if (db->canonical) {
            mer.canonicalize();
        }

        uint64_t val = 0;
        if (db->hash->get_val_for_key(mer, &val)) {
            return val;
        }
        return 0;
    } catch (...) {
        return 0;
    }
}

uint32_t jf_kmer_length(const JfDatabase* db) {
    if (!db) return 0;
    return db->kmer_length;
}

bool jf_canonical(const char* kmer, uint32_t k, char* out) {
    if (!kmer || !out) return false;

    try {
        jellyfish::mer_dna::k(k);
        jellyfish::mer_dna mer(kmer);
        mer.canonicalize();

        std::string s = mer.to_str();
        if (s.size() != k) return false;
        std::memcpy(out, s.c_str(), k);
        out[k] = '\0';
        return true;
    } catch (...) {
        return false;
    }
}

} // extern "C"
