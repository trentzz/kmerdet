#ifndef JF_WRAPPER_H
#define JF_WRAPPER_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/// Opaque handle to a jellyfish database
typedef struct JfDatabase JfDatabase;

/// Open a jellyfish database file. Returns NULL on failure.
JfDatabase* jf_open(const char* path);

/// Close a jellyfish database and free resources.
void jf_close(JfDatabase* db);

/// Query the count of a k-mer string. Returns 0 if not found.
uint64_t jf_query(const JfDatabase* db, const char* kmer);

/// Get the k-mer length from the database.
uint32_t jf_kmer_length(const JfDatabase* db);

/// Canonicalize a k-mer string (min of forward and reverse complement).
/// Writes result to `out` buffer which must be at least `k+1` bytes.
/// Returns true on success.
bool jf_canonical(const char* kmer, uint32_t k, char* out);

#ifdef __cplusplus
}
#endif

#endif // JF_WRAPPER_H
