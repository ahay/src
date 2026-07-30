/* Pipe-context helpers for sfmdioread / sfmdiowrite.  See mdio_pipe.hh. */

#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "absl/strings/cord.h"
#include "mdio_pipe.hh"
#include "tensorstore/array.h"
#include "tensorstore/kvstore/kvstore.h"
#include "tensorstore/kvstore/operations.h"
#include "tensorstore/tensorstore.h"

/* ------------------------------------------------------------------ */
/* Tiny self-contained SHA-256 (public-domain style compact impl).    */
/* No OpenSSL / no tensorstore hash dependency.                       */
/* ------------------------------------------------------------------ */
namespace {

typedef unsigned char  mdio_u8;
typedef unsigned int   mdio_u32;
typedef unsigned long long mdio_u64;

struct Sha256Ctx {
    mdio_u32 state[8];
    mdio_u64 bitlen;
    mdio_u8  data[64];
    size_t   datalen;
};

static mdio_u32 rotr(mdio_u32 x, mdio_u32 n) { return (x >> n) | (x << (32 - n)); }

static void sha256_transform(Sha256Ctx* ctx, const mdio_u8 data[64])
{
    static const mdio_u32 k[64] = {
        0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
        0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
        0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
        0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
        0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
        0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
        0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
        0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
    };
    mdio_u32 m[64], a, b, c, d, e, f, g, h, t1, t2;
    for (int i = 0, j = 0; i < 16; i++, j += 4)
        m[i] = ((mdio_u32)data[j] << 24) | ((mdio_u32)data[j+1] << 16) |
               ((mdio_u32)data[j+2] << 8) | ((mdio_u32)data[j+3]);
    for (int i = 16; i < 64; i++) {
        mdio_u32 s0 = rotr(m[i-15], 7) ^ rotr(m[i-15], 18) ^ (m[i-15] >> 3);
        mdio_u32 s1 = rotr(m[i-2], 17) ^ rotr(m[i-2], 19) ^ (m[i-2] >> 10);
        m[i] = m[i-16] + s0 + m[i-7] + s1;
    }
    a = ctx->state[0]; b = ctx->state[1]; c = ctx->state[2]; d = ctx->state[3];
    e = ctx->state[4]; f = ctx->state[5]; g = ctx->state[6]; h = ctx->state[7];
    for (int i = 0; i < 64; i++) {
        mdio_u32 S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        mdio_u32 ch = (e & f) ^ ((~e) & g);
        t1 = h + S1 + ch + k[i] + m[i];
        mdio_u32 S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        mdio_u32 maj = (a & b) ^ (a & c) ^ (b & c);
        t2 = S0 + maj;
        h = g; g = f; f = e; e = d + t1;
        d = c; c = b; b = a; a = t1 + t2;
    }
    ctx->state[0] += a; ctx->state[1] += b; ctx->state[2] += c; ctx->state[3] += d;
    ctx->state[4] += e; ctx->state[5] += f; ctx->state[6] += g; ctx->state[7] += h;
}

static void sha256_init(Sha256Ctx* ctx)
{
    ctx->datalen = 0;
    ctx->bitlen  = 0;
    ctx->state[0] = 0x6a09e667; ctx->state[1] = 0xbb67ae85;
    ctx->state[2] = 0x3c6ef372; ctx->state[3] = 0xa54ff53a;
    ctx->state[4] = 0x510e527f; ctx->state[5] = 0x9b05688c;
    ctx->state[6] = 0x1f83d9ab; ctx->state[7] = 0x5be0cd19;
}

static void sha256_update(Sha256Ctx* ctx, const mdio_u8* data, size_t len)
{
    for (size_t i = 0; i < len; i++) {
        ctx->data[ctx->datalen++] = data[i];
        if (ctx->datalen == 64) {
            sha256_transform(ctx, ctx->data);
            ctx->bitlen += 512;
            ctx->datalen = 0;
        }
    }
}

static void sha256_final(Sha256Ctx* ctx, mdio_u8 hash[32])
{
    size_t i = ctx->datalen;
    if (ctx->datalen < 56) {
        ctx->data[i++] = 0x80;
        while (i < 56) ctx->data[i++] = 0x00;
    } else {
        ctx->data[i++] = 0x80;
        while (i < 64) ctx->data[i++] = 0x00;
        sha256_transform(ctx, ctx->data);
        memset(ctx->data, 0, 56);
    }
    ctx->bitlen += (mdio_u64) ctx->datalen * 8;
    ctx->data[63] = (mdio_u8)(ctx->bitlen);
    ctx->data[62] = (mdio_u8)(ctx->bitlen >> 8);
    ctx->data[61] = (mdio_u8)(ctx->bitlen >> 16);
    ctx->data[60] = (mdio_u8)(ctx->bitlen >> 24);
    ctx->data[59] = (mdio_u8)(ctx->bitlen >> 32);
    ctx->data[58] = (mdio_u8)(ctx->bitlen >> 40);
    ctx->data[57] = (mdio_u8)(ctx->bitlen >> 48);
    ctx->data[56] = (mdio_u8)(ctx->bitlen >> 56);
    sha256_transform(ctx, ctx->data);
    for (i = 0; i < 4; i++) {
        hash[i]      = (mdio_u8)((ctx->state[0] >> (24 - i * 8)) & 0xff);
        hash[i + 4]  = (mdio_u8)((ctx->state[1] >> (24 - i * 8)) & 0xff);
        hash[i + 8]  = (mdio_u8)((ctx->state[2] >> (24 - i * 8)) & 0xff);
        hash[i + 12] = (mdio_u8)((ctx->state[3] >> (24 - i * 8)) & 0xff);
        hash[i + 16] = (mdio_u8)((ctx->state[4] >> (24 - i * 8)) & 0xff);
        hash[i + 20] = (mdio_u8)((ctx->state[5] >> (24 - i * 8)) & 0xff);
        hash[i + 24] = (mdio_u8)((ctx->state[6] >> (24 - i * 8)) & 0xff);
        hash[i + 28] = (mdio_u8)((ctx->state[7] >> (24 - i * 8)) & 0xff);
    }
}

static std::string sha256_hex(const std::string& msg)
{
    Sha256Ctx ctx;
    mdio_u8 hash[32];
    sha256_init(&ctx);
    sha256_update(&ctx, (const mdio_u8*) msg.data(), msg.size());
    sha256_final(&ctx, hash);
    char hex[65];
    for (int i = 0; i < 32; i++)
        snprintf(hex + i * 2, 3, "%02x", hash[i]);
    return std::string(hex);
}

/* Dataset-level name from mdio-cpp metadata (handles nested / flat layouts). */
static std::string dataset_name(mdio::Dataset& ds)
{
    const nlohmann::json& meta = ds.getMetadata();
    if (meta.contains("name") && meta["name"].is_string())
        return meta["name"].get<std::string>();
    if (meta.contains("metadata") && meta["metadata"].is_object() &&
        meta["metadata"].contains("name") && meta["metadata"]["name"].is_string())
        return meta["metadata"]["name"].get<std::string>();
    return "";
}

/* Canonical descriptor — must match between reader and writer. */
static std::string canonical_descriptor(mdio::Dataset& ds,
                                        const std::string& datavar,
                                        const std::string& dtype)
{
    std::ostringstream os;
    os.setf(std::ios::fixed);
    os.precision(17);

    os << dataset_name(ds) << '\n';
    os << datavar << '\n';
    os << dtype << '\n';

    std::vector<MdioAxis> axes = mdio_axes(ds, datavar);
    for (size_t i = 0; i < axes.size(); i++) {
        const MdioAxis& a = axes[i];
        os << a.label << '\t' << a.size << '\t'
           << a.o << '\t' << a.d << '\t' << a.unit << '\n';
    }
    return os.str();
}

} /* anonymous namespace */

/* ------------------------------------------------------------------ */
/* Public API                                                         */
/* ------------------------------------------------------------------ */

std::string mdio_pipe_fingerprint(mdio::Dataset& ds,
                                  const std::string& datavar,
                                  const std::string& dtype)
{
    return std::string("sha256:") +
           sha256_hex(canonical_descriptor(ds, datavar, dtype));
}

void mdio_pipe_stamp(sf_file out,
                     const std::string& source,
                     const std::string& fingerprint,
                     const std::string& datavar,
                     const std::string& dtype,
                     const std::string& contract,
                     const std::string& version,
                     const std::string& selection)
{
    sf_putstring(out, MDIO_PIPE_KEY_SOURCE,      source.c_str());
    sf_putstring(out, MDIO_PIPE_KEY_FINGERPRINT, fingerprint.c_str());
    sf_putstring(out, MDIO_PIPE_KEY_DATA_VAR,    datavar.c_str());
    sf_putstring(out, MDIO_PIPE_KEY_DTYPE,       dtype.c_str());
    sf_putstring(out, MDIO_PIPE_KEY_CONTRACT,    contract.c_str());
    sf_putstring(out, MDIO_PIPE_KEY_VERSION,     version.c_str());
    if (!selection.empty())
        sf_putstring(out, MDIO_PIPE_KEY_SELECTION, selection.c_str());
}

void mdio_pipe_read_context(sf_file in, MdioPipeContext& ctx)
{
    ctx.source.clear();
    ctx.fingerprint.clear();
    ctx.data_variable.clear();
    ctx.dtype.clear();
    ctx.contract.clear();
    ctx.version.clear();
    ctx.selection.clear();
    ctx.present = false;

    char* s;
    if ((s = sf_histstring(in, MDIO_PIPE_KEY_SOURCE))) {
        ctx.source = s; free(s); ctx.present = true;
    }
    if ((s = sf_histstring(in, MDIO_PIPE_KEY_FINGERPRINT))) {
        ctx.fingerprint = s; free(s);
    }
    if ((s = sf_histstring(in, MDIO_PIPE_KEY_DATA_VAR))) {
        ctx.data_variable = s; free(s);
    }
    if ((s = sf_histstring(in, MDIO_PIPE_KEY_DTYPE))) {
        ctx.dtype = s; free(s);
    }
    if ((s = sf_histstring(in, MDIO_PIPE_KEY_CONTRACT))) {
        ctx.contract = s; free(s);
    }
    if ((s = sf_histstring(in, MDIO_PIPE_KEY_VERSION))) {
        ctx.version = s; free(s);
    }
    if ((s = sf_histstring(in, MDIO_PIPE_KEY_SELECTION))) {
        ctx.selection = s; free(s);
    }
}

bool mdio_pipe_context_from_parent(mdio::Dataset& ds,
                                   const std::string& parent_uri,
                                   const std::string& datavar,
                                   MdioPipeContext& ctx)
{
    std::string dtype = mdio_variable_dtype(ds, datavar);
    if (dtype != MDIO_PIPE_DTYPE_FLOAT32) return false;
    ctx.source        = parent_uri;
    ctx.data_variable = datavar;
    ctx.dtype         = dtype;
    ctx.contract      = MDIO_PIPE_CONTRACT_SAME_GEOM;
    ctx.version       = MDIO_PIPE_CONTEXT_VERSION;
    ctx.fingerprint   = mdio_pipe_fingerprint(ds, datavar, ctx.dtype);
    ctx.present       = true;
    return true;
}

/* ------------------------------------------------------------------ */
/* Phase 2/3 finalizer helpers                                        */
/* ------------------------------------------------------------------ */

namespace fs = std::filesystem;

bool mdio_pipe_is_remote_uri(const std::string& uri)
{
    return uri.rfind("s3://", 0) == 0 ||
           uri.rfind("gs://", 0) == 0 ||
           uri.rfind("gcs://", 0) == 0 ||
           uri.rfind("http://", 0) == 0 ||
           uri.rfind("https://", 0) == 0;
}

/* Open the root kvstore for an MDIO URI (local file, s3://, gs://). */
static bool open_root_kvstore(const std::string& uri,
                              tensorstore::KvStore& out,
                              std::string& err)
{
    auto f = mdio::internal::dataset_kvs_store(uri);
    if (!f.status().ok()) {
        err = "cannot open kvstore for \"" + uri + "\": " +
              f.status().ToString();
        return false;
    }
    out = f.value();
    return true;
}

/* Delete every key under a dataset URI (local or remote). */
static bool kvstore_delete_all(const std::string& uri, std::string& err)
{
    tensorstore::KvStore kvs;
    if (!open_root_kvstore(uri, kvs, err)) return false;
    auto del = tensorstore::kvstore::DeleteRange(kvs, {}).result();
    if (!del.ok()) {
        err = "kvstore delete failed for \"" + uri +
              "\": " + del.status().ToString();
        return false;
    }
    return true;
}

static bool is_zarr_var_metadata_name(const std::string& name)
{
    return name == "zarr.json" || name == ".zarray" || name == ".zattrs" ||
           name == ".zgroup" || name == ".zmetadata";
}

/* True when kvstore key is an amplitude (skip_var) chunk payload rather than
   the variable's metadata.  Keys look like "amplitude/zarr.json" or
   "amplitude/c/0/1/..." (leading slash optional). */
static bool kvstore_key_is_skipped_chunk(const std::string& key,
                                         const std::string& skip_var)
{
    if (skip_var.empty()) return false;
    std::string k = key;
    if (!k.empty() && k[0] == '/') k.erase(0, 1);
    const std::string prefix = skip_var + "/";
    if (k.compare(0, prefix.size(), prefix) != 0) return false;
    const std::string rest = k.substr(prefix.size());
    const size_t slash = rest.find('/');
    if (slash == std::string::npos)
        return !is_zarr_var_metadata_name(rest);
    /* Nested under skip_var/ → chunk payload (e.g. c/...). */
    return true;
}

/* Byte-copy keys from src → dst, optionally skipping skip_var chunk payloads. */
static bool kvstore_copy_all(const tensorstore::KvStore& src,
                             const tensorstore::KvStore& dst,
                             std::string& err,
                             const std::string& skip_var = std::string())
{
    auto list = tensorstore::kvstore::ListFuture(src).result();
    if (!list.ok()) {
        err = "kvstore list failed: " + list.status().ToString();
        return false;
    }
    for (const auto& entry : *list) {
        if (kvstore_key_is_skipped_chunk(entry.key, skip_var)) continue;
        auto rd = tensorstore::kvstore::Read(src, entry.key).result();
        if (!rd.ok()) {
            err = "kvstore read \"" + entry.key +
                  "\" failed: " + rd.status().ToString();
            return false;
        }
        if (!rd->has_value()) continue;
        auto wr = tensorstore::kvstore::Write(dst, entry.key, rd->value)
                      .result();
        if (!wr.ok()) {
            err = "kvstore write \"" + entry.key +
                  "\" failed: " + wr.status().ToString();
            return false;
        }
    }
    return true;
}

/* Recursive local copy; when at_root and a directory named skip_var is
   found, copy only its metadata files (not chunk payloads under c/). */
static bool copy_tree_skip_var(const fs::path& src, const fs::path& dst,
                               const std::string& skip_var, bool at_root,
                               std::string& err)
{
    std::error_code ec;
    fs::create_directories(dst, ec);
    if (ec) {
        err = "mkdir \"" + dst.string() + "\": " + ec.message();
        return false;
    }
    for (auto& ent : fs::directory_iterator(src, ec)) {
        if (ec) {
            err = "list \"" + src.string() + "\": " + ec.message();
            return false;
        }
        const std::string name = ent.path().filename().string();
        const fs::path dest = dst / name;
        if (ent.is_directory(ec)) {
            if (at_root && !skip_var.empty() && name == skip_var) {
                fs::create_directories(dest, ec);
                if (ec) {
                    err = "mkdir \"" + dest.string() + "\": " + ec.message();
                    return false;
                }
                std::error_code ec2;
                for (auto& f : fs::directory_iterator(ent.path(), ec2)) {
                    if (ec2 || !f.is_regular_file(ec2)) continue;
                    const std::string fn = f.path().filename().string();
                    if (!is_zarr_var_metadata_name(fn)) continue;
                    fs::copy_file(f.path(), dest / fn,
                                  fs::copy_options::overwrite_existing, ec2);
                    if (ec2) {
                        err = "copy \"" + f.path().string() +
                              "\": " + ec2.message();
                        return false;
                    }
                }
            } else if (!copy_tree_skip_var(ent.path(), dest, skip_var,
                                           /*at_root=*/false, err)) {
                return false;
            }
        } else if (ent.is_regular_file(ec)) {
            fs::copy_file(ent.path(), dest,
                          fs::copy_options::overwrite_existing, ec);
            if (ec) {
                err = "copy \"" + ent.path().string() + "\": " + ec.message();
                return false;
            }
        }
    }
    return true;
}

void mdio_pipe_remove_path(const std::string& path)
{
    if (mdio_pipe_is_remote_uri(path)) {
        std::string err;
        kvstore_delete_all(path, err); /* best-effort */
        return;
    }
    std::error_code ec;
    fs::remove_all(path, ec);
}

bool mdio_clone_store_skip_var(const std::string& parent_uri,
                               const std::string& child_path,
                               const std::string& skip_var,
                               std::string& err)
{
    /* Fast path: local → local selective filesystem copy. */
    if (!mdio_pipe_is_remote_uri(parent_uri) &&
        !mdio_pipe_is_remote_uri(child_path)) {
        std::error_code ec;
        fs::path src(parent_uri);
        fs::path dst(child_path);
        if (!fs::exists(src, ec) || !fs::is_directory(src, ec)) {
            err = "parent MDIO path does not exist or is not a directory: " +
                  parent_uri;
            return false;
        }

        mdio_pipe_remove_path(child_path);
        fs::create_directories(dst.parent_path(), ec);
        if (skip_var.empty()) {
            fs::copy(src, dst,
                     fs::copy_options::recursive |
                         fs::copy_options::overwrite_existing,
                     ec);
            if (ec) {
                err = "clone failed: " + ec.message();
                mdio_pipe_remove_path(child_path);
                return false;
            }
            return true;
        }
        if (!copy_tree_skip_var(src, dst, skip_var, /*at_root=*/true, err)) {
            mdio_pipe_remove_path(child_path);
            return false;
        }
        return true;
    }

    /* Cloud / mixed: kvstore List + Read + Write (byte-identical clone,
       optionally skipping skip_var chunk payloads). */
    mdio_pipe_remove_path(child_path);

    if (!mdio_pipe_is_remote_uri(child_path)) {
        std::error_code ec;
        fs::create_directories(fs::path(child_path).parent_path(), ec);
    }

    tensorstore::KvStore src, dst;
    if (!open_root_kvstore(parent_uri, src, err)) return false;
    if (!open_root_kvstore(child_path, dst, err)) return false;
    if (!kvstore_copy_all(src, dst, err, skip_var)) {
        mdio_pipe_remove_path(child_path);
        return false;
    }
    return true;
}

bool mdio_clone_store(const std::string& parent_uri,
                      const std::string& child_path,
                      std::string& err)
{
    return mdio_clone_store_skip_var(parent_uri, child_path,
                                     /*skip_var=*/std::string(), err);
}

bool mdio_pipe_publish(const std::string& tmp_path,
                       const std::string& final_path,
                       std::string& err)
{
    /* Local rename (atomic on same filesystem). */
    if (!mdio_pipe_is_remote_uri(tmp_path) &&
        !mdio_pipe_is_remote_uri(final_path)) {
        std::error_code ec;
        mdio_pipe_remove_path(final_path);
        fs::create_directories(fs::path(final_path).parent_path(), ec);
        fs::rename(tmp_path, final_path, ec);
        if (ec) {
            /* Cross-device rename can fail; fall back to copy + remove. */
            ec.clear();
            fs::copy(tmp_path, final_path,
                     fs::copy_options::recursive |
                         fs::copy_options::overwrite_existing,
                     ec);
            if (ec) {
                err = "publish rename/copy failed: " + ec.message();
                return false;
            }
            mdio_pipe_remove_path(tmp_path);
        }
        return true;
    }

    /* Cloud publish: copy tmp → final, then delete tmp.
       S3/GCS have no atomic directory rename — there is a short window where
       both prefixes may exist. */
    mdio_pipe_remove_path(final_path);
    tensorstore::KvStore src, dst;
    if (!open_root_kvstore(tmp_path, src, err)) return false;
    if (!open_root_kvstore(final_path, dst, err)) return false;
    if (!kvstore_copy_all(src, dst, err)) return false;
    mdio_pipe_remove_path(tmp_path);
    return true;
}

float mdio_pipe_dead_fill(void)
{
    return std::numeric_limits<float>::quiet_NaN();
}

bool mdio_pipe_read_trace_mask(mdio::Dataset& ds,
                               std::vector<char>& live,
                               const std::string& mask_name)
{
    if (!ds.variables.contains_key(mask_name)) {
        /* No mask → treat every cell as live.  Caller sizes live. */
        for (size_t i = 0; i < live.size(); i++) live[i] = 1;
        return true;
    }

    std::vector<double> vals;
    if (!mdio_read_field(ds, mask_name, vals)) return false;

    live.resize(vals.size());
    for (size_t i = 0; i < vals.size(); i++)
        live[i] = (vals[i] != 0.0) ? 1 : 0;
    return true;
}

void mdio_pipe_apply_dead_fill(float* buf, long first_trace, long ntr, long ns,
                               const std::vector<char>& live)
{
    if (!buf || ntr <= 0 || ns <= 0) return;
    const float fill = mdio_pipe_dead_fill();
    const long nmask = (long) live.size();
    for (long t = 0; t < ntr; t++) {
        const long gt = first_trace + t;
        if (gt >= nmask || live[(size_t) gt]) continue;
        float* tr = buf + t * ns;
        for (long i = 0; i < ns; i++) tr[i] = fill;
    }
}

void MdioStatsAcc::update(const float* buf, long n)
{
    if (!buf || n <= 0) return;
    for (long i = 0; i < n; i++) {
        float x = buf[i];
        if (!std::isfinite(x)) continue;
        double xd = (double) x;
        if (std::fabs(xd) <= MDIO_PIPE_STATS_ZERO_ATOL) continue;
        if (!have) { vmin = vmax = xd; have = true; }
        else {
            if (xd < vmin) vmin = xd;
            if (xd > vmax) vmax = xd;
        }
        sum += xd;
        sum_sq += xd * xd;
        count++;
    }
}

std::string MdioStatsAcc::finalize_json() const
{
    nlohmann::json j;
    if (count == 0) {
        j["count"] = 0;
        j["min"] = 0.0;
        j["max"] = 0.0;
        j["sum"] = 0.0;
        j["sumSquares"] = 0.0;
    } else {
        j["count"] = count;
        j["min"] = vmin;
        j["max"] = vmax;
        j["sum"] = sum;
        j["sumSquares"] = sum_sq;
    }
    j["histogram"] = nlohmann::json::object({
        {"binCenters", nlohmann::json::array()},
        {"counts",     nlohmann::json::array()}
    });
    return j.dump();
}

std::string mdio_pipe_statsv1_json(const float* buf, long total)
{
    MdioStatsAcc acc;
    acc.update(buf, total);
    return acc.finalize_json();
}

bool mdio_pipe_restamp_stats(mdio::Dataset& ds,
                             const std::string& datavar,
                             const std::string& stats_json,
                             std::string& err)
{
    /* Use at() (not get<>) so UpdateAttributes mutates the shared
       UserAttributes pointer held by the collection.  Do NOT variables.add
       afterwards — a re-insert resets attributesAddress and CommitMetadata
       then sees "No variables were modified." */
    auto vr = ds.variables.at(datavar);
    if (!vr.status().ok()) {
        err = "cannot open data variable \"" + datavar + "\" for stats restamp";
        return false;
    }
    auto var = vr.value();

    nlohmann::json stats;
    try {
        stats = nlohmann::json::parse(stats_json);
    } catch (const std::exception& e) {
        err = std::string("statsV1 JSON parse failed: ") + e.what();
        return false;
    }

    nlohmann::json attrs = var.GetAttributes();
    /* GetAttributes shape (USER_GUIDE): top-level statsV1 / unitsV1 /
       attributes.  Always set top-level statsV1. */
    attrs["statsV1"] = stats;

    auto upd = var.UpdateAttributes(attrs);
    if (!upd.status().ok()) {
        err = "UpdateAttributes(statsV1) failed: " + upd.status().ToString();
        return false;
    }

    auto commit = ds.CommitMetadata();
    if (!commit.status().ok()) {
        err = "CommitMetadata failed: " + commit.status().ToString();
        return false;
    }
    return true;
}

/* Locate and parse the root metadata document of a dataset URI via kvstore
   (works for local file, s3://, gs://).  Prefers zarr.json (v3) then .zattrs
   (v2).  Sets meta_key to the key that was read. */
static bool load_root_meta_kv(const std::string& child_path,
                              tensorstore::KvStore& kvs,
                              std::string& meta_key,
                              nlohmann::json& root,
                              std::string& err)
{
    if (!open_root_kvstore(child_path, kvs, err)) return false;

    const char* candidates[] = {"zarr.json", ".zattrs", nullptr};
    for (int i = 0; candidates[i]; i++) {
        auto rd = tensorstore::kvstore::Read(kvs, candidates[i]).result();
        if (!rd.ok() || !rd->has_value()) continue;
        try {
            std::string text = std::string(rd->value);
            root = nlohmann::json::parse(text);
            meta_key = candidates[i];
            return true;
        } catch (const std::exception& e) {
            err = std::string("failed to parse ") + candidates[i] + ": " +
                  e.what();
            return false;
        }
    }
    err = "no root zarr.json / .zattrs under " + child_path;
    return false;
}

static std::string utc_now_created_on(void)
{
    /* Match Python writer: "%Y-%m-%d %H:%M:%S.%f+00:00" */
    using namespace std::chrono;
    auto now = system_clock::now();
    std::time_t t = system_clock::to_time_t(now);
    auto us = duration_cast<microseconds>(now.time_since_epoch()) % 1000000;
    struct tm g;
    gmtime_r(&t, &g);
    char buf[64];
    std::snprintf(buf, sizeof(buf),
                  "%04d-%02d-%02d %02d:%02d:%02d.%06lld+00:00",
                  g.tm_year + 1900, g.tm_mon + 1, g.tm_mday,
                  g.tm_hour, g.tm_min, g.tm_sec,
                  (long long) us.count());
    return std::string(buf);
}

bool mdio_pipe_patch_provenance(const std::string& child_path,
                                const nlohmann::json* processing,
                                const char* created_on,
                                std::string& err)
{
    tensorstore::KvStore kvs;
    std::string meta_key;
    nlohmann::json root;
    if (!load_root_meta_kv(child_path, kvs, meta_key, root, err)) return false;

    const std::string created =
        (created_on && *created_on) ? std::string(created_on)
                                   : utc_now_created_on();

    /* Zarr v3 nests dataset attrs under root["attributes"]; v2 .zattrs is
       already the attr document. */
    nlohmann::json* attrs = &root;
    if (meta_key == "zarr.json") {
        if (!root.contains("attributes") || !root["attributes"].is_object())
            root["attributes"] = nlohmann::json::object();
        attrs = &root["attributes"];
    }

    (*attrs)["createdOn"] = created;

    if (processing != NULL) {
        nlohmann::json nested;
        if (attrs->contains("attributes") && (*attrs)["attributes"].is_object())
            nested = (*attrs)["attributes"];
        else
            nested = nlohmann::json::object();

        nlohmann::json proc_list = nlohmann::json::array();
        if (nested.contains("processing")) {
            if (nested["processing"].is_array())
                proc_list = nested["processing"];
            else
                proc_list.push_back(nested["processing"]);
        }
        proc_list.push_back(*processing);
        nested["processing"] = proc_list;
        (*attrs)["attributes"] = nested;
    }

    try {
        std::string text = root.dump(2) + "\n";
        auto wr = tensorstore::kvstore::Write(kvs, meta_key, absl::Cord(text))
                      .result();
        if (!wr.ok()) {
            err = "provenance write failed: " + wr.status().ToString();
            return false;
        }
    } catch (const std::exception& e) {
        err = std::string("provenance write failed: ") + e.what();
        return false;
    }
    return true;
}

bool mdio_pipe_assert_same_geometry(sf_file in,
                                    const std::vector<MdioAxis>& parent_axes,
                                    std::string& err)
{
    off_t n[SF_MAX_DIM];
    /* sf_largefiledims returns the highest axis with n#>1, so a trailing
       RSF unit axis (from a size-1 MDIO slow axis like sail_line=1) does
       not inflate the reported rank.  Compare against the parent's
       effective rank the same way. */
    int rsf_dim = sf_largefiledims(in, n);
    const int prank = (int) parent_axes.size();
    if (prank < 1) {
        err = "parent has no axes";
        return false;
    }

    int expect_dim = 1;
    for (int ri = 1; ri <= prank && ri <= SF_MAX_DIM; ri++) {
        int a = prank - ri; /* MDIO axis for RSF axis ri */
        if (parent_axes[(size_t) a].size > 1) expect_dim = ri;
    }
    if (rsf_dim != expect_dim) {
        err = "RSF effective rank " + std::to_string(rsf_dim) +
              " != parent effective rank " + std::to_string(expect_dim) +
              " (parent has " + std::to_string(prank) + " axes; "
              "geometry-changing ops need Phase 4 adapters)";
        return false;
    }

    const double atol = 1e-5;
    for (int a = 0; a < prank; a++) {
        int ri = prank - a; /* RSF axis number (1-based) */
        const MdioAxis& pax = parent_axes[(size_t) a];
        long rsf_n = (ri >= 1 && ri <= SF_MAX_DIM) ? (long) n[ri - 1] : 1;
        if (rsf_n != pax.size) {
            err = "axis " + pax.label + " size mismatch: RSF n" +
                  std::to_string(ri) + "=" + std::to_string(rsf_n) +
                  " vs parent " + std::to_string(pax.size) +
                  " (geometry-changing ops need Phase 4 adapters)";
            return false;
        }

        /* Size-1 axes often lack o#/d# in the RSF header after rank trim;
           only enforce origin/sampling when the axis contributes to rank. */
        if (pax.size <= 1) continue;

        char key[8];
        float fo = 0.f, fd = 1.f;
        snprintf(key, sizeof(key), "o%d", ri);
        double rsf_o = sf_histfloat(in, key, &fo) ? (double) fo : 0.0;
        snprintf(key, sizeof(key), "d%d", ri);
        double rsf_d = sf_histfloat(in, key, &fd) ? (double) fd : 1.0;

        if (std::fabs(rsf_o - pax.o) > atol * (1.0 + std::fabs(pax.o)) &&
            std::fabs(rsf_o - pax.o) > 1e-3) {
            err = "axis " + pax.label + " origin mismatch: RSF o" +
                  std::to_string(ri) + "=" + std::to_string(rsf_o) +
                  " vs parent " + std::to_string(pax.o) +
                  " (geometry-changing ops need Phase 4 adapters)";
            return false;
        }
        if (std::fabs(rsf_d - pax.d) > atol * (1.0 + std::fabs(pax.d)) &&
            std::fabs(rsf_d - pax.d) > 1e-6) {
            err = "axis " + pax.label + " sampling mismatch: RSF d" +
                  std::to_string(ri) + "=" + std::to_string(rsf_d) +
                  " vs parent " + std::to_string(pax.d) +
                  " (geometry-changing ops need Phase 4 adapters)";
            return false;
        }
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* Phase 4 geometry adapters                                          */
/* ------------------------------------------------------------------ */

void mdio_pipe_selection_finalize(MdioSelection& sel)
{
    sel.same_geometry = true;
    sel.sample_changed = false;
    sel.spatial_changed = false;
    sel.exact = true;
    const int prank = (int) sel.axes.size();
    for (int a = 0; a < prank; a++) {
        const MdioAxisSel& ax = sel.axes[(size_t) a];
        if (ax.kind == MDIO_SEL_RESAMPLED) sel.exact = false;
        if (ax.kind != MDIO_SEL_PRESERVED) {
            sel.same_geometry = false;
            if (a == prank - 1) sel.sample_changed = true;
            else                sel.spatial_changed = true;
        }
    }

    if (sel.same_geometry) {
        sel.contract_name = MDIO_PIPE_CONTRACT_SAME_GEOM;
    } else if (sel.sample_changed && !sel.spatial_changed) {
        /* Sample-axis only: truncate (step==1), integer decimate, or
           interpolating resample.  Both decimate and RESAMPLED map to the
           "resample" fidelity contract name. */
        bool resample = false;
        if (prank >= 1) {
            MdioAxisSelKind k = sel.axes[(size_t) (prank - 1)].kind;
            resample = (k == MDIO_SEL_DECIMATED || k == MDIO_SEL_RESAMPLED);
        }
        sel.contract_name = resample ? MDIO_PIPE_CONTRACT_RESAMPLE
                                     : MDIO_PIPE_CONTRACT_TRUNCATE;
    } else if (!sel.sample_changed && sel.spatial_changed) {
        bool decim = false;
        for (int a = 0; a < prank - 1; a++)
            if (sel.axes[(size_t) a].kind == MDIO_SEL_DECIMATED) {
                decim = true; break;
            }
        sel.contract_name = decim ? MDIO_PIPE_CONTRACT_DECIMATE
                                  : MDIO_PIPE_CONTRACT_SLICE;
    } else {
        /* Combined sample+spatial change: prefer the stricter name for
           provenance; writer still applies the full selection. */
        bool decim = false;
        for (int a = 0; a < prank; a++) {
            MdioAxisSelKind k = sel.axes[(size_t) a].kind;
            if (k == MDIO_SEL_DECIMATED || k == MDIO_SEL_RESAMPLED) {
                decim = true; break;
            }
        }
        sel.contract_name = decim ? MDIO_PIPE_CONTRACT_DECIMATE
                                  : MDIO_PIPE_CONTRACT_SLICE;
    }
}

std::string mdio_pipe_selection_encode(const MdioSelection& sel)
{
    std::ostringstream os;
    for (size_t i = 0; i < sel.axes.size(); i++) {
        if (i) os << ';';
        const MdioAxisSel& a = sel.axes[i];
        os << a.label << ':' << a.start << ':' << a.step << ':' << a.count;
    }
    return os.str();
}

bool mdio_pipe_selection_parse(const std::string& text,
                               const std::vector<MdioAxis>& parent_axes,
                               MdioSelection& sel,
                               std::string& err)
{
    sel = MdioSelection();
    if (text.empty()) {
        err = "empty mdio_sel";
        return false;
    }
    const int prank = (int) parent_axes.size();
    sel.axes.resize((size_t) prank);

    /* Split on ';' */
    std::vector<std::string> parts;
    {
        std::string cur;
        for (size_t i = 0; i <= text.size(); i++) {
            if (i == text.size() || text[i] == ';') {
                if (!cur.empty()) parts.push_back(cur);
                cur.clear();
            } else {
                cur.push_back(text[i]);
            }
        }
    }
    if ((int) parts.size() != prank) {
        err = "mdio_sel has " + std::to_string(parts.size()) +
              " axes, parent has " + std::to_string(prank);
        return false;
    }

    for (int a = 0; a < prank; a++) {
        const std::string& p = parts[(size_t) a];
        /* label:start:step:count */
        size_t c1 = p.find(':');
        size_t c2 = (c1 == std::string::npos) ? std::string::npos
                                              : p.find(':', c1 + 1);
        size_t c3 = (c2 == std::string::npos) ? std::string::npos
                                              : p.find(':', c2 + 1);
        if (c1 == std::string::npos || c2 == std::string::npos ||
            c3 == std::string::npos) {
            err = "bad mdio_sel token \"" + p +
                  "\" (want label:start:step:count)";
            return false;
        }
        MdioAxisSel ax;
        ax.label = p.substr(0, c1);
        ax.start = std::atol(p.c_str() + c1 + 1);
        ax.step  = std::atol(p.c_str() + c2 + 1);
        ax.count = std::atol(p.c_str() + c3 + 1);
        if (ax.label != parent_axes[(size_t) a].label) {
            err = "mdio_sel label \"" + ax.label + "\" != parent \"" +
                  parent_axes[(size_t) a].label + "\"";
            return false;
        }
        const MdioAxis& pax = parent_axes[(size_t) a];
        const long n_p = pax.size;
        if (ax.step < 1 || ax.count < 1 || ax.start < 0 ||
            ax.start + (ax.count - 1) * ax.step > n_p - 1) {
            err = "mdio_sel out of range for axis " + ax.label;
            return false;
        }
        if (ax.start == 0 && ax.step == 1 && ax.count == n_p)
            ax.kind = MDIO_SEL_PRESERVED;
        else if (ax.step == 1)
            ax.kind = MDIO_SEL_SUBSET;
        else
            ax.kind = MDIO_SEL_DECIMATED;
        ax.o = pax.o + (double) ax.start * pax.d;
        ax.d = pax.d * (double) ax.step;
        sel.axes[(size_t) a] = ax;
    }
    mdio_pipe_selection_finalize(sel);
    return true;
}

bool mdio_pipe_selection_from_window(const std::vector<MdioAxis>& parent_axes,
                                     const std::vector<long>& f,
                                     const std::vector<long>& jp,
                                     const std::vector<long>& m,
                                     MdioSelection& sel,
                                     std::string& err)
{
    const int prank = (int) parent_axes.size();
    if ((int) f.size() != prank || (int) jp.size() != prank ||
        (int) m.size() != prank) {
        err = "window arrays rank mismatch";
        return false;
    }
    sel = MdioSelection();
    sel.axes.resize((size_t) prank);
    for (int a = 0; a < prank; a++) {
        MdioAxisSel ax;
        const MdioAxis& pax = parent_axes[(size_t) a];
        ax.label = pax.label;
        ax.start = f[(size_t) a];
        ax.step  = jp[(size_t) a] < 1 ? 1 : jp[(size_t) a];
        ax.count = m[(size_t) a];
        const long n_p = pax.size;
        if (ax.start < 0 || ax.count < 1 ||
            ax.start + (ax.count - 1) * ax.step > n_p - 1) {
            err = "window out of range for axis " + ax.label;
            return false;
        }
        if (ax.start == 0 && ax.step == 1 && ax.count == n_p)
            ax.kind = MDIO_SEL_PRESERVED;
        else if (ax.step == 1)
            ax.kind = MDIO_SEL_SUBSET;
        else
            ax.kind = MDIO_SEL_DECIMATED;
        ax.o = pax.o + (double) ax.start * pax.d;
        ax.d = pax.d * (double) ax.step;
        sel.axes[(size_t) a] = ax;
    }
    mdio_pipe_selection_finalize(sel);
    return true;
}

bool mdio_pipe_detect_selection(sf_file in,
                                const std::vector<MdioAxis>& parent_axes,
                                MdioSelection& sel,
                                std::string& err)
{
    off_t n[SF_MAX_DIM];
    int rsf_dim = sf_largefiledims(in, n);
    const int prank = (int) parent_axes.size();
    if (prank < 1) {
        err = "parent has no axes";
        return false;
    }

    /* Effective parent rank (trailing size-1 slow axes may be dropped from
       RSF headers).  Size-1 parent axes that RSF omitted are treated as
       preserved. */
    int expect_dim = 1;
    for (int ri = 1; ri <= prank && ri <= SF_MAX_DIM; ri++) {
        int a = prank - ri;
        if (parent_axes[(size_t) a].size > 1) expect_dim = ri;
    }
    if (rsf_dim != expect_dim) {
        err = "RSF effective rank " + std::to_string(rsf_dim) +
              " != parent effective rank " + std::to_string(expect_dim) +
              " (unsupported geometry change)";
        return false;
    }

    const double rtol = 1e-4;
    sel = MdioSelection();
    sel.axes.resize((size_t) prank);
    sel.exact = true;

    for (int a = 0; a < prank; a++) {
        int ri = prank - a;
        const MdioAxis& pax = parent_axes[(size_t) a];
        MdioAxisSel ax;
        ax.label = pax.label;

        long rsf_n = (ri >= 1 && ri <= SF_MAX_DIM) ? (long) n[ri - 1] : 1;
        /* Size-1 / RSF-trimmed axes: preserved. */
        if (pax.size <= 1) {
            ax.start = 0;
            ax.step  = 1;
            ax.count = pax.size < 1 ? 1 : pax.size;
            ax.kind  = MDIO_SEL_PRESERVED;
            ax.o = pax.o;
            ax.d = pax.d;
            sel.axes[(size_t) a] = ax;
            continue;
        }

        char key[8];
        float fo = 0.f, fd = 1.f;
        snprintf(key, sizeof(key), "o%d", ri);
        double rsf_o = sf_histfloat(in, key, &fo) ? (double) fo : pax.o;
        snprintf(key, sizeof(key), "d%d", ri);
        double rsf_d = sf_histfloat(in, key, &fd) ? (double) fd : pax.d;

        if (pax.d == 0.0) {
            err = "parent axis " + pax.label + " has d=0";
            return false;
        }

        double start_f = (rsf_o - pax.o) / pax.d;
        double step_f  = rsf_d / pax.d;
        long start = (long) std::llround(start_f);
        long step  = (long) std::llround(step_f);
        if (step < 1) step = 1;

        auto close_enough = [rtol](double a, double b) {
            return std::fabs(a - b) <= rtol * (1.0 + std::fabs(b)) ||
                   std::fabs(a - b) <= 1e-6;
        };

        if (!close_enough(start_f, (double) start) ||
            !close_enough(step_f, (double) step) ||
            !close_enough(pax.o + start * pax.d, rsf_o) ||
            !close_enough(pax.d * step, rsf_d)) {
            /* Non-integer grid: only the sample axis may be an interpolating
               resample onto a regular in-range grid (sfremap1).  Spatial
               axes still fail closed. */
            if (a != prank - 1) {
                err = "axis " + pax.label +
                      " is not an exact integer selection of the parent "
                      "(RSF o=" + std::to_string(rsf_o) +
                      " d=" + std::to_string(rsf_d) +
                      " n=" + std::to_string(rsf_n) +
                      "; parent o=" + std::to_string(pax.o) +
                      " d=" + std::to_string(pax.d) +
                      " n=" + std::to_string(pax.size) +
                      "). Non-integer geometry on a non-sample axis is not "
                      "supported (fail closed).";
                sel.exact = false;
                return false;
            }
            if (!(rsf_d > 0.0) || rsf_n < 1) {
                err = "resampled sample axis \"" + pax.label +
                      "\" requires d > 0 and n >= 1 (got d=" +
                      std::to_string(rsf_d) + " n=" + std::to_string(rsf_n) +
                      ")";
                sel.exact = false;
                return false;
            }
            const double parent_last = pax.o + (double) (pax.size - 1) * pax.d;
            const double child_last  = rsf_o + (double) (rsf_n - 1) * rsf_d;
            const double atol =
                1e-6 + rtol * (1.0 + std::fabs(parent_last) + std::fabs(pax.o));
            if (rsf_o < pax.o - atol || child_last > parent_last + atol) {
                err = "resampled sample axis \"" + pax.label +
                      "\" extrapolates past the parent range "
                      "(child o=" + std::to_string(rsf_o) +
                      " d=" + std::to_string(rsf_d) +
                      " n=" + std::to_string(rsf_n) +
                      " last=" + std::to_string(child_last) +
                      "; parent o=" + std::to_string(pax.o) +
                      " d=" + std::to_string(pax.d) +
                      " n=" + std::to_string(pax.size) +
                      " last=" + std::to_string(parent_last) +
                      "). No extrapolation past the last parent sample.";
                sel.exact = false;
                return false;
            }
            ax.start = 0;
            ax.step  = 0;
            ax.count = rsf_n;
            ax.kind  = MDIO_SEL_RESAMPLED;
            ax.o = rsf_o;
            ax.d = rsf_d;
            sel.exact = false;
            sel.axes[(size_t) a] = ax;
            continue;
        }

        if (start < 0 || rsf_n < 1 ||
            start + (rsf_n - 1) * step > pax.size - 1) {
            err = "axis " + pax.label + " selection out of parent range "
                  "(start=" + std::to_string(start) +
                  " step=" + std::to_string(step) +
                  " count=" + std::to_string(rsf_n) +
                  " parent_n=" + std::to_string(pax.size) + ")";
            return false;
        }

        ax.start = start;
        ax.step  = step;
        ax.count = rsf_n;
        if (start == 0 && step == 1 && rsf_n == pax.size)
            ax.kind = MDIO_SEL_PRESERVED;
        else if (step == 1)
            ax.kind = MDIO_SEL_SUBSET;
        else
            ax.kind = MDIO_SEL_DECIMATED;
        ax.o = pax.o + (double) start * pax.d;
        ax.d = pax.d * (double) step;
        sel.axes[(size_t) a] = ax;
    }

    mdio_pipe_selection_finalize(sel);
    /* Sample-axis integer decimation maps to "resample" contract name but is
       still an exact integer selection of parent samples (every Nth).  True
       interpolating resample is MDIO_SEL_RESAMPLED above. */
    return true;
}

std::vector<mdio::RangeDescriptor<mdio::Index> >
mdio_pipe_selection_ranges(const MdioSelection& sel)
{
    std::vector<mdio::RangeDescriptor<mdio::Index> > ranges;
    ranges.reserve(sel.axes.size());
    for (size_t i = 0; i < sel.axes.size(); i++) {
        const MdioAxisSel& ax = sel.axes[i];
        mdio::RangeDescriptor<mdio::Index> r = {
            ax.label.c_str(),
            (mdio::Index) ax.start,
            (mdio::Index) (ax.start + ax.count * ax.step),
            (mdio::Index) ax.step};
        ranges.push_back(r);
    }
    return ranges;
}

std::vector<long> mdio_pipe_child_sizes(const MdioSelection& sel)
{
    std::vector<long> out(sel.axes.size());
    for (size_t i = 0; i < sel.axes.size(); i++)
        out[i] = sel.axes[i].count;
    return out;
}

bool mdio_pipe_chunk_shape(const std::string& store_uri,
                           const std::string& datavar,
                           std::vector<long>& chunk,
                           std::string& err)
{
    chunk.clear();
    nlohmann::json meta;

    if (!mdio_pipe_is_remote_uri(store_uri)) {
        fs::path zj = fs::path(store_uri) / datavar / "zarr.json";
        if (!fs::exists(zj)) {
            /* Zarr v2 fallback. */
            zj = fs::path(store_uri) / datavar / ".zarray";
        }
        if (!fs::exists(zj)) {
            err = "no zarr.json / .zarray for \"" + datavar + "\" under " +
                  store_uri;
            return false;
        }
        try {
            std::ifstream in(zj);
            in >> meta;
        } catch (const std::exception& e) {
            err = std::string("parse chunk metadata: ") + e.what();
            return false;
        }
    } else {
        tensorstore::KvStore kvs;
        if (!open_root_kvstore(store_uri, kvs, err)) return false;
        const std::string candidates[] = {
            datavar + "/zarr.json",
            datavar + "/.zarray",
        };
        bool got = false;
        for (const std::string& key : candidates) {
            auto rd = tensorstore::kvstore::Read(kvs, key).result();
            if (!rd.ok() || !rd->has_value()) continue;
            try {
                meta = nlohmann::json::parse(std::string(rd->value));
                got = true;
                break;
            } catch (const std::exception& e) {
                err = std::string("parse ") + key + ": " + e.what();
                return false;
            }
        }
        if (!got) {
            err = "no zarr.json / .zarray for \"" + datavar + "\" under " +
                  store_uri;
            return false;
        }
    }

    const nlohmann::json* arr = NULL;
    if (meta.contains("chunk_grid") && meta["chunk_grid"].is_object() &&
        meta["chunk_grid"].contains("configuration") &&
        meta["chunk_grid"]["configuration"].is_object() &&
        meta["chunk_grid"]["configuration"].contains("chunk_shape") &&
        meta["chunk_grid"]["configuration"]["chunk_shape"].is_array()) {
        arr = &meta["chunk_grid"]["configuration"]["chunk_shape"];
    } else if (meta.contains("chunks") && meta["chunks"].is_array()) {
        arr = &meta["chunks"]; /* Zarr v2 */
    }
    if (!arr) {
        err = datavar + ": no chunk_shape / chunks in metadata";
        return false;
    }
    chunk.reserve(arr->size());
    for (const auto& v : *arr) {
        if (v.is_number_integer())
            chunk.push_back(v.get<long>());
        else if (v.is_number_unsigned())
            chunk.push_back((long) v.get<unsigned long long>());
        else {
            err = datavar + ": non-integer chunk extent";
            chunk.clear();
            return false;
        }
    }
    return true;
}

int mdio_pipe_aligned_panel_traces(int panel_cap, long n_fast,
                                   long chunk_along)
{
    if (n_fast < 1) n_fast = 1;
    if (panel_cap < 1) panel_cap = 64;
    long cap = panel_cap;
    if (cap > n_fast) cap = n_fast;

    if (chunk_along < 1) return (int) cap;

    /* Prefer a multiple of the chunk extent along the streamed axis, not
       exceeding the user cap.  If the cap is smaller than one chunk, use
       one chunk when it fits in n_fast (full-chunk I/O wins over a tiny
       partial panel); otherwise keep the cap. */
    if (cap >= chunk_along) {
        long aligned = (cap / chunk_along) * chunk_along;
        if (aligned < 1) aligned = chunk_along;
        if (aligned > n_fast) aligned = n_fast;
        return (int) aligned;
    }
    if (chunk_along <= n_fast) return (int) chunk_along;
    return (int) cap;
}

bool mdio_pipe_plan_blocks(const std::vector<long>& sizes,
                           const std::vector<long>& chunk,
                           long budget_floats,
                           MdioBlockPlan& plan)
{
    const int rank = (int) sizes.size();
    if (rank < 2 || (int) chunk.size() != rank) return false;
    for (int a = 0; a < rank; a++)
        if (sizes[(size_t) a] < 1 || chunk[(size_t) a] < 1) return false;
    if (budget_floats < 1) budget_floats = 1;

    /* A chunk extent wider than its axis is one truncated chunk, so clamp
       before asking whether a single index covers it. */
    std::vector<long> eff((size_t) rank);
    for (int a = 0; a < rank; a++)
        eff[(size_t) a] = (chunk[(size_t) a] < sizes[(size_t) a]) ?
                          chunk[(size_t) a] : sizes[(size_t) a];

    /* tail[a] = floats in every axis faster than a, taken in full. */
    std::vector<long> tail((size_t) rank, 1);
    for (int a = rank - 2; a >= 0; a--)
        tail[(size_t) a] = tail[(size_t) (a + 1)] * sizes[(size_t) (a + 1)];

    /* Never pivot on the sample axis: a partial trace would break the
       trace-oriented RSF handoff (and dead-trace fill).  Viable pivots are
       the prefix 0..K of the spatial axes, where K is the first axis with
       an effective chunk extent > 1 — beyond it, a single index would only
       ever cover part of a chunk. */
    int last_pivot = rank - 2;
    for (int a = 0; a < rank - 2; a++) {
        if (eff[(size_t) a] > 1) { last_pivot = a; break; }
    }

    int  best = -1;
    long best_block = 0;
    for (int k = 0; k <= last_pivot; k++) {
        long one = eff[(size_t) k] * tail[(size_t) k];
        if (one > budget_floats) continue;
        /* Slowest fit wins: the widest whole-chunk span per Write. */
        long units = budget_floats / one;
        long extent = units * eff[(size_t) k];
        if (extent > sizes[(size_t) k]) extent = sizes[(size_t) k];
        best = k;
        best_block = extent;
        break;
    }

    if (best < 0) {
        /* Nothing fits: fall back to one chunk along the finest viable
           pivot and let the caller decide whether the overshoot is worth
           the alignment. */
        plan.pivot        = last_pivot;
        plan.pivot_block  = eff[(size_t) last_pivot];
        plan.block_floats = plan.pivot_block * tail[(size_t) last_pivot];
        plan.aligned      = plan.block_floats <= budget_floats * 8;
        return true;
    }

    plan.pivot        = best;
    plan.pivot_block  = best_block;
    plan.block_floats = best_block * tail[(size_t) best];
    plan.aligned      = true;
    return true;
}

void mdio_pipe_decode_index(long lin, int ndim, const long* shape, long* idx)
{
    for (int a = ndim - 1; a >= 0; a--) {
        idx[a] = lin % shape[a];
        lin /= shape[a];
    }
}

void mdio_pipe_resolve_block_loop(const std::string& store_uri,
                                  const std::string& datavar,
                                  const std::vector<long>& sizes,
                                  int panel_cap, int blockmb,
                                  MdioBlockLoop& loop,
                                  std::string& note)
{
    const int rank = (int) sizes.size();
    if (panel_cap < 1) panel_cap = 64;
    if (blockmb < 1) blockmb = MDIO_PIPE_BLOCK_MB;

    std::vector<long> chunk;
    std::string cerr;
    if (!mdio_pipe_chunk_shape(store_uri, datavar, chunk, cerr) ||
        (int) chunk.size() != rank) {
        chunk.clear();
        if (!cerr.empty())
            note = "chunk shape unavailable (" + cerr + "); using panel sizing";
    }

    /* Block along the slowest axis whose blocks still span whole chunks, so
       each chunk is touched once instead of once per index of every
       coarser-chunked axis. */
    MdioBlockPlan plan;
    loop.aligned = rank >= 2 && !chunk.empty() &&
        mdio_pipe_plan_blocks(sizes, chunk,
                              (long) blockmb * 1024 * 1024 / 4, plan) &&
        plan.aligned;

    if (loop.aligned) {
        loop.pivot       = plan.pivot;
        loop.pivot_block = plan.pivot_block;
    } else if (rank >= 2) {
        long chunk_along = chunk.empty() ? 0 : chunk[(size_t) (rank - 2)];
        loop.pivot       = rank - 2;
        loop.pivot_block = mdio_pipe_aligned_panel_traces(
            panel_cap, sizes[(size_t) (rank - 2)], chunk_along);
        if (loop.pivot_block < 1) loop.pivot_block = 1;
        if (note.empty())
            note = "no whole-chunk block fits " + std::to_string(blockmb) +
                   " MiB; falling back to unaligned panels of " +
                   std::to_string(loop.pivot_block) +
                   " (slow: expect chunk read-modify-write)";
    } else {
        loop.pivot       = 0;
        loop.pivot_block = sizes[0];    /* single axis: one block */
    }

    /* Axes slower than the pivot are walked one index at a time; axes faster
       than it are always taken in full. */
    loop.n_outer = 1;
    loop.outer_shape.assign((size_t) (loop.pivot > 0 ? loop.pivot : 1), 1);
    for (int a = 0; a < loop.pivot; a++) {
        loop.outer_shape[(size_t) a] = sizes[(size_t) a];
        loop.n_outer *= sizes[(size_t) a];
    }
    loop.n_pivot = sizes[(size_t) loop.pivot];
    loop.tail_floats = 1;
    for (int a = loop.pivot + 1; a < rank; a++)
        loop.tail_floats *= sizes[(size_t) a];
    loop.block_floats = loop.pivot_block * loop.tail_floats;
}

/* Forward decls — used by rematerialize before their definitions. */
static void cap_chunk_shape(nlohmann::json& meta);
static void normalize_structured_dtype(nlohmann::json& meta,
                                       const nlohmann::json* parent_meta);

/* Edit one variable's zarr.json shape + capped chunk_shape to match sel.
   Returns false if the file is missing (ok for some vars) via missing_ok. */
static bool reshape_var_metadata(const fs::path& child_root,
                                 const fs::path& parent_root,
                                 const std::string& name,
                                 const MdioSelection& sel,
                                 bool missing_ok,
                                 std::string& err)
{
    fs::path zj = child_root / name / "zarr.json";
    if (!fs::exists(zj)) {
        if (missing_ok) return true;
        err = "missing zarr.json for \"" + name + "\"";
        return false;
    }
    std::ifstream in(zj);
    nlohmann::json meta;
    try { in >> meta; }
    catch (const std::exception& e) {
        err = std::string("parse ") + zj.string() + ": " + e.what();
        return false;
    }
    in.close();

    /* Restore canonical dtype/codecs from parent before editing shape. */
    fs::path pz = parent_root / name / "zarr.json";
    nlohmann::json parent_meta;
    const nlohmann::json* pmeta = NULL;
    if (fs::exists(pz)) {
        std::ifstream pin(pz);
        if (pin.good()) {
            try { pin >> parent_meta; pmeta = &parent_meta; }
            catch (...) { pmeta = NULL; }
        }
    }
    normalize_structured_dtype(meta, pmeta);
    if (pmeta) {
        if (pmeta->contains("codecs")) meta["codecs"] = (*pmeta)["codecs"];
        if (pmeta->contains("chunk_key_encoding"))
            meta["chunk_key_encoding"] = (*pmeta)["chunk_key_encoding"];
        if (pmeta->contains("fill_value"))
            meta["fill_value"] = (*pmeta)["fill_value"];
        if (pmeta->contains("data_type"))
            meta["data_type"] = (*pmeta)["data_type"];
    }

    if (!meta.contains("shape") || !meta["shape"].is_array()) {
        err = name + ": no shape in zarr.json";
        return false;
    }
    if (!meta.contains("dimension_names") ||
        !meta["dimension_names"].is_array()) {
        err = name + ": no dimension_names in zarr.json";
        return false;
    }

    auto& shape = meta["shape"];
    auto& dims = meta["dimension_names"];
    size_t n = shape.size();
    if (dims.size() < n) n = dims.size();

    bool changed = false;
    for (size_t i = 0; i < n; i++) {
        if (!dims[i].is_string()) continue;
        std::string lab = dims[i].get<std::string>();
        for (size_t a = 0; a < sel.axes.size(); a++) {
            if (sel.axes[a].label == lab &&
                sel.axes[a].kind != MDIO_SEL_PRESERVED) {
                shape[i] = sel.axes[a].count;
                changed = true;
                break;
            }
        }
    }
    if (!changed) return true;

    /* Cap / set chunk_shape to the new shape (single-chunk vars). */
    if (!meta.contains("chunk_grid") || !meta["chunk_grid"].is_object())
        meta["chunk_grid"] = {{"name", "regular"},
                              {"configuration", {{"chunk_shape", shape}}}};
    else {
        meta["chunk_grid"]["name"] = "regular";
        meta["chunk_grid"]["configuration"]["chunk_shape"] = shape;
    }
    /* For multi-dim, keep parent chunk on preserved dims and cap on changed. */
    if (pmeta && pmeta->contains("chunk_grid") &&
        (*pmeta)["chunk_grid"].contains("configuration") &&
        (*pmeta)["chunk_grid"]["configuration"].contains("chunk_shape") &&
        (*pmeta)["chunk_grid"]["configuration"]["chunk_shape"].is_array()) {
        const auto& pc =
            (*pmeta)["chunk_grid"]["configuration"]["chunk_shape"];
        nlohmann::json cs = nlohmann::json::array();
        for (size_t i = 0; i < n; i++) {
            long s = shape[i].get<long>();
            long c = s;
            if (i < pc.size() && pc[i].is_number_integer()) {
                c = pc[i].get<long>();
                if (c > s) c = s;
                if (c < 1) c = s;
            }
            cs.push_back(c);
        }
        meta["chunk_grid"]["configuration"]["chunk_shape"] = cs;
    }

    std::ofstream out(zj, std::ios::trunc);
    if (!out.good()) {
        err = "cannot write " + zj.string();
        return false;
    }
    out << meta.dump(4) << "\n";

    /* Drop stale chunk files so the next Write creates correctly-sized ones. */
    fs::path cdir = child_root / name / "c";
    if (fs::exists(cdir)) {
        std::error_code ec;
        fs::remove_all(cdir, ec);
    }
    return true;
}

/* Does this variable intersect any changed selection axis? */
static bool var_needs_reshape(const nlohmann::json& meta,
                              const MdioSelection& sel)
{
    if (!meta.contains("dimension_names") ||
        !meta["dimension_names"].is_array())
        return false;
    for (auto& d : meta["dimension_names"]) {
        if (!d.is_string()) continue;
        std::string lab = d.get<std::string>();
        for (size_t a = 0; a < sel.axes.size(); a++)
            if (sel.axes[a].label == lab &&
                sel.axes[a].kind != MDIO_SEL_PRESERVED)
                return true;
    }
    return false;
}

/* Write sample coordinate as o + i*d.  Keep an integer parent dtype only when
   every target value is integral; otherwise float64 (never float32) — matches
   reference/resample.py::_new_sample_coord so Axis A regularity passes.
   Used by rematerialize when the sample axis is MDIO_SEL_RESAMPLED. */
static bool write_synthesized_sample_coord(const std::string& child_path,
                                           const std::string& label,
                                           const MdioAxisSel& sax,
                                           const std::string& parent_uri,
                                           std::string& err)
{
    if (sax.count < 1 || !(sax.d > 0.0)) {
        err = "synthesized sample coord \"" + label +
              "\" requires count >= 1 and d > 0";
        return false;
    }

    fs::path zj = fs::path(child_path) / label / "zarr.json";
    if (!fs::exists(zj)) {
        err = "missing sample coord zarr.json at " + zj.string();
        return false;
    }
    nlohmann::json meta;
    try {
        std::ifstream in(zj);
        in >> meta;
    } catch (const std::exception& e) {
        err = std::string("parse sample coord meta: ") + e.what();
        return false;
    }

    /* Parent dtype (for the integral-preserve decision). */
    std::string parent_dtype;
    {
        fs::path pz = fs::path(parent_uri) / label / "zarr.json";
        if (fs::exists(pz)) {
            try {
                std::ifstream pin(pz);
                nlohmann::json pm;
                pin >> pm;
                if (pm.contains("data_type")) {
                    if (pm["data_type"].is_string())
                        parent_dtype = pm["data_type"].get<std::string>();
                    else if (pm["data_type"].is_object())
                        parent_dtype = pm["data_type"].value("name", "");
                }
            } catch (...) {
                parent_dtype.clear();
            }
        }
    }

    bool all_integral = true;
    for (long i = 0; i < sax.count; i++) {
        double v = sax.o + (double) i * sax.d;
        if (std::fabs(v - std::round(v)) > 1e-9) {
            all_integral = false;
            break;
        }
    }
    const bool parent_is_int =
        parent_dtype == "int8" || parent_dtype == "int16" ||
        parent_dtype == "int32" || parent_dtype == "int64" ||
        parent_dtype == "uint8" || parent_dtype == "uint16" ||
        parent_dtype == "uint32" || parent_dtype == "uint64";
    const bool use_int = parent_is_int && all_integral;
    const std::string out_dtype = use_int ? parent_dtype : "float64";

    /* Ensure on-disk dtype matches before opening for Write. */
    if (!meta.contains("data_type") ||
        (meta["data_type"].is_string() &&
         meta["data_type"].get<std::string>() != out_dtype) ||
        (meta["data_type"].is_object() &&
         meta["data_type"].value("name", "") != out_dtype)) {
        meta["data_type"] = out_dtype;
        std::ofstream out(zj, std::ios::trunc);
        if (!out) {
            err = "cannot rewrite sample coord dtype at " + zj.string();
            return false;
        }
        out << meta.dump(4) << "\n";
        std::error_code ec;
        fs::remove_all(fs::path(child_path) / label / "c", ec);
    }

    auto cf = mdio::Dataset::Open(child_path, mdio::constants::kOpen);
    if (!cf.status().ok()) {
        err = "reopen child for synthesized coord: " + cf.status().ToString();
        return false;
    }
    mdio::Dataset child = cf.value();
    auto cvr = child.variables.at(label);
    if (!cvr.status().ok()) {
        err = "open synthesized coord \"" + label + "\": " +
              cvr.status().ToString();
        return false;
    }
    auto cvar = cvr.value();
    auto store = cvar.get_store();

    if (use_int) {
        /* Write through the variable's native integer dtype via from_variable
           when possible; fall back to float64 path if reopen typed wrong. */
        auto try_i32 = child.variables.get<mdio::dtypes::int32_t>(label);
        if (try_i32.status().ok() &&
            (out_dtype == "int32" || parent_dtype == "int32")) {
            auto var = try_i32.value();
            auto vdr = mdio::from_variable<mdio::dtypes::int32_t>(var);
            if (!vdr.status().ok()) {
                err = "from_variable int32 coord: " + vdr.status().ToString();
                return false;
            }
            auto vd = vdr.value();
            long n = (long) vd.num_samples();
            if (n != sax.count) {
                err = "synthesized coord length mismatch";
                return false;
            }
            auto off = vd.get_flattened_offset();
            auto* p = static_cast<mdio::dtypes::int32_t*>(
                vd.get_data_accessor().data());
            for (long i = 0; i < n; i++)
                p[off + i] = (mdio::dtypes::int32_t) std::llround(
                    sax.o + (double) i * sax.d);
            if (!var.Write(vd).status().ok()) {
                err = "Write synthesized int32 sample coord failed";
                return false;
            }
            return true;
        }
    }

    /* float64 (or non-int32 integral — still write float64 for simplicity of
       the non-int32 integral edge cases; Axis A accepts either). */
    {
        /* Force float64 store via tensorstore after dtype rewrite above. */
        auto rd = tensorstore::Read(store).result();
        /* Store may be empty (no chunks) — allocate and Write instead. */
        std::vector<double> vals((size_t) sax.count);
        for (long i = 0; i < sax.count; i++)
            vals[(size_t) i] = sax.o + (double) i * sax.d;

        /* Prefer typed Variable Write when dtype is float64. */
        auto try_f64 = child.variables.get<mdio::dtypes::float64_t>(label);
        if (try_f64.status().ok()) {
            auto var = try_f64.value();
            auto vdr = mdio::from_variable<mdio::dtypes::float64_t>(var);
            if (!vdr.status().ok()) {
                err = "from_variable float64 coord: " + vdr.status().ToString();
                return false;
            }
            auto vd = vdr.value();
            long n = (long) vd.num_samples();
            if (n != sax.count) {
                err = "synthesized float64 coord length " + std::to_string(n) +
                      " != " + std::to_string(sax.count);
                return false;
            }
            auto off = vd.get_flattened_offset();
            auto* p = static_cast<double*>(vd.get_data_accessor().data());
            for (long i = 0; i < n; i++) p[off + i] = vals[(size_t) i];
            if (!var.Write(vd).status().ok()) {
                err = "Write synthesized float64 sample coord failed";
                return false;
            }
            return true;
        }

        /* Last resort: int32 parent kept integral. */
        auto try_i32 = child.variables.get<mdio::dtypes::int32_t>(label);
        if (try_i32.status().ok() && use_int) {
            auto var = try_i32.value();
            auto vdr = mdio::from_variable<mdio::dtypes::int32_t>(var);
            if (!vdr.status().ok()) {
                err = vdr.status().ToString();
                return false;
            }
            auto vd = vdr.value();
            long n = (long) vd.num_samples();
            auto off = vd.get_flattened_offset();
            auto* p = static_cast<mdio::dtypes::int32_t*>(
                vd.get_data_accessor().data());
            for (long i = 0; i < n && i < sax.count; i++)
                p[off + i] = (mdio::dtypes::int32_t) std::llround(vals[(size_t) i]);
            if (!var.Write(vd).status().ok()) {
                err = "Write synthesized int sample coord failed";
                return false;
            }
            return true;
        }

        err = "cannot open sample coord \"" + label +
              "\" as float64 or int32 after dtype=" + out_dtype +
              " rewrite (store dtype probe failed)";
        (void) rd;
        (void) store;
        return false;
    }
}

static bool rematerialize_child_variables(mdio::Dataset& parent,
                                          const std::string& parent_uri,
                                          const std::string& child_path,
                                          const MdioSelection& sel,
                                          const std::string& datavar,
                                          std::string& err)
{
    fs::path child_root(child_path);
    fs::path parent_root(parent_uri);
    if (parent_uri.find("://") != std::string::npos) {
        err = "reshaped child build requires a local parent path";
        return false;
    }

    /* First pass: reshape metadata + wipe chunks for every affected var
       (including amplitude). */
    std::vector<std::string> to_copy;
    std::error_code ec;
    for (auto& ent : fs::directory_iterator(child_root, ec)) {
        if (ec || !ent.is_directory()) continue;
        std::string name = ent.path().filename().string();
        fs::path zj = ent.path() / "zarr.json";
        if (!fs::exists(zj)) continue;
        std::ifstream in(zj);
        nlohmann::json meta;
        try { in >> meta; } catch (...) { continue; }
        if (!var_needs_reshape(meta, sel)) continue;
        if (!reshape_var_metadata(child_root, parent_root, name, sel,
                                  false, err))
            return false;
        if (name != datavar) to_copy.push_back(name);
    }

    /* Reopen child with new shapes and Copy sliced parent data. */
    auto cf = mdio::Dataset::Open(child_path, mdio::constants::kOpen);
    if (!cf.status().ok()) {
        err = "reopen after reshape: " + cf.status().ToString();
        return false;
    }
    mdio::Dataset child = cf.value();

    const int prank = (int) sel.axes.size();
    const MdioAxisSel* sample_ax =
        (prank >= 1) ? &sel.axes[(size_t) (prank - 1)] : NULL;

    for (const std::string& name : to_copy) {
        auto pvr = parent.variables.at(name);
        auto cvr = child.variables.at(name);
        if (!pvr.status().ok() || !cvr.status().ok()) {
            err = "open \"" + name + "\" for rematerialize";
            return false;
        }
        auto pvar = pvr.value();
        auto cvar = cvr.value();

        auto labels = pvar.dimensions().labels();
        bool touches_resampled = false;
        std::vector<mdio::RangeDescriptor<mdio::Index> > src_ranges;
        for (size_t i = 0; i < labels.size(); i++) {
            std::string lab(labels[i].begin(), labels[i].end());
            if (lab.empty()) continue;
            const MdioAxisSel* ax = NULL;
            for (size_t a = 0; a < sel.axes.size(); a++)
                if (sel.axes[a].label == lab) { ax = &sel.axes[a]; break; }
            if (!ax) continue;
            if (ax->kind == MDIO_SEL_RESAMPLED) {
                touches_resampled = true;
                break;
            }
            mdio::RangeDescriptor<mdio::Index> src = {
                ax->label.c_str(),
                (mdio::Index) ax->start,
                (mdio::Index) (ax->start + ax->count * ax->step),
                (mdio::Index) ax->step};
            src_ranges.push_back(src);
        }

        if (touches_resampled) {
            /* Only the sample-coordinate variable is rewritten across a
               resampled sample axis; anything else with that dim fails. */
            if (sample_ax && name == sample_ax->label &&
                sample_ax->kind == MDIO_SEL_RESAMPLED) {
                if (!write_synthesized_sample_coord(child_path, name,
                                                    *sample_ax, parent_uri,
                                                    err))
                    return false;
                continue;
            }
            err = "cannot rematerialize \"" + name +
                  "\" across a resampled sample axis (only the sample "
                  "coordinate is synthesized)";
            return false;
        }

        auto src_sl = pvar.slice(src_ranges);
        if (!src_sl.status().ok()) {
            err = "slice parent \"" + name + "\": " +
                  src_sl.status().ToString();
            return false;
        }
        auto copy = tensorstore::Copy(src_sl.value().get_store(),
                                      cvar.get_store());
        if (!copy.status().ok()) {
            err = "rematerialize Copy \"" + name + "\": " +
                  copy.status().ToString();
            return false;
        }
    }

    /* Metadata was rewritten on disk in reshape_var_metadata; Copy writes
       chunk payloads only.  CommitMetadata would fail with "No variables
       were modified" here — skip it. */
    (void) child;
    return true;
}

/* Cap chunk_shape[i] to shape[i].  TensorStore Resize updates shape but
   leaves parent chunk sizes, which can exceed the new shape. */
static void cap_chunk_shape(nlohmann::json& meta)
{
    if (!meta.contains("shape") || !meta["shape"].is_array()) return;
    nlohmann::json* chunk = NULL;
    if (meta.contains("chunk_grid") && meta["chunk_grid"].is_object() &&
        meta["chunk_grid"].contains("configuration") &&
        meta["chunk_grid"]["configuration"].is_object() &&
        meta["chunk_grid"]["configuration"].contains("chunk_shape") &&
        meta["chunk_grid"]["configuration"]["chunk_shape"].is_array())
        chunk = &meta["chunk_grid"]["configuration"]["chunk_shape"];
    if (!chunk) return;
    const auto& shape = meta["shape"];
    size_t n = shape.size();
    if (chunk->size() < n) n = chunk->size();
    for (size_t i = 0; i < n; i++) {
        if (!(*chunk)[i].is_number_integer() || !shape[i].is_number_integer())
            continue;
        long c = (*chunk)[i].get<long>();
        long s = shape[i].get<long>();
        if (c > s && s > 0) (*chunk)[i] = s;
    }
}

/* Convert TensorStore's rewritten {"name":"struct", fields:[{name,data_type}]}
   back to the MDIO/Python Zarr-v3 form {"name":"structured",
   fields:[[name,dtype],...]} that multidimio accepts. */
static void normalize_structured_dtype(nlohmann::json& meta,
                                       const nlohmann::json* parent_meta)
{
    if (!meta.contains("data_type") || !meta["data_type"].is_object()) return;
    auto& dt = meta["data_type"];
    std::string name = dt.value("name", "");
    if (name != "struct" && name != "structured") return;

    /* Prefer the parent's canonical dtype / fill when available.  Do NOT
       copy codecs — the child may have rewritten chunk payloads under a
       different codec (e.g. uncompressed bytes after sample-geometry stamp). */
    if (parent_meta && parent_meta->contains("data_type")) {
        meta["data_type"] = (*parent_meta)["data_type"];
        if (parent_meta->contains("fill_value"))
            meta["fill_value"] = (*parent_meta)["fill_value"];
        if (parent_meta->contains("chunk_key_encoding") &&
            !meta.contains("chunk_key_encoding"))
            meta["chunk_key_encoding"] = (*parent_meta)["chunk_key_encoding"];
        return;
    }

    if (name != "struct") return;
    if (!dt.contains("configuration") ||
        !dt["configuration"].contains("fields") ||
        !dt["configuration"]["fields"].is_array())
        return;

    nlohmann::json fields = nlohmann::json::array();
    for (auto& f : dt["configuration"]["fields"]) {
        if (f.is_array() && f.size() >= 2) {
            fields.push_back(f);
        } else if (f.is_object() && f.contains("name") &&
                   f.contains("data_type")) {
            fields.push_back(nlohmann::json::array({f["name"], f["data_type"]}));
        }
    }
    meta["data_type"] = {
        {"name", "structured"},
        {"configuration", {{"fields", fields}}}
    };

    /* Object-form fill_value is not accepted by Python MDIO; use base64
       zeros sized to the struct (sum of known scalar widths). */
    if (meta.contains("fill_value") && meta["fill_value"].is_object()) {
        long nbytes = 0;
        for (auto& f : fields) {
            if (!f.is_array() || f.size() < 2) continue;
            std::string t = f[1].get<std::string>();
            if (t == "int8" || t == "uint8" || t == "bool") nbytes += 1;
            else if (t == "int16" || t == "uint16") nbytes += 2;
            else if (t == "int32" || t == "uint32" || t == "float32")
                nbytes += 4;
            else if (t == "int64" || t == "uint64" || t == "float64")
                nbytes += 8;
            else nbytes += 4;
        }
        static const char* b64 =
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        std::string raw((size_t) nbytes, '\0');
        std::string enc;
        enc.reserve(((size_t) nbytes + 2) / 3 * 4);
        for (size_t i = 0; i < raw.size(); i += 3) {
            unsigned n = (unsigned) (unsigned char) raw[i] << 16;
            if (i + 1 < raw.size())
                n |= (unsigned) (unsigned char) raw[i + 1] << 8;
            if (i + 2 < raw.size())
                n |= (unsigned) (unsigned char) raw[i + 2];
            enc.push_back(b64[(n >> 18) & 63]);
            enc.push_back(b64[(n >> 12) & 63]);
            enc.push_back(i + 1 < raw.size() ? b64[(n >> 6) & 63] : '=');
            enc.push_back(i + 2 < raw.size() ? b64[n & 63] : '=');
        }
        meta["fill_value"] = enc;
    }
}

/* After Resize / SelectField Write, TensorStore may rewrite zarr.json into a
   form Python MDIO cannot open (struct vs structured) and leave chunk_shape
   larger than shape.  Repair every variable directory in-place. */
bool mdio_pipe_repair_zarr_metadata(const std::string& child_path,
                                    const std::string& parent_uri,
                                    std::string& err,
                                    bool cap_chunks)
{
    fs::path child_root(child_path);
    fs::path parent_root;
    bool have_parent = false;
    if (!parent_uri.empty() && parent_uri.find("://") == std::string::npos) {
        parent_root = fs::path(parent_uri);
        have_parent = fs::is_directory(parent_root);
    }

    std::error_code ec;
    if (!fs::is_directory(child_root, ec)) {
        err = "repair: not a directory: " + child_path;
        return false;
    }

    for (auto& ent : fs::directory_iterator(child_root, ec)) {
        if (ec) break;
        if (!ent.is_directory()) continue;
        fs::path zj = ent.path() / "zarr.json";
        if (!fs::exists(zj)) continue;

        std::ifstream in(zj);
        if (!in.good()) {
            err = "repair: cannot read " + zj.string();
            return false;
        }
        nlohmann::json meta;
        try {
            in >> meta;
        } catch (const std::exception& e) {
            err = std::string("repair: parse ") + zj.string() + ": " + e.what();
            return false;
        }
        in.close();

        nlohmann::json parent_meta;
        const nlohmann::json* pmeta = NULL;
        if (have_parent) {
            fs::path pz = parent_root / ent.path().filename() / "zarr.json";
            if (fs::exists(pz)) {
                std::ifstream pin(pz);
                if (pin.good()) {
                    try {
                        pin >> parent_meta;
                        pmeta = &parent_meta;
                    } catch (...) {
                        pmeta = NULL;
                    }
                }
            }
        }

        normalize_structured_dtype(meta, pmeta);
        /* Cap only after all TensorStore writes — editing chunk_shape while
           the store is still open for panel writes breaks subsequent Writes. */
        if (cap_chunks) cap_chunk_shape(meta);

        /* Ensure chunk_key_encoding keeps the default separator. */
        if (meta.contains("chunk_key_encoding") &&
            meta["chunk_key_encoding"].is_object() &&
            !meta["chunk_key_encoding"].contains("configuration")) {
            meta["chunk_key_encoding"]["configuration"] = {
                {"separator", "/"}};
        }

        std::ofstream out(zj, std::ios::trunc);
        if (!out.good()) {
            err = "repair: cannot write " + zj.string();
            return false;
        }
        out << meta.dump(4) << "\n";
    }
    return true;
}

bool mdio_pipe_build_child_store(mdio::Dataset& parent,
                                 const std::string& parent_uri,
                                 const MdioSelection& sel,
                                 const std::string& datavar,
                                 const std::string& tmp_path,
                                 std::string& err)
{
    if (!mdio_clone_store(parent_uri, tmp_path, err)) return false;

    if (sel.same_geometry) return true;

    /* Reshape metadata + rewrite correctly-sized chunks from parent.isel.
       (Resize alone leaves oversized single-chunk payloads that Python
       MDIO cannot open.)  Leave dtype as TensorStore wrote it — the
       writer repairs to Python-MDIO "structured" just before publish. */
    if (!rematerialize_child_variables(parent, parent_uri, tmp_path, sel,
                                       datavar, err)) {
        mdio_pipe_remove_path(tmp_path);
        return false;
    }
    return true;
}

/* The sample-axis geometry a windowed child has to restamp into its SEG-Y
   metadata, carried in SEG-Y's own units rather than the coordinate's. */
struct SampleGeometryStamp {
    long ns;        /* samples per trace after windowing           */
    int  dt_us;     /* sample interval, microseconds               */
    int  shift_ms;  /* added to delay_recording_time, milliseconds */
};

/* Byte offset of every field within one structured-header record, parsed from
   the zarr.json `data_type` field list; `itemsize` returns the record stride.

   Zarr v3 structured dtypes are packed without padding, so a field's offset is
   just the running sum of the widths ahead of it.  That makes the width table
   below load-bearing: a single wrong width silently shifts every later field,
   so an unrecognized type fails closed instead of misplacing the stamp. */
static bool structured_field_offsets(const nlohmann::json& meta,
                                     std::unordered_map<std::string, size_t>& off,
                                     size_t& itemsize,
                                     std::string& err)
{
    off.clear();
    itemsize = 0;
    if (!meta.contains("data_type") || !meta["data_type"].is_object()) {
        err = "headers: no data_type";
        return false;
    }
    const auto& dt = meta["data_type"];
    if (!dt.contains("configuration") ||
        !dt["configuration"].contains("fields") ||
        !dt["configuration"]["fields"].is_array()) {
        err = "headers: no fields";
        return false;
    }

    size_t cursor = 0;
    for (auto& f : dt["configuration"]["fields"]) {
        /* Two spellings occur in the wild: [name, type] pairs and objects. */
        std::string name, typ;
        if (f.is_array() && f.size() >= 2) {
            name = f[0].get<std::string>();
            typ  = f[1].get<std::string>();
        } else if (f.is_object()) {
            name = f.value("name", "");
            typ  = f.value("data_type", "");
        } else continue;

        size_t w;
        if      (typ == "int8"  || typ == "uint8"  || typ == "bool")    w = 1;
        else if (typ == "int16" || typ == "uint16")                     w = 2;
        else if (typ == "int32" || typ == "uint32" || typ == "float32") w = 4;
        else if (typ == "int64" || typ == "uint64" || typ == "float64") w = 8;
        else {
            err = "headers: unsupported field type \"" + typ + "\" for \"" +
                  name + "\" (would misplace every later field)";
            return false;
        }
        off[name] = cursor;
        cursor += w;
    }

    itemsize = cursor;
    return itemsize > 0;
}

/* Patch samples_per_trace / sample_interval / delay_recording_time through
   the structured headers store's open_as_void trailing byte dimension.

   mdio-cpp opens structured arrays as rank-N+1 byte views.  Reading and writing
   that view lets TensorStore keep dtype, codec, chunk grid, chunk-key encoding
   and edge-chunk handling — avoiding SelectField writes (which zero siblings)
   and the old single-chunk hand-written Zarr path (Bug A/B). */
static bool stamp_headers_sample_geometry(const std::string& child_path,
                                          const std::string& /*parent_uri*/,
                                          const MdioSelection& /*sel*/,
                                          const SampleGeometryStamp& geom,
                                          const std::vector<char>& live,
                                          std::string& err)
{
    const fs::path child_hdr = fs::path(child_path) / "headers";
    const fs::path zj = child_hdr / "zarr.json";
    if (!fs::exists(zj)) return true;  /* no headers */

    nlohmann::json meta;
    try {
        std::ifstream in(zj);
        in >> meta;
    } catch (const std::exception& e) {
        err = std::string("headers zarr.json: ") + e.what();
        return false;
    }

    std::unordered_map<std::string, size_t> foff;
    size_t itemsize = 0;
    if (!structured_field_offsets(meta, foff, itemsize, err)) return false;

    auto need = [&](const char* name) -> bool {
        if (foff.find(name) != foff.end()) return true;
        err = std::string("headers missing required field \"") + name +
              "\" for sample-geometry stamp at " + child_hdr.string();
        return false;
    };
    if (!need("samples_per_trace") || !need("sample_interval") ||
        !need("delay_recording_time"))
        return false;

    /* Field widths: SEG-Y sample-geometry fields are int16. */
    auto width_of = [&](const char* name) -> size_t {
        size_t o = foff[name];
        size_t next = itemsize;
        for (const auto& kv : foff) {
            if (kv.second > o && kv.second < next) next = kv.second;
        }
        return next - o;
    };
    if (width_of("samples_per_trace") != 2 ||
        width_of("sample_interval") != 2 ||
        width_of("delay_recording_time") != 2) {
        err = "sample-geometry stamp expects int16 fields; widths are "
              "samples_per_trace=" + std::to_string(width_of("samples_per_trace")) +
              " sample_interval=" + std::to_string(width_of("sample_interval")) +
              " delay_recording_time=" +
              std::to_string(width_of("delay_recording_time")) +
              " at " + child_hdr.string();
        return false;
    }

    auto cf = mdio::Dataset::Open(child_path, mdio::constants::kOpen);
    if (!cf.status().ok()) {
        err = "reopen child for header byte-view stamp: " +
              cf.status().ToString();
        return false;
    }
    auto cvr = cf.value().variables.at("headers");
    if (!cvr.status().ok()) {
        err = "open headers for byte-view stamp: " + cvr.status().ToString();
        return false;
    }
    auto hdr = cvr.value();
    auto store = hdr.get_store();

    /* Expect open_as_void: dtype byte/uint8, trailing unlabeled byte axis. */
    {
        std::string dt(store.dtype().name());
        if (dt != "byte" && dt != "uint8") {
            err = "headers store is not an open_as_void byte view (dtype=" +
                  dt + "); cannot patch sample-geometry fields in place";
            return false;
        }
    }
    const int rank = (int) store.rank();
    if (rank < 2) {
        err = "headers byte-view rank " + std::to_string(rank) +
              " < 2 at " + child_hdr.string();
        return false;
    }
    auto shape = store.domain().shape();
    const long byte_dim = (long) shape[rank - 1];
    if (byte_dim != (long) itemsize) {
        err = "headers trailing byte dim " + std::to_string(byte_dim) +
              " != structured itemsize " + std::to_string(itemsize) +
              " at " + child_hdr.string();
        return false;
    }
    long ntr = 1;
    for (int i = 0; i < rank - 1; i++) {
        long s = (long) shape[i];
        if (s < 1) {
            err = "headers spatial shape has non-positive entry";
            return false;
        }
        if (ntr > (std::numeric_limits<long>::max() / s)) {
            err = "headers spatial shape product overflows";
            return false;
        }
        ntr *= s;
    }
    if (ntr < 1) {
        err = "headers byte-view has zero traces at " + child_hdr.string();
        return false;
    }

    auto rd = tensorstore::Read(store).result();
    if (!rd.ok()) {
        err = "tensorstore::Read headers byte view failed: " +
              rd.status().ToString();
        return false;
    }
    auto arr = rd.value();
    if ((long) arr.num_elements() != ntr * (long) itemsize) {
        err = "headers byte-view element count " +
              std::to_string((long) arr.num_elements()) + " != ntr*itemsize " +
              std::to_string(ntr * (long) itemsize);
        return false;
    }

    unsigned char* base = static_cast<unsigned char*>(arr.data());
    const size_t off_ns  = foff["samples_per_trace"];
    const size_t off_dt  = foff["sample_interval"];
    const size_t off_del = foff["delay_recording_time"];
    const int16_t ns16 = (int16_t) geom.ns;
    const int16_t dt16 = (int16_t) geom.dt_us;

    for (long i = 0; i < ntr; i++) {
        unsigned char* rec = base + (size_t) i * itemsize;
        std::memcpy(rec + off_ns, &ns16, 2);
        std::memcpy(rec + off_dt, &dt16, 2);

        int16_t delay = 0;
        std::memcpy(&delay, rec + off_del, 2);
        const bool is_live = (size_t) i >= live.size() || live[(size_t) i];
        if (is_live) {
            long dv = (long) delay + (long) geom.shift_ms;
            if (dv < -32768) dv = -32768;
            if (dv >  32767) dv =  32767;
            delay = (int16_t) dv;
            std::memcpy(rec + off_del, &delay, 2);
        }
    }

    auto wr = tensorstore::Write(arr, store).result();
    if (!wr.ok()) {
        err = "tensorstore::Write headers byte view failed: " +
              wr.status().ToString();
        return false;
    }
    return true;
}

/* Patch segy_file_header.attributes.binaryHeader in zarr.json directly, rather
   than through CommitMetadata, which has side effects on sibling variables.
   A parent without that attribute block is left alone. */
static bool stamp_binary_header_ns_dt(const std::string& child_path,
                                      const SampleGeometryStamp& geom,
                                      std::string& err)
{
    fs::path zj = fs::path(child_path) / "segy_file_header" / "zarr.json";
    if (!fs::exists(zj)) return true;
    std::ifstream in(zj);
    nlohmann::json meta;
    try { in >> meta; }
    catch (const std::exception& e) {
        err = std::string("segy_file_header zarr.json: ") + e.what();
        return false;
    }
    in.close();

    nlohmann::json* attrs = NULL;
    if (meta.contains("attributes") && meta["attributes"].is_object())
        attrs = &meta["attributes"];
    if (!attrs) return true;
    nlohmann::json* binary = NULL;
    if (attrs->contains("binaryHeader") &&
        (*attrs)["binaryHeader"].is_object())
        binary = &(*attrs)["binaryHeader"];
    if (!binary) return true;
    (*binary)["samples_per_trace"] = (int) geom.ns;
    (*binary)["sample_interval"]   = geom.dt_us;

    std::ofstream out(zj, std::ios::trunc);
    if (!out.good()) {
        err = "cannot write " + zj.string();
        return false;
    }
    out << meta.dump(4) << "\n";
    return true;
}

bool mdio_pipe_stamp_sample_geometry(mdio::Dataset& child,
                                     const std::string& child_path,
                                     const std::string& parent_uri,
                                     const std::string& datavar,
                                     const std::vector<MdioAxis>& parent_axes,
                                     const MdioSelection& sel,
                                     std::string& err)
{
    if (!sel.sample_changed) return true;
    const int prank = (int) sel.axes.size();
    if (prank < 1) return true;
    const MdioAxisSel& sax = sel.axes[(size_t) (prank - 1)];
    const MdioAxis& pax = parent_axes[(size_t) (prank - 1)];

    /* Derive SEG-Y sample geometry from the child grid (o/d), not from integer
       start/step.  Matches reference/_sampling.py::stamp_sample_geometry /
       shift_recording_delay.  For an integer selection, sax.o/d were set to
       pax.o+start*pax.d and pax.d*step, so this is numerically identical. */
    SampleGeometryStamp geom;
    geom.ns       = sax.count;
    geom.dt_us    = (int) std::llround(sax.d * 1000.0);
    geom.shift_ms = (int) std::llround(sax.o - pax.o);

    if (!stamp_binary_header_ns_dt(child_path, geom, err))
        return false;

    /* Dead traces are stamped differently, so the mask is required, not
       optional — a missing one would silently shift dead traces. */
    std::vector<char> live;
    if (!mdio_pipe_read_trace_mask(child, live)) {
        err = "failed to read child trace_mask for sample-geometry stamp at " +
              child_path + " (datavar=" + datavar +
              "); refusing to guess live/dead because a missing mask would "
              "silently shift dead-trace delays";
        return false;
    }

    if (!stamp_headers_sample_geometry(child_path, parent_uri, sel,
                                       geom, live, err))
        return false;

    /* TensorStore may rewrite headers dtype spelling; repair back to
       Python-MDIO form (no chunk cap — shapes unchanged). */
    if (!mdio_pipe_repair_zarr_metadata(child_path, parent_uri, err,
                                        /*cap_chunks=*/false))
        return false;

    (void) datavar;
    (void) child;
    return true;
}

