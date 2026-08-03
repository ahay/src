/* Pipe helpers for sfmdioread / sfmdiowrite. */

#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unistd.h>
#include <vector>

#include "absl/strings/cord.h"
#include "mdio_pipe.hh"
#include "tensorstore/array.h"
#include "tensorstore/kvstore/key_range.h"
#include "tensorstore/kvstore/kvstore.h"
#include "tensorstore/kvstore/operations.h"
#include "tensorstore/tensorstore.h"

namespace {

/* 128-bit FNV-1a parent descriptor digest (pipe mismatch check, not crypto). */
static std::string digest_hex(const std::string& msg)
{
    unsigned long long h[2] = {0xcbf29ce484222325ULL, 0x9ae16a3b2f90404fULL};
    for (size_t i = 0; i < msg.size(); i++)
        for (int l = 0; l < 2; l++) {
            h[l] ^= (unsigned char) msg[i];
            h[l] *= 0x100000001b3ULL;
        }
    char hex[33];
    snprintf(hex, sizeof(hex), "%016llx%016llx", h[0], h[1]);
    return std::string(hex);
}

/* Dataset name from mdio-cpp metadata. */
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

/* Canonical descriptor — must match reader and writer. */
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

std::string mdio_pipe_fingerprint(mdio::Dataset& ds,
                                  const std::string& datavar,
                                  const std::string& dtype)
{
    return "mdio1:" + digest_hex(canonical_descriptor(ds, datavar, dtype));
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

static void read_hist_key(sf_file in, const char* key, std::string& out)
{
    out.clear();
    char* s = sf_histstring(in, key);
    if (s) { out = s; free(s); }
}

void mdio_pipe_read_context(sf_file in, MdioPipeContext& ctx)
{
    read_hist_key(in, MDIO_PIPE_KEY_SOURCE,      ctx.source);
    read_hist_key(in, MDIO_PIPE_KEY_FINGERPRINT, ctx.fingerprint);
    read_hist_key(in, MDIO_PIPE_KEY_DATA_VAR,    ctx.data_variable);
    read_hist_key(in, MDIO_PIPE_KEY_DTYPE,       ctx.dtype);
    read_hist_key(in, MDIO_PIPE_KEY_CONTRACT,    ctx.contract);
    read_hist_key(in, MDIO_PIPE_KEY_VERSION,     ctx.version);
    read_hist_key(in, MDIO_PIPE_KEY_SELECTION,   ctx.selection);
    ctx.present = !ctx.source.empty();
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

/* Stats zero tolerance (mdio_compat STATS_ZERO_ATOL). */
#define MDIO_PIPE_STATS_ZERO_ATOL 1e-8

/* Open root kvstore for MDIO URI (local path or s3:// / gs://). */
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

/* Join kvstore key segments (no leading slash). */
static std::string join_key(const std::string& a, const std::string& b)
{
    if (a.empty()) return b;
    if (b.empty()) return a;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

/* Delete all keys under dataset URI. */
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

/* Delete keys with the given prefix (e.g. "amplitude/c/"). */
static bool kv_delete_prefix(tensorstore::KvStore& kvs,
                             const std::string& prefix,
                             std::string& err)
{
    auto del = tensorstore::kvstore::DeleteRange(
                   kvs, tensorstore::KeyRange::Prefix(prefix))
                   .result();
    if (!del.ok()) {
        err = "kvstore DeleteRange(\"" + prefix +
              "\") failed: " + del.status().ToString();
        return false;
    }
    return true;
}

/* Read JSON from a kvstore key. found=false when the key is missing. */
static bool kv_read_json(tensorstore::KvStore& kvs, const std::string& key,
                         nlohmann::json& out, bool& found, std::string& err)
{
    found = false;
    auto rd = tensorstore::kvstore::Read(kvs, key).result();
    if (!rd.ok()) {
        err = "kvstore read \"" + key + "\": " + rd.status().ToString();
        return false;
    }
    if (!rd->has_value()) return true;
    try {
        out = nlohmann::json::parse(std::string(rd->value));
    } catch (const std::exception& e) {
        err = "parse \"" + key + "\": " + e.what();
        return false;
    }
    found = true;
    return true;
}

/* Write JSON to a kvstore key (pretty-printed). */
static bool kv_write_json(tensorstore::KvStore& kvs, const std::string& key,
                          const nlohmann::json& meta, std::string& err)
{
    try {
        std::string text = meta.dump(4) + "\n";
        auto wr =
            tensorstore::kvstore::Write(kvs, key, absl::Cord(text)).result();
        if (!wr.ok()) {
            err = "kvstore write \"" + key + "\": " + wr.status().ToString();
            return false;
        }
    } catch (const std::exception& e) {
        err = "serialize \"" + key + "\": " + e.what();
        return false;
    }
    return true;
}

/* Variable zarr.json or .zarray metadata. */
static bool kv_read_var_meta(tensorstore::KvStore& kvs, const std::string& var,
                             nlohmann::json& out, std::string& err)
{
    bool found = false;
    for (const char* leaf : {"zarr.json", ".zarray"}) {
        if (!kv_read_json(kvs, join_key(var, leaf), out, found, err))
            return false;
        if (found) return true;
    }
    err = "no zarr.json / .zarray for \"" + var + "\"";
    return false;
}

/* Write variable metadata as zarr.json (v3 layout used by the pipe). */
static bool kv_write_var_meta(tensorstore::KvStore& kvs, const std::string& var,
                              const nlohmann::json& meta, std::string& err)
{
    return kv_write_json(kvs, join_key(var, "zarr.json"), meta, err);
}

static bool is_zarr_var_metadata_name(const std::string& name)
{
    return name == "zarr.json" || name == ".zarray" || name == ".zattrs" ||
           name == ".zgroup" || name == ".zmetadata";
}

/* skip_var chunk payload key (not zarr.json etc.). */
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
    /* nested under skip_var/ -> chunk payload */
    return true;
}

/* Copy kvstore keys; optionally skip skip_var chunk payloads. */
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

bool mdio_pipe_is_remote_uri(const std::string& uri)
{
    const std::string::size_type p = uri.find("://");
    if (p == std::string::npos) return false;       /* bare local path */
    return uri.compare(0, 7, "file://") != 0;        /* file:// is local */
}

void mdio_pipe_remove_path(const std::string& path)
{
    std::string err;
    kvstore_delete_all(path, err); /* best-effort */
}

bool mdio_clone_store_skip_var(const std::string& parent_uri,
                               const std::string& child_path,
                               const std::string& skip_var,
                               std::string& err)
{
    mdio_pipe_remove_path(child_path);

    tensorstore::KvStore src, dst;
    if (!open_root_kvstore(parent_uri, src, err)) return false;
    if (!open_root_kvstore(child_path, dst, err)) return false;
    if (!kvstore_copy_all(src, dst, err, skip_var)) {
        mdio_pipe_remove_path(child_path);
        return false;
    }
    return true;
}

static bool mdio_clone_store(const std::string& parent_uri,
                             const std::string& child_path,
                             std::string& err)
{
    return mdio_clone_store_skip_var(parent_uri, child_path,
                                     /*skip_var=*/std::string(), err);
}

bool mdio_pipe_materialize_local(const std::string& uri,
                                 std::string& local,
                                 std::string& err)
{
    local.clear();
    if (!mdio_pipe_is_remote_uri(uri)) {
        local = uri;
        return true;
    }

    /* kvstore open already splits s3://bucket/... correctly; Dataset::Open
       on the remote URI does not (Zarr-V3 per-variable specs drop bucket).
       Stage a full local clone so Dataset::Open sees a file store. */
    const char* base = getenv("TMPDIR");
    std::string dir = (NULL != base && '\0' != *base) ? base : "/tmp";
    if (!dir.empty() && '/' == dir.back()) dir.pop_back();

    static unsigned seq = 0;
    local = dir + "/sfmdio_open." + std::to_string((long) getpid()) + "." +
            std::to_string(++seq) + ".mdio";

    if (!mdio_clone_store(uri, local, err)) {
        mdio_pipe_remove_path(local);
        local.clear();
        return false;
    }
    return true;
}

bool mdio_pipe_publish(const std::string& tmp_path,
                       const std::string& final_path,
                       std::string& err)
{
    /* Copy tmp->final then delete tmp. Object stores have no atomic dir
       rename; local uses the same path for uniformity. */
    mdio_pipe_remove_path(final_path);
    tensorstore::KvStore src, dst;
    if (!open_root_kvstore(tmp_path, src, err)) return false;
    if (!open_root_kvstore(final_path, dst, err)) return false;
    if (!kvstore_copy_all(src, dst, err)) return false;
    mdio_pipe_remove_path(tmp_path);
    return true;
}

/* NaN dead-trace fill (mdio float fill_value). */
static float mdio_pipe_dead_fill(void)
{
    return std::numeric_limits<float>::quiet_NaN();
}

bool mdio_pipe_read_trace_mask(mdio::Dataset& ds,
                               std::vector<char>& live,
                               const std::string& mask_name)
{
    if (!ds.variables.contains_key(mask_name)) {
        /* no mask -> all live */
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

bool mdio_pipe_restamp_stats(mdio::Dataset& ds,
                             const std::string& datavar,
                             const std::string& stats_json,
                             std::string& err)
{
    /* at() so UpdateAttributes mutates shared attrs; do not variables.add */
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
    /* top-level statsV1 (USER_GUIDE layout) */
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

/* Read root zarr.json or .zattrs via kvstore. */
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
    /* Python writer createdOn format */
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

    /* v3: root["attributes"]; v2: .zattrs is attrs doc */
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

/* Set selection flags + contract from per-axis kinds. */
static void mdio_pipe_selection_finalize(MdioSelection& sel)
{
    const int prank = (int) sel.axes.size();
    sel.same_geometry = true;
    sel.sample_changed = false;
    sel.spatial_changed = false;
    sel.exact = true;
    bool resampled = false;

    for (int a = 0; a < prank; a++) {
        const MdioAxisSelKind k = sel.axes[(size_t) a].kind;
        if (k == MDIO_SEL_PRESERVED) continue;
        sel.same_geometry = false;
        if (k == MDIO_SEL_RESAMPLED) sel.exact = false;
        if (k != MDIO_SEL_SUBSET) resampled = true;
        if (a == prank - 1) sel.sample_changed = true;
        else                sel.spatial_changed = true;
    }

    if (sel.same_geometry)
        sel.contract_name = MDIO_PIPE_CONTRACT_SAME_GEOM;
    else if (!sel.spatial_changed)
        sel.contract_name = resampled ? MDIO_PIPE_CONTRACT_RESAMPLE
                                      : MDIO_PIPE_CONTRACT_TRUNCATE;
    else
        sel.contract_name = resampled ? MDIO_PIPE_CONTRACT_DECIMATE
                                      : MDIO_PIPE_CONTRACT_SLICE;
}

/* Classify integer axis selection; derive child o/d. */
static void set_integer_axis_sel(const MdioAxis& pax, long start, long step,
                                 long count, MdioAxisSel& ax)
{
    ax.label = pax.label;
    ax.start = start;
    ax.step  = step;
    ax.count = count;
    if (start == 0 && step == 1 && count == pax.size)
        ax.kind = MDIO_SEL_PRESERVED;
    else if (step == 1)
        ax.kind = MDIO_SEL_SUBSET;
    else
        ax.kind = MDIO_SEL_DECIMATED;
    ax.o = pax.o + (double) start * pax.d;
    ax.d = pax.d * (double) step;
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

    std::vector<std::string> parts;
    std::istringstream is(text);
    for (std::string part; std::getline(is, part, ';'); )
        if (!part.empty()) parts.push_back(part);
    if ((int) parts.size() != prank) {
        err = "mdio_sel has " + std::to_string(parts.size()) +
              " axes, parent has " + std::to_string(prank);
        return false;
    }

    for (int a = 0; a < prank; a++) {
        const MdioAxis& pax = parent_axes[(size_t) a];
        char label[256];
        long start, step, count;
        if (4 != sscanf(parts[(size_t) a].c_str(), "%255[^:]:%ld:%ld:%ld",
                        label, &start, &step, &count)) {
            err = "bad mdio_sel token \"" + parts[(size_t) a] +
                  "\" (want label:start:step:count)";
            return false;
        }
        if (pax.label != label) {
            err = "mdio_sel label \"" + std::string(label) + "\" != parent \"" +
                  pax.label + "\"";
            return false;
        }
        if (step < 1 || count < 1 || start < 0 ||
            start + (count - 1) * step > pax.size - 1) {
            err = "mdio_sel out of range for axis " + pax.label;
            return false;
        }
        set_integer_axis_sel(pax, start, step, count, sel.axes[(size_t) a]);
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
        const MdioAxis& pax = parent_axes[(size_t) a];
        const long start = f[(size_t) a];
        const long step  = jp[(size_t) a] < 1 ? 1 : jp[(size_t) a];
        const long count = m[(size_t) a];
        if (start < 0 || count < 1 ||
            start + (count - 1) * step > pax.size - 1) {
            err = "window out of range for axis " + pax.label;
            return false;
        }
        set_integer_axis_sel(pax, start, step, count, sel.axes[(size_t) a]);
    }
    mdio_pipe_selection_finalize(sel);
    return true;
}

/* Grid-comparison tolerance shared by every detected axis. */
static const double MDIO_SEL_RTOL = 1e-4;
static bool mdio_grid_close(double x, double y)
{
    return std::fabs(x - y) <= MDIO_SEL_RTOL * (1.0 + std::fabs(y)) ||
           std::fabs(x - y) <= 1e-6;
}

/* RSF effective rank: parent rank minus trailing size-1 axes RSF may drop. */
static int mdio_effective_rank(const std::vector<MdioAxis>& parent_axes)
{
    const int prank = (int) parent_axes.size();
    int expect_dim = 1;
    for (int ri = 1; ri <= prank && ri <= SF_MAX_DIM; ri++) {
        int a = prank - ri;
        if (parent_axes[(size_t) a].size > 1) expect_dim = ri;
    }
    return expect_dim;
}

/* Classify one RSF axis against its parent axis (preserved / integer subset /
   decimated / resampled sample axis). Fills ax; false + err when unsupported. */
static bool mdio_detect_one_axis(sf_file in,
                                 const std::vector<MdioAxis>& parent_axes,
                                 const off_t* n, int a,
                                 MdioAxisSel& ax, std::string& err)
{
    const int prank = (int) parent_axes.size();
    const int ri = prank - a;
    const MdioAxis& pax = parent_axes[(size_t) a];

    /* size-1 / RSF-trimmed -> preserved */
    if (pax.size <= 1) {
        set_integer_axis_sel(pax, 0, 1, pax.size < 1 ? 1 : pax.size, ax);
        return true;
    }
    if (pax.d == 0.0) {
        err = "parent axis " + pax.label + " has d=0";
        return false;
    }

    const long rsf_n = (ri <= SF_MAX_DIM) ? (long) n[ri - 1] : 1;
    char key[8];
    float fo = 0.f, fd = 1.f;
    snprintf(key, sizeof(key), "o%d", ri);
    const double rsf_o = sf_histfloat(in, key, &fo) ? (double) fo : pax.o;
    snprintf(key, sizeof(key), "d%d", ri);
    const double rsf_d = sf_histfloat(in, key, &fd) ? (double) fd : pax.d;

    const std::string grid =
        " (child o=" + std::to_string(rsf_o) +
        " d=" + std::to_string(rsf_d) + " n=" + std::to_string(rsf_n) +
        "; parent o=" + std::to_string(pax.o) +
        " d=" + std::to_string(pax.d) +
        " n=" + std::to_string(pax.size) + ")";

    const double start_f = (rsf_o - pax.o) / pax.d;
    const double step_f  = rsf_d / pax.d;
    long start = (long) std::llround(start_f);
    long step  = (long) std::llround(step_f);
    if (step < 1) step = 1;

    if (mdio_grid_close(start_f, (double) start) &&
        mdio_grid_close(step_f, (double) step) &&
        mdio_grid_close(pax.o + start * pax.d, rsf_o) &&
        mdio_grid_close(pax.d * step, rsf_d)) {
        if (start < 0 || rsf_n < 1 ||
            start + (rsf_n - 1) * step > pax.size - 1) {
            err = "axis " + pax.label + " selects outside the parent" + grid;
            return false;
        }
        set_integer_axis_sel(pax, start, step, rsf_n, ax);
        return true;
    }

    /* off-grid: sample axis only may resample in-range (sfremap1) */
    if (a != prank - 1) {
        err = "axis " + pax.label + " is not an exact integer selection "
              "of the parent" + grid;
        return false;
    }
    if (!(rsf_d > 0.0) || rsf_n < 1) {
        err = "resampled sample axis \"" + pax.label +
              "\" needs d > 0 and n >= 1" + grid;
        return false;
    }
    const double parent_last = pax.o + (double) (pax.size - 1) * pax.d;
    const double child_last  = rsf_o + (double) (rsf_n - 1) * rsf_d;
    const double atol =
        1e-6 + MDIO_SEL_RTOL * (1.0 + std::fabs(parent_last) + std::fabs(pax.o));
    if (rsf_o < pax.o - atol || child_last > parent_last + atol) {
        err = "resampled sample axis \"" + pax.label +
              "\" extrapolates past the parent range" + grid;
        return false;
    }
    ax.label = pax.label;
    ax.start = 0;
    ax.step  = 0;
    ax.count = rsf_n;
    ax.kind  = MDIO_SEL_RESAMPLED;
    ax.o = rsf_o;
    ax.d = rsf_d;
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

    const int expect_dim = mdio_effective_rank(parent_axes);
    if (rsf_dim != expect_dim) {
        err = "RSF effective rank " + std::to_string(rsf_dim) +
              " != parent effective rank " + std::to_string(expect_dim) +
              " (unsupported geometry change)";
        return false;
    }

    sel = MdioSelection();
    sel.axes.resize((size_t) prank);
    for (int a = 0; a < prank; a++)
        if (!mdio_detect_one_axis(in, parent_axes, n, a,
                                  sel.axes[(size_t) a], err))
            return false;

    mdio_pipe_selection_finalize(sel);
    return true;
}

std::vector<long> mdio_pipe_child_sizes(const MdioSelection& sel)
{
    std::vector<long> out(sel.axes.size());
    for (size_t i = 0; i < sel.axes.size(); i++)
        out[i] = sel.axes[i].count;
    return out;
}

/* Chunk shape via mdio-cpp; Variable::get_chunk_shape() resolves the Zarr v2
   ("chunks") vs v3 ("chunk_grid") layout for us. */
static bool mdio_pipe_chunk_shape(const std::string& store_uri,
                                  const std::string& datavar,
                                  std::vector<long>& chunk,
                                  std::string& err)
{
    chunk.clear();
    auto dsr = mdio::Dataset::Open(store_uri, mdio::constants::kOpen);
    if (!dsr.status().ok()) {
        err = "cannot open \"" + store_uri + "\": " + dsr.status().ToString();
        return false;
    }
    mdio::Dataset ds = dsr.value();

    auto vr = ds.variables.at(datavar);
    if (!vr.status().ok()) {
        err = "no variable \"" + datavar + "\" under " + store_uri;
        return false;
    }

    auto csr = vr.value().get_chunk_shape();
    if (!csr.status().ok()) {
        err = datavar + ": " + csr.status().ToString();
        return false;
    }
    for (auto c : csr.value()) chunk.push_back((long) c);
    return true;
}

/* Pick pivot axis + block size; false if no whole-chunk block fits budget. */
static bool plan_blocks(const std::vector<long>& sizes,
                        const std::vector<long>& chunk,
                        long budget_floats,
                        int& pivot, long& pivot_block)
{
    const int rank = (int) sizes.size();
    if (rank < 2 || (int) chunk.size() != rank) return false;
    for (int a = 0; a < rank; a++)
        if (sizes[(size_t) a] < 1 || chunk[(size_t) a] < 1) return false;
    if (budget_floats < 1) budget_floats = 1;

    /* clamp chunk extent to axis size before fit test */
    std::vector<long> eff((size_t) rank);
    for (int a = 0; a < rank; a++)
        eff[(size_t) a] = std::min(chunk[(size_t) a], sizes[(size_t) a]);

    /* tail[a] = floats in all faster axes */
    std::vector<long> tail((size_t) rank, 1);
    for (int a = rank - 2; a >= 0; a--)
        tail[(size_t) a] = tail[(size_t) (a + 1)] * sizes[(size_t) (a + 1)];

    /* never pivot on sample axis; stop at first eff chunk > 1 */
    int last_pivot = rank - 2;
    for (int a = 0; a < rank - 2; a++)
        if (eff[(size_t) a] > 1) { last_pivot = a; break; }

    /* slowest axis that fits -> widest whole-chunk Write span */
    for (int k = 0; k <= last_pivot; k++) {
        const long one = eff[(size_t) k] * tail[(size_t) k];
        if (one > budget_floats) continue;
        pivot = k;
        pivot_block = std::min((budget_floats / one) * eff[(size_t) k],
                               sizes[(size_t) k]);
        return true;
    }

    /* one chunk over budget (<=8x): cheaper than R-M-W per chunk */
    const long one = eff[(size_t) last_pivot] * tail[(size_t) last_pivot];
    if (one > budget_floats * 8) return false;
    pivot = last_pivot;
    pivot_block = eff[(size_t) last_pivot];
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
                                  int panel_cap,
                                  int blockmb,
                                  MdioBlockLoop& loop,
                                  std::string& note)
{
    const int rank = (int) sizes.size();
    if (blockmb < 1) blockmb = MDIO_PIPE_BLOCK_MB;
    const long budget = (long) blockmb * 1024 * 1024 / 4;

    std::vector<long> chunk;
    std::string cerr;
    if (!mdio_pipe_chunk_shape(store_uri, datavar, chunk, cerr) ||
        (int) chunk.size() != rank) {
        chunk.clear();
        if (!cerr.empty()) note = "chunk shape unavailable (" + cerr + ")";
    }

    loop.aligned = rank >= 2 && !chunk.empty() &&
        plan_blocks(sizes, chunk, budget, loop.pivot, loop.pivot_block);

    if (!loop.aligned && rank >= 2) {
        /* straddle chunks on fastest spatial axis (slow R-M-W path) */
        const long ns = std::max(1L, sizes[(size_t) (rank - 1)]);
        loop.pivot       = rank - 2;
        loop.pivot_block = std::max(1L, std::min(budget / ns,
                                                 sizes[(size_t) (rank - 2)]));
        /* panel= caps traces per fallback block along the fastest spatial axis */
        if (panel_cap >= 1 && loop.pivot_block > (long) panel_cap)
            loop.pivot_block = (long) panel_cap;
        note += (note.empty() ? "" : "; ");
        note += "no whole-chunk block fits " + std::to_string(blockmb) +
                " MiB; blocking " + std::to_string(loop.pivot_block) +
                " along the fastest spatial axis (slow: expect chunk "
                "read-modify-write)";
    } else if (!loop.aligned) {
        loop.pivot       = 0;
        loop.pivot_block = sizes[0];    /* single axis */
    }

    /* slower axes stepped; faster axes full extent */
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

/* data_type as string or object.name. */
static std::string json_dtype_name(const nlohmann::json& meta)
{
    if (!meta.contains("data_type")) return "";
    const nlohmann::json& dt = meta["data_type"];
    if (dt.is_string()) return dt.get<std::string>();
    if (dt.is_object()) return dt.value("name", "");
    return "";
}

/* forward decls for rematerialize */
static void cap_chunk_shape(nlohmann::json& meta);
static void normalize_structured_dtype(nlohmann::json& meta,
                                       const nlohmann::json* parent_meta);

/* Changed selection axis for label, or NULL if preserved. */
static const MdioAxisSel* changed_axis_for(const MdioSelection& sel,
                                           const std::string& label)
{
    for (size_t a = 0; a < sel.axes.size(); a++)
        if (sel.axes[a].label == label &&
            sel.axes[a].kind != MDIO_SEL_PRESERVED)
            return &sel.axes[a];
    return NULL;
}

/* Trim variable zarr.json shape/chunks to selection (kvstore path). */
static bool reshape_var_metadata(tensorstore::KvStore& child_kvs,
                                 tensorstore::KvStore* parent_kvs,
                                 const std::string& name,
                                 const MdioSelection& sel,
                                 std::string& err)
{
    nlohmann::json meta;
    if (!kv_read_var_meta(child_kvs, name, meta, err)) {
        err = "missing or unreadable zarr.json for \"" + name + "\": " + err;
        return false;
    }

    /* restore parent dtype/codecs before shape edit */
    nlohmann::json parent_meta;
    const nlohmann::json* pmeta = NULL;
    if (parent_kvs) {
        std::string perr;
        if (kv_read_var_meta(*parent_kvs, name, parent_meta, perr))
            pmeta = &parent_meta;
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
        const MdioAxisSel* ax =
            changed_axis_for(sel, dims[i].get<std::string>());
        if (!ax) continue;
        shape[i] = ax->count;
        changed = true;
    }
    if (!changed) return true;

    /* cap chunk_shape to new shape */
    if (!meta.contains("chunk_grid") || !meta["chunk_grid"].is_object())
        meta["chunk_grid"] = {{"name", "regular"},
                              {"configuration", {{"chunk_shape", shape}}}};
    else {
        meta["chunk_grid"]["name"] = "regular";
        meta["chunk_grid"]["configuration"]["chunk_shape"] = shape;
    }
    /* multi-dim: parent chunk on preserved dims, cap on changed */
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

    if (!kv_write_var_meta(child_kvs, name, meta, err)) return false;

    /* remove stale chunks before resized Write */
    return kv_delete_prefix(child_kvs, join_key(name, "c") + "/", err);
}

/* Variable touches any changed selection axis? */
static bool var_needs_reshape(const nlohmann::json& meta,
                              const MdioSelection& sel)
{
    if (!meta.contains("dimension_names") ||
        !meta["dimension_names"].is_array())
        return false;
    for (auto& d : meta["dimension_names"])
        if (d.is_string() && changed_axis_for(sel, d.get<std::string>()))
            return true;
    return false;
}

/* Write 1-D var; false (err untouched) if wrong MdioT. */
template <typename MdioT>
static bool write_1d_variable(mdio::Dataset& ds, const std::string& label,
                              const std::vector<double>& vals,
                              std::string& err)
{
    auto vr = ds.variables.get<MdioT>(label);
    if (!vr.status().ok()) return false;
    auto var = vr.value();
    auto vdr = mdio::from_variable<MdioT>(var);
    if (!vdr.status().ok()) {
        err = "from_variable \"" + label + "\": " + vdr.status().ToString();
        return false;
    }
    auto vd = vdr.value();
    if ((long) vd.num_samples() != (long) vals.size()) {
        err = "coordinate \"" + label + "\" holds " +
              std::to_string((long) vd.num_samples()) + " values, expected " +
              std::to_string(vals.size());
        return false;
    }
    auto off = vd.get_flattened_offset();
    MdioT* p = static_cast<MdioT*>(vd.get_data_accessor().data());
    for (size_t i = 0; i < vals.size(); i++)
        p[off + i] = static_cast<MdioT>(
            std::is_integral<MdioT>::value ? (double) std::llround(vals[i])
                                           : vals[i]);
    if (!var.Write(vd).status().ok()) {
        err = "Write coordinate \"" + label + "\" failed";
        return false;
    }
    return true;
}

/* Resampled sample coord o+i*d; int dtype if all integral else float64. */
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

    std::vector<double> vals((size_t) sax.count);
    bool all_integral = true;
    for (long i = 0; i < sax.count; i++) {
        vals[(size_t) i] = sax.o + (double) i * sax.d;
        if (std::fabs(vals[(size_t) i] - std::round(vals[(size_t) i])) > 1e-9)
            all_integral = false;
    }

    tensorstore::KvStore child_kvs, parent_kvs;
    if (!open_root_kvstore(child_path, child_kvs, err)) return false;

    nlohmann::json meta, pmeta;
    if (!kv_read_var_meta(child_kvs, label, meta, err)) {
        err = "missing or unreadable sample coord metadata for \"" + label +
              "\": " + err;
        return false;
    }
    bool have_parent = false;
    std::string perr;
    if (open_root_kvstore(parent_uri, parent_kvs, perr) &&
        kv_read_var_meta(parent_kvs, label, pmeta, perr))
        have_parent = true;
    const std::string parent_dtype =
        have_parent ? json_dtype_name(pmeta) : std::string();
    const bool keep_int = all_integral &&
        (parent_dtype.rfind("int", 0) == 0 || parent_dtype.rfind("uint", 0) == 0);
    const std::string out_dtype = keep_int ? parent_dtype : "float64";

    /* set on-disk dtype before Open for Write */
    if (json_dtype_name(meta) != out_dtype) {
        meta["data_type"] = out_dtype;
        if (!kv_write_var_meta(child_kvs, label, meta, err)) return false;
        if (!kv_delete_prefix(child_kvs, join_key(label, "c") + "/", err))
            return false;
    }

    auto cf = mdio::Dataset::Open(child_path, mdio::constants::kOpen);
    if (!cf.status().ok()) {
        err = "reopen child for synthesized coord: " + cf.status().ToString();
        return false;
    }
    mdio::Dataset child = cf.value();

    if (write_1d_variable<mdio::dtypes::float64_t>(child, label, vals, err))
        return true;
    if (err.empty() &&
        write_1d_variable<mdio::dtypes::int32_t>(child, label, vals, err))
        return true;
    if (err.empty())
        err = "cannot open sample coord \"" + label + "\" as float64 or "
              "int32 after setting its dtype to " + out_dtype;
    return false;
}

static bool rematerialize_child_variables(mdio::Dataset& parent,
                                          const std::string& parent_uri,
                                          const std::string& child_path,
                                          const MdioSelection& sel,
                                          const std::string& datavar,
                                          std::string& err)
{
    tensorstore::KvStore child_kvs, parent_kvs;
    if (!open_root_kvstore(child_path, child_kvs, err)) return false;
    tensorstore::KvStore* parent_kvs_ptr = NULL;
    std::string perr;
    if (open_root_kvstore(parent_uri, parent_kvs, perr))
        parent_kvs_ptr = &parent_kvs;

    /* Enumerate vars via Dataset (works for local and remote). */
    auto cf0 = mdio::Dataset::Open(child_path, mdio::constants::kOpen);
    if (!cf0.status().ok()) {
        err = "open child for reshape: " + cf0.status().ToString();
        return false;
    }
    mdio::Dataset child0 = cf0.value();

    /* pass 1: reshape metadata + wipe chunks (incl. amplitude / headers) */
    std::vector<std::string> reshape_names = child0.variables.get_keys();
    for (const std::string& hk : child0.header_variables.get_keys())
        reshape_names.push_back(hk);

    std::vector<std::string> to_copy;
    for (const std::string& name : reshape_names) {
        nlohmann::json meta;
        std::string merr;
        if (!kv_read_var_meta(child_kvs, name, meta, merr)) continue;
        if (!var_needs_reshape(meta, sel)) continue;
        if (!reshape_var_metadata(child_kvs, parent_kvs_ptr, name, sel, err))
            return false;
        /* payload rematerialize is for regular variables only */
        if (name != datavar && child0.variables.contains_key(name))
            to_copy.push_back(name);
    }

    /* reopen child; Copy sliced parent data */
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
            /* resampled sample axis: synthesize sample coord only */
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

    /* metadata on disk; Copy is payloads only — skip CommitMetadata */
    (void) child;
    return true;
}

/* Cap chunk_shape[i] to shape[i] after Resize. */
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

/* struct -> structured spelling for Python MDIO. */
static void normalize_structured_dtype(nlohmann::json& meta,
                                       const nlohmann::json* parent_meta)
{
    if (!meta.contains("data_type") || !meta["data_type"].is_object()) return;
    auto& dt = meta["data_type"];
    std::string name = dt.value("name", "");
    if (name != "struct" && name != "structured") return;

    /* prefer parent dtype/fill; do not copy codecs (child may differ) */
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

    /* object fill_value -> base64 zero struct for Python MDIO */
    if (meta.contains("fill_value") && meta["fill_value"].is_object()) {
        size_t nbytes = 0;
        for (auto& f : fields) {
            if (!f.is_array() || f.size() < 2) continue;
            size_t w = mdio_dtype_width(f[1].get<std::string>());
            nbytes += w ? w : 4;
        }
        std::string enc;
        for (size_t i = 0; i < nbytes; i += 3)
            enc += (nbytes - i >= 3) ? "AAAA" :
                   (nbytes - i == 2 ? "AAA=" : "AA==");
        meta["fill_value"] = enc;
    }
}

/* Repair per-var zarr.json after TensorStore writes (Python MDIO readability). */
bool mdio_pipe_repair_zarr_metadata(const std::string& child_path,
                                    const std::string& parent_uri,
                                    std::string& err,
                                    bool cap_chunks)
{
    tensorstore::KvStore child_kvs;
    if (!open_root_kvstore(child_path, child_kvs, err)) return false;

    tensorstore::KvStore parent_kvs;
    tensorstore::KvStore* parent_kvs_ptr = NULL;
    if (!parent_uri.empty()) {
        std::string perr;
        if (open_root_kvstore(parent_uri, parent_kvs, perr))
            parent_kvs_ptr = &parent_kvs;
    }

    auto cf = mdio::Dataset::Open(child_path, mdio::constants::kOpen);
    if (!cf.status().ok()) {
        err = "repair: cannot open \"" + child_path +
              "\": " + cf.status().ToString();
        return false;
    }
    mdio::Dataset child = cf.value();

    std::vector<std::string> names = child.variables.get_keys();
    for (const std::string& hk : child.header_variables.get_keys())
        names.push_back(hk);

    for (const std::string& name : names) {
        nlohmann::json meta;
        if (!kv_read_var_meta(child_kvs, name, meta, err)) {
            err = "repair: cannot read metadata for \"" + name + "\": " + err;
            return false;
        }

        nlohmann::json parent_meta;
        const nlohmann::json* pmeta = NULL;
        if (parent_kvs_ptr) {
            std::string perr;
            if (kv_read_var_meta(*parent_kvs_ptr, name, parent_meta, perr))
                pmeta = &parent_meta;
        }

        normalize_structured_dtype(meta, pmeta);
        /* cap chunks only after all Writes complete */
        if (cap_chunks) cap_chunk_shape(meta);

        /* default chunk_key_encoding separator */
        if (meta.contains("chunk_key_encoding") &&
            meta["chunk_key_encoding"].is_object() &&
            !meta["chunk_key_encoding"].contains("configuration")) {
            meta["chunk_key_encoding"]["configuration"] = {
                {"separator", "/"}};
        }

        if (!kv_write_var_meta(child_kvs, name, meta, err)) {
            err = "repair: cannot write metadata for \"" + name + "\": " + err;
            return false;
        }
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

    /* reshape + isel copy; dtype repaired before publish */
    if (!rematerialize_child_variables(parent, parent_uri, tmp_path, sel,
                                       datavar, err)) {
        mdio_pipe_remove_path(tmp_path);
        return false;
    }
    return true;
}

/* SEG-Y sample geometry for header restamp. */
struct SampleGeometryStamp {
    long ns;        /* samples per trace after windowing           */
    int  dt_us;     /* sample interval, microseconds               */
    int  shift_ms;  /* added to delay_recording_time, milliseconds */
};

/* Struct field byte offsets from zarr.json data_type (fail-closed). */
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
        /* [name,type] or object field spellings */
        std::string name, typ;
        if (f.is_array() && f.size() >= 2) {
            name = f[0].get<std::string>();
            typ  = f[1].get<std::string>();
        } else if (f.is_object()) {
            name = f.value("name", "");
            typ  = f.value("data_type", "");
        } else continue;

        const size_t w = mdio_dtype_width(typ);
        if (0 == w) {
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

/* Patch SEG-Y fields via headers byte view (preserves dtype/codec/chunks). */
static bool stamp_headers_sample_geometry(const std::string& child_path,
                                          const SampleGeometryStamp& geom,
                                          const std::vector<char>& live,
                                          std::string& err)
{
    tensorstore::KvStore child_kvs;
    if (!open_root_kvstore(child_path, child_kvs, err)) return false;

    nlohmann::json meta;
    std::string merr;
    if (!kv_read_var_meta(child_kvs, "headers", meta, merr))
        return true; /* no headers */

    std::unordered_map<std::string, size_t> foff;
    size_t itemsize = 0;
    if (!structured_field_offsets(meta, foff, itemsize, err)) return false;

    /* required int16 SEG-Y geometry fields */
    static const char* kGeomFields[] = {"samples_per_trace", "sample_interval",
                                        "delay_recording_time"};
    for (const char* name : kGeomFields) {
        auto it = foff.find(name);
        if (it == foff.end()) {
            err = std::string("headers missing required field \"") + name +
                  "\" for sample-geometry stamp at " + child_path;
            return false;
        }
        size_t next = itemsize;
        for (const auto& kv : foff)
            if (kv.second > it->second && kv.second < next) next = kv.second;
        if (next - it->second != 2) {
            err = std::string("sample-geometry stamp expects int16 \"") + name +
                  "\"; found " + std::to_string(next - it->second) +
                  " bytes at " + child_path;
            return false;
        }
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

    /* open_as_void byte/uint8 + trailing byte axis */
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
              " < 2 at " + child_path;
        return false;
    }
    auto shape = store.domain().shape();
    const long byte_dim = (long) shape[rank - 1];
    if (byte_dim != (long) itemsize) {
        err = "headers trailing byte dim " + std::to_string(byte_dim) +
              " != structured itemsize " + std::to_string(itemsize) +
              " at " + child_path;
        return false;
    }
    long ntr = 1;
    for (int i = 0; i < rank - 1; i++) ntr *= (long) shape[i];
    if (ntr < 1) {
        err = "headers byte-view has no traces at " + child_path;
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

/* Patch binaryHeader in segy_file_header zarr.json (avoid CommitMetadata). */
static bool stamp_binary_header_ns_dt(const std::string& child_path,
                                      const SampleGeometryStamp& geom,
                                      std::string& err)
{
    tensorstore::KvStore child_kvs;
    if (!open_root_kvstore(child_path, child_kvs, err)) return false;

    nlohmann::json meta;
    std::string merr;
    if (!kv_read_var_meta(child_kvs, "segy_file_header", meta, merr))
        return true; /* optional */

    if (!meta.contains("attributes") || !meta["attributes"].is_object())
        return true;
    nlohmann::json& attrs = meta["attributes"];
    if (!attrs.contains("binaryHeader") || !attrs["binaryHeader"].is_object())
        return true;
    attrs["binaryHeader"]["samples_per_trace"] = (int) geom.ns;
    attrs["binaryHeader"]["sample_interval"]   = geom.dt_us;

    return kv_write_var_meta(child_kvs, "segy_file_header", meta, err);
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

    /* SEG-Y geometry from child o/d (works for resample too) */
    SampleGeometryStamp geom;
    geom.ns       = sax.count;
    geom.dt_us    = (int) std::llround(sax.d * 1000.0);
    geom.shift_ms = (int) std::llround(sax.o - pax.o);

    if (!stamp_binary_header_ns_dt(child_path, geom, err))
        return false;

    /* need trace_mask — dead traces keep original delay */
    std::vector<char> live;
    if (!mdio_pipe_read_trace_mask(child, live)) {
        err = "failed to read child trace_mask for the sample-geometry stamp "
              "at " + child_path + " (datavar=" + datavar + ")";
        return false;
    }

    if (!stamp_headers_sample_geometry(child_path, geom, live, err))
        return false;

    /* repair headers dtype spelling; no chunk cap (shapes unchanged) */
    return mdio_pipe_repair_zarr_metadata(child_path, parent_uri, err,
                                          /*cap_chunks=*/false);
}

