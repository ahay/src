/* Convert an RSF dataset to MDIO.

Writes an RSF amplitude file to the MDIO (Zarr) format through the mdio-cpp
library.  The default path is the faithful parent-clone finalizer: require
mdio_* pipe context (or like=/context= recovery), reopen the parent, verify
the fingerprint, detect geometry (same-geometry, truncate/slice/decimate, or
interpolating sample-axis resample), clone or reshape the parent store, stream-
overwrite only the data variable from the RSF stream (with dead-trace fill),
recompute statsV1, stamp sample-geometry headers when the sample axis changed,
append provenance, and publish.

Sample-axis windows (f1/n1/j1 from sfmdioread) and interpolating resample
composed as `sfmdioread | sfremap1 n1=... d1=... o1=... | sfmdiowrite` are
supported.  sfremap1 is ENO grid remapping, not anti-aliased resampling — for
signal-processing-quality downsampling put a low-pass stage before sfremap1.

Pass reduced=y to opt into the legacy reduced-schema path (lossy; will not pass
mdio_conv Axis A/B checks).  With reduced=y and no upstream pipe context (no
sfmdioread, like=, or context=), it originates a fresh MDIO from the RSF stream
alone — the only way to seed a parent store from Madagascar, e.g. for examples
and smoke tests.

The output may be a local path, or a gs:// or s3:// URL (when mdio-cpp was
built with the corresponding cloud drivers).  Cloud publish copies a temp
prefix then deletes it (S3/GCS have no atomic directory rename).

Amplitude is streamed in bounded blocks that each span whole Zarr chunks, so
every chunk is compressed and written exactly once.  Blocking only the fastest
spatial axis instead (the fallback) makes TensorStore read-modify-write and
re-compress a whole chunk once per index of every coarser-chunked axis.  See
blockmb= and panel= to tune the block sizing.

For processing=, prefer a path to a .json file over inline JSON: a shell eats
the quotes, and unquoted {...} triggers brace expansion.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "mdio2segy.hh"
#include "mdio_pipe.hh"

/* Sanitize an RSF axis label into a valid, unique MDIO dimension name. */
static std::string clean_name(const char* label, int rsfaxis,
                              std::vector<std::string>& used)
{
    std::string s;
    if (label) {
        for (const char* p = label; *p; p++)
            s += (isalnum((unsigned char) *p) || *p == '_') ? *p : '_';
    }
    if (s.empty()) { char b[16]; snprintf(b, sizeof(b), "dim%d", rsfaxis); s = b; }

    std::string base = s;
    int suffix = 1;
    bool clash = true;
    while (clash) {
        clash = false;
        for (size_t i = 0; i < used.size(); i++)
            if (used[i] == s) { clash = true; break; }
        if (clash) { char b[16]; snprintf(b, sizeof(b), "_%d", suffix++); s = base + b; }
    }
    used.push_back(s);
    return s;
}

static bool write_float_buf(mdio::Dataset& ds, const std::string& name,
                            const float* buf, long total)
{
    auto vr = ds.variables.get<mdio::dtypes::float32_t>(name);
    if (!vr.status().ok()) return false;
    auto var = vr.value();
    auto vdr = mdio::from_variable<mdio::dtypes::float32_t>(var);
    if (!vdr.status().ok()) return false;
    auto vd = vdr.value();
    long n = (long) vd.num_samples();
    auto off = vd.get_flattened_offset();
    float* p = static_cast<float*>(vd.get_data_accessor().data());
    for (long i = 0; i < n && i < total; i++) p[off + i] = buf[i];
    return var.Write(vd).status().ok();
}

/* Flatten the RSF stream geometry: total = product of all axes, ns = fast
   (sample) axis length, ntr = total / ns. */
static void rsf_stream_dims(sf_file in, long& ns, long& ntr, long& total,
                            int& dim)
{
    off_t n[SF_MAX_DIM];
    dim = sf_largefiledims(in, n);
    ns = (long) n[0];
    total = 1;
    for (int i = 0; i < dim; i++) total *= (long) n[i];
    ntr = total / ns;
}

/* Removes a temp MDIO store on scope exit unless commit() was called.
   Keeps the fail-closed finalizer from ever leaving a partial child on disk
   without repeating the cleanup at every error branch. */
struct TmpStoreGuard {
    std::string path;
    bool        armed;
    explicit TmpStoreGuard(std::string p) : path(std::move(p)), armed(true) {}
    ~TmpStoreGuard() { if (armed) mdio_pipe_remove_path(path); }
    void commit() { armed = false; }
    TmpStoreGuard(const TmpStoreGuard&) = delete;
    TmpStoreGuard& operator=(const TmpStoreGuard&) = delete;
};

static bool write_int_var(mdio::Dataset& ds, const std::string& name,
                          const std::vector<int>& col)
{
    auto vr = ds.variables.get<mdio::dtypes::int32_t>(name);
    if (!vr.status().ok()) return false;
    auto var = vr.value();
    auto vdr = mdio::from_variable<mdio::dtypes::int32_t>(var);
    if (!vdr.status().ok()) return false;
    auto vd = vdr.value();
    long n = (long) vd.num_samples();
    auto off = vd.get_flattened_offset();
    mdio::dtypes::int32_t* p =
        static_cast<mdio::dtypes::int32_t*>(vd.get_data_accessor().data());
    for (long i = 0; i < n && i < (long) col.size(); i++)
        p[off + i] = (mdio::dtypes::int32_t) col[(size_t) i];
    return var.Write(vd).status().ok();
}

/* The three dataset-valued parameters, read exactly once each so that sfdoc
   lists each with its own description.  context= falls back to headers= then
   like= for pipe-context recovery; the reduced path uses headers= then like=
   as its header-copy source. */
struct DatasetParams {
    char* context;
    char* headers;
    char* like;

    const char* recovery_source() const
    {
        if (NULL != context) return context;
        if (NULL != headers) return headers;
        return like;
    }
    const char* header_source() const
    {
        return (NULL != headers) ? headers : like;
    }
};

static DatasetParams get_dataset_params(void)
{
    DatasetParams p;

    p.context = sf_getstring("context");
    /* MDIO parent to recover mdio_* pipe context from, when an intermediate
       filter dropped the keys; falls back to headers= then like= */

    p.headers = sf_getstring("headers");
    /* MDIO dataset to copy text/binary/trace headers from (reduced=y) */

    p.like = sf_getstring("like");
    /* alias for headers=, and last fallback for context= */

    return p;
}

/* Recover pipe context from a parent MDIO when an intermediate filter stripped
   the mdio_* keys off the stream. */
static void recover_pipe_context(const char* source, const char* dataname,
                                 MdioPipeContext& ctx, bool verb)
{
    auto sf = mdio::Dataset::Open(std::string(source), mdio::constants::kOpen);
    if (!sf.status().ok())
        sf_error("Cannot open context recovery source \"%s\"", source);
    mdio::Dataset rds = sf.value();

    std::string rvar = mdio_data_variable(rds, dataname);
    if (rvar.empty())
        sf_error("No data variable in context source \"%s\"", source);
    if (!mdio_pipe_context_from_parent(rds, std::string(source), rvar, ctx))
        sf_error("Context recovery from \"%s\" failed "
                 "(float32 data variable required)", source);
    if (verb) sf_warning("Recovered mdio_* context from \"%s\"", source);
}

/* Fill in defaults and reject anything this build cannot honour, so the
   finalizer never runs on a context it only half understands. */
static void validate_pipe_context(MdioPipeContext& ctx)
{
    if (ctx.version.empty())
        ctx.version = MDIO_PIPE_CONTEXT_VERSION;
    if (ctx.version != MDIO_PIPE_CONTEXT_VERSION)
        sf_error("Unsupported mdio_context_version=\"%s\" (need \"%s\")",
                 ctx.version.c_str(), MDIO_PIPE_CONTEXT_VERSION);

    if (ctx.dtype.empty()) ctx.dtype = MDIO_PIPE_DTYPE_FLOAT32;
    if (ctx.dtype != MDIO_PIPE_DTYPE_FLOAT32)
        sf_error("Unsupported mdio_data_dtype=\"%s\"; the lossless pipe "
                 "carries float32 only (the RSF stream between sfmdioread and "
                 "sfmdiowrite is SF_FLOAT)", ctx.dtype.c_str());

    if (ctx.contract.empty())
        ctx.contract = MDIO_PIPE_CONTRACT_SAME_GEOM;
    if (ctx.contract != MDIO_PIPE_CONTRACT_SAME_GEOM &&
        ctx.contract != MDIO_PIPE_CONTRACT_TRUNCATE &&
        ctx.contract != MDIO_PIPE_CONTRACT_SLICE &&
        ctx.contract != MDIO_PIPE_CONTRACT_DECIMATE &&
        ctx.contract != MDIO_PIPE_CONTRACT_RESAMPLE)
        sf_error("Unsupported mdio_contract=\"%s\"", ctx.contract.c_str());

    if (ctx.data_variable.empty())
        sf_error("mdio_data_variable missing from pipe context");
    if (ctx.source.empty() || ctx.fingerprint.empty())
        sf_error("Incomplete mdio_* pipe context "
                 "(need mdio_source and mdio_source_fingerprint)");
}

/* Reopen the parent named by the context and prove it is the dataset the
   stream was read from, on both dtype and metadata fingerprint. */
static mdio::Dataset open_verified_parent(const MdioPipeContext& ctx, bool verb)
{
    auto pf = mdio::Dataset::Open(ctx.source, mdio::constants::kOpen);
    if (!pf.status().ok())
        sf_error("Cannot open parent MDIO \"%s\" for fingerprint check",
                 ctx.source.c_str());
    mdio::Dataset parent = pf.value();

    /* Compare against the parent's actual dtype rather than re-probing for
       float32: that is what makes mdio_data_dtype a claim the parent has to
       answer for, instead of a constant both sides echo. */
    std::string pdtype = mdio_variable_dtype(parent, ctx.data_variable);
    if (pdtype != ctx.dtype)
        sf_error("Parent data variable \"%s\" has dtype \"%s\", "
                 "but the stream declares mdio_data_dtype=\"%s\"",
                 ctx.data_variable.c_str(),
                 pdtype.empty() ? "unknown" : pdtype.c_str(),
                 ctx.dtype.c_str());

    std::string expect =
        mdio_pipe_fingerprint(parent, ctx.data_variable, ctx.dtype);
    if (expect != ctx.fingerprint)
        sf_error("Parent fingerprint mismatch for \"%s\":\n"
                 "  stream  %s\n"
                 "  parent  %s",
                 ctx.source.c_str(), ctx.fingerprint.c_str(), expect.c_str());
    if (verb) sf_warning("Parent fingerprint OK (%s)", ctx.fingerprint.c_str());

    return parent;
}

/* Prefer the stream's actual n#/o#/d# geometry.  The mdio_sel stamp is only
   supporting context: use it when it agrees with the stream, otherwise the
   stream wins (a middle filter such as sfremap1 may legitimately have changed
   the grid).  Fails closed when the stream is not a valid selection. */
static bool selections_agree(const MdioSelection& a, const MdioSelection& b)
{
    if (a.axes.size() != b.axes.size()) return false;
    const double rtol = 1e-4;
    auto close = [rtol](double x, double y) {
        return std::fabs(x - y) <= rtol * (1.0 + std::fabs(y)) ||
               std::fabs(x - y) <= 1e-6;
    };
    for (size_t i = 0; i < a.axes.size(); i++) {
        const MdioAxisSel& aa = a.axes[i];
        const MdioAxisSel& bb = b.axes[i];
        if (aa.label != bb.label || aa.count != bb.count) return false;
        if (aa.kind == MDIO_SEL_RESAMPLED || bb.kind == MDIO_SEL_RESAMPLED) {
            if (aa.kind != bb.kind) return false;
            if (!close(aa.o, bb.o) || !close(aa.d, bb.d)) return false;
        } else {
            if (aa.start != bb.start || aa.step != bb.step) return false;
        }
    }
    return true;
}

static MdioSelection resolve_selection(sf_file in, const MdioPipeContext& ctx,
                                       const std::vector<MdioAxis>& parent_axes,
                                       bool verb)
{
    MdioSelection stream_sel;
    std::string gerr;
    if (!mdio_pipe_detect_selection(in, parent_axes, stream_sel, gerr))
        sf_error("%s", gerr.c_str());

    MdioSelection sel = stream_sel;
    if (!ctx.selection.empty()) {
        MdioSelection stamp_sel;
        std::string perr;
        if (mdio_pipe_selection_parse(ctx.selection, parent_axes, stamp_sel,
                                      perr)) {
            if (selections_agree(stamp_sel, stream_sel)) {
                sel = stamp_sel;
                if (verb)
                    sf_warning("mdio_sel stamp agrees with stream geometry");
            } else if (verb) {
                sf_warning("mdio_sel stamp disagrees with stream "
                           "(stamp=%s stream=%s); using stream geometry",
                           mdio_pipe_selection_encode(stamp_sel).c_str(),
                           mdio_pipe_selection_encode(stream_sel).c_str());
            }
        } else if (verb) {
            sf_warning("mdio_sel stamp unusable (%s); using stream geometry",
                       perr.c_str());
        }
    }

    if (verb)
        sf_warning("Selection OK contract=%s sel=%s",
                   sel.contract_name.c_str(),
                   mdio_pipe_selection_encode(sel).c_str());
    return sel;
}

/* Clone (same-geometry) or reshape (every other contract) the parent into a
   temp store, ready for the amplitude overwrite. */
static void build_child_store(const MdioPipeContext& ctx,
                              const MdioSelection& sel,
                              const std::string& datavar,
                              const std::string& tmp_path,
                              mdio::Dataset& parent)
{
    std::string err;
    if (sel.same_geometry) {
        /* Skip copying the parent's amplitude chunks — the stream overwrites
           every sample anyway. */
        if (!mdio_clone_store_skip_var(ctx.source, tmp_path, datavar, err))
            sf_error("Parent clone failed: %s", err.c_str());
    } else {
        if (!mdio_pipe_build_child_store(parent, ctx.source, sel, datavar,
                                         tmp_path, err))
            sf_error("Reshaped child build failed: %s", err.c_str());
    }
}

/* The child's (possibly sliced) trace mask, all-live when it has none. */
static std::vector<char> child_trace_mask(mdio::Dataset& child, long ntr)
{
    std::vector<char> live((size_t) ntr, 1);
    if (!mdio_pipe_read_trace_mask(child, live))
        sf_error("Failed to read child trace_mask");
    if ((long) live.size() != ntr) {
        if (child.variables.contains_key("trace_mask"))
            sf_error("trace_mask length %ld != ntr %ld",
                     (long) live.size(), ntr);
        live.assign((size_t) ntr, 1);
    }
    return live;
}

/* Stream the RSF amplitudes over the child in chunk-aligned blocks, applying
   dead-trace fill and accumulating statsV1 in the same pass.  Returns the
   statsV1 JSON. */
static std::string stream_amplitude_to_child(
    mdio::Dataset& child, const std::string& store_uri,
    const std::string& datavar, const std::vector<MdioAxis>& axes,
    const std::vector<long>& sizes, sf_file in, long ns, long ntr, bool verb)
{
    const int rank = (int) axes.size();

    int panel_cap;
    if (!sf_getint("panel", &panel_cap)) panel_cap = 64;
    /* max traces per streamed block along the fastest spatial axis; fallback
       sizing only, used when no whole-chunk block fits blockmb= */

    int blockmb;
    if (!sf_getint("blockmb", &blockmb)) blockmb = 128;
    /* streaming block budget in MiB, and roughly the peak amplitude buffer;
       blocks span whole Zarr chunks, so a larger budget buys more concurrent
       chunk compression */
    static_assert(MDIO_PIPE_BLOCK_MB == 128,
                  "sfdoc needs a literal default; keep it equal to the macro");

    MdioBlockLoop loop;
    std::string note;
    mdio_pipe_resolve_block_loop(store_uri, datavar, sizes,
                                 panel_cap, blockmb, loop, note);
    if (verb && !note.empty()) sf_warning("%s", note.c_str());

    std::vector<char> live = child_trace_mask(child, ntr);

    /* Hoist the amplitude Variable once for the whole block loop. */
    MdioFloat32Var amp_var;
    if (!mdio_get_float32_var(child, datavar, amp_var))
        sf_error("Cannot open float32 variable \"%s\" on child",
                 datavar.c_str());

    /* The sample axis is always faster than the pivot, so one pivot step is a
       whole number of traces (0 in the degenerate single-axis case, which has
       exactly one block starting at trace 0). */
    const long tail_traces = (rank >= 2) ? (loop.tail_floats / ns) : 0;

    std::vector<float> pbuf((size_t) loop.block_floats);
    std::vector<long> oidx(loop.outer_shape.size(), 0);
    MdioStatsAcc stats;

    if (verb)
        sf_warning("Streaming write %ld traces in blocks of %ld along \"%s\" "
                   "(%ld floats/block)",
                   ntr, loop.pivot_block,
                   axes[(size_t) loop.pivot].label.c_str(), loop.block_floats);

    for (long o = 0; o < loop.n_outer; o++) {
        if (loop.pivot > 0)
            mdio_pipe_decode_index(o, loop.pivot, &loop.outer_shape[0],
                                   &oidx[0]);

        for (long ps = 0; ps < loop.n_pivot; ps += loop.pivot_block) {
            long pc = loop.pivot_block;
            if (ps + pc > loop.n_pivot) pc = loop.n_pivot - ps;
            const long block_n = pc * loop.tail_floats;

            sf_floatread(pbuf.data(), block_n, in);
            mdio_pipe_apply_dead_fill(pbuf.data(),
                                      (o * loop.n_pivot + ps) * tail_traces,
                                      block_n / ns, ns, live);

            std::vector<mdio::RangeDescriptor<mdio::Index> > pslices;
            pslices.reserve((size_t) rank);
            for (int a = 0; a < rank; a++) {
                long start, stop;
                if (a < loop.pivot) {
                    start = oidx[(size_t) a];
                    stop  = start + 1;
                } else if (a == loop.pivot) {
                    start = ps;
                    stop  = ps + pc;
                } else {
                    start = 0;
                    stop  = sizes[(size_t) a];
                }
                mdio::RangeDescriptor<mdio::Index> r = {
                    axes[(size_t) a].label.c_str(),
                    (mdio::Index) start,
                    (mdio::Index) stop,
                    1};
                pslices.push_back(r);
            }

            if (!mdio_write_float_block(amp_var, pslices, pbuf.data(),
                                        block_n))
                sf_error("Failed to write data block (var \"%s\", %ld floats)",
                         datavar.c_str(), block_n);

            stats.update(pbuf.data(), block_n);
        }
    }

    return stats.finalize_json();
}

/* Resolve the processing= provenance entry, synthesizing one when the geometry
   changed and the caller omitted it.  Returns NULL when there is nothing to
   append beyond createdOn. */
static nlohmann::json* resolve_processing(const MdioSelection& sel,
                                          nlohmann::json& proc)
{
    char* registered = sf_getstring("processing");
    /* provenance JSON object, or path to a .json file, appended under
       attributes.processing */
    if (NULL != registered) free(registered);

    /* sf_getstring yields only the text between the first and last quote,
       which mangles inline JSON, so take the raw parameter instead.  The call
       above is what makes sfdoc list processing= at all. */
    const char* raw = sf_simtab_get(sf_getpars(), "processing");

    if (NULL != raw && raw[0] != '\0') {
        std::string text = raw;
        std::ifstream fs(text);
        if (fs.good()) {
            std::ostringstream oss;
            oss << fs.rdbuf();
            text = oss.str();
        }
        try {
            proc = nlohmann::json::parse(text);
        } catch (const std::exception& e) {
            sf_error("processing= is not valid JSON (inline or file): %s",
                     e.what());
        }
        return &proc;
    }

    if (!sel.same_geometry) {
        proc = nlohmann::json::object({
            {"op", sel.contract_name},
            {"selection", mdio_pipe_selection_encode(sel)}
        });
        return &proc;
    }
    return NULL;
}

/* Everything between the last amplitude write and the publish: statsV1,
   sample-geometry headers, zarr metadata repair, and provenance. */
static void finalize_child(mdio::Dataset& child, const std::string& tmp_path,
                           const MdioPipeContext& ctx,
                           const std::string& datavar,
                           const std::vector<MdioAxis>& parent_axes,
                           const MdioSelection& sel,
                           const std::string& stats_json, bool verb)
{
    std::string err;

    if (!mdio_pipe_restamp_stats(child, datavar, stats_json, err))
        sf_error("statsV1 restamp failed: %s", err.c_str());
    if (verb) sf_warning("statsV1 restamped");

    if (sel.sample_changed) {
        if (!mdio_pipe_stamp_sample_geometry(child, tmp_path, ctx.source,
                                             datavar, parent_axes, sel, err))
            sf_error("sample-geometry stamp failed (child=%s parent=%s "
                     "datavar=%s contract=%s selection=%s): %s",
                     tmp_path.c_str(), ctx.source.c_str(),
                     datavar.c_str(), sel.contract_name.c_str(),
                     mdio_pipe_selection_encode(sel).c_str(), err.c_str());
        if (verb) sf_warning("sample-geometry headers stamped");
    }

    /* CommitMetadata / SelectField may rewrite structured zarr.json; repair
       before publish so Python MDIO can open the child. */
    if (!sel.same_geometry &&
        !mdio_pipe_repair_zarr_metadata(tmp_path, ctx.source, err,
                                        /*cap_chunks=*/true))
        sf_error("zarr metadata repair failed: %s", err.c_str());

    nlohmann::json proc;
    nlohmann::json* proc_ptr = resolve_processing(sel, proc);
    if (!mdio_pipe_patch_provenance(tmp_path, proc_ptr, NULL, err))
        sf_error("provenance patch failed: %s", err.c_str());
}

/* The faithful finalizer: clone (or reshape) the parent, stream the amplitude
   over it, restamp the derived metadata, then publish atomically.  Fails
   closed — TmpStoreGuard removes the temp store unless every stage passed. */
static void write_faithful_child(sf_file in, const char* out_path,
                                 const std::string& datavar,
                                 const MdioPipeContext& ctx,
                                 mdio::Dataset& parent,
                                 const std::vector<MdioAxis>& parent_axes,
                                 const MdioSelection& sel, bool verb)
{
    const std::string tmp_path =
        std::string(out_path) + ".tmp." + std::to_string((long) getpid());

    build_child_store(ctx, sel, datavar, tmp_path, parent);
    TmpStoreGuard tmp_guard(tmp_path);

    auto cf = mdio::Dataset::Open(tmp_path, mdio::constants::kOpen);
    if (!cf.status().ok())
        sf_error("Cannot open cloned child \"%s\": %s",
                 tmp_path.c_str(), cf.status().ToString().c_str());
    mdio::Dataset child = cf.value();

    /* The RSF stream must hold exactly the child the selection implies. */
    long ns, ntr, total;
    int  rsf_dim;
    rsf_stream_dims(in, ns, ntr, total, rsf_dim);

    std::vector<long> child_sizes = mdio_pipe_child_sizes(sel);
    long ctotal = 1;
    for (size_t a = 0; a < parent_axes.size(); a++) ctotal *= child_sizes[a];
    if (ctotal != total)
        sf_error("RSF sample count %ld != child product %ld", total, ctotal);

    std::string stats_json =
        stream_amplitude_to_child(child, tmp_path, datavar, parent_axes,
                                  child_sizes, in, ns, ntr, verb);

    finalize_child(child, tmp_path, ctx, datavar, parent_axes, sel,
                   stats_json, verb);

    tmp_guard.commit();
    std::string err;
    if (!mdio_pipe_publish(tmp_path, std::string(out_path), err)) {
        mdio_pipe_remove_path(tmp_path);
        sf_error("publish failed: %s", err.c_str());
    }
    if (verb) sf_warning("Faithful child published at \"%s\"", out_path);
}

/* ====================================================================
   Legacy reduced-schema create path (reduced=y).
   ==================================================================== */

/* Which headers hdrcopy= asks to take from the headers=/like= source. */
struct HeaderCopy { bool text, binary, trace; };

static HeaderCopy get_header_copy(void)
{
    const char* hdrcopy = sf_getstring("hdrcopy");
    /* what to copy from headers=/like=: text, binary, trace, or all */
    if (NULL == hdrcopy) hdrcopy = "all";

    const std::string hc = hdrcopy;
    const bool all = hc.find("all") != std::string::npos;
    HeaderCopy c;
    c.text   = all || hc.find("text")   != std::string::npos;
    c.binary = all || hc.find("binary") != std::string::npos;
    c.trace  = all || hc.find("trace")  != std::string::npos;
    return c;
}

/* MDIO axes (slowest first, sample axis last) from the RSF stream header —
   the reverse of RSF axis order. */
static std::vector<MdioAxis> rsf_axes_from_stream(sf_file in, int dim,
                                                  const off_t* n)
{
    std::vector<MdioAxis> axes((size_t) dim);
    std::vector<std::string> used;
    for (int a = 0; a < dim; a++) {
        const int ri = dim - a;          /* RSF axis number (1-based) */
        char key[8];
        MdioAxis ax;
        ax.size = (long) n[ri - 1];
        float fo = 0., fd = 1.;
        snprintf(key, sizeof(key), "o%d", ri);
        ax.o = sf_histfloat(in, key, &fo) ? (double) fo : 0.;
        snprintf(key, sizeof(key), "d%d", ri);
        ax.d = sf_histfloat(in, key, &fd) ? (double) fd : 1.;
        snprintf(key, sizeof(key), "label%d", ri);
        char* lab = sf_histstring(in, key);
        ax.label = clean_name(lab, ri, used);
        if (lab) free(lab);
        snprintf(key, sizeof(key), "unit%d", ri);
        char* un = sf_histstring(in, key);
        if (un) { ax.unit = un; free(un); }
        ax.sample = (ri == 1);
        axes[(size_t) a] = ax;
    }
    return axes;
}

/* Per-key trace-header columns from tfile=, keeping only keys that carry a
   nonzero value somewhere. */
static void read_tfile_columns(long ntr,
                               std::vector<std::string>& keynames,
                               std::vector<std::vector<int> >& columns)
{
    char* tname = sf_getstring("tfile");
    /* input trace header file (from sfsegyread or sfsegyheader) */
    if (NULL == tname) {
        segy_init(SF_NKEYS, NULL);
        return;
    }

    sf_file hin = sf_input("tfile");
    int nk;
    if (!sf_histint(hin, "n1", &nk)) sf_error("No n1= in tfile");
    const long thtr = sf_leftsize(hin, 1);
    if (thtr != ntr)
        sf_warning("tfile has %ld traces, data has %ld", thtr, ntr);
    const long use = (thtr < ntr) ? thtr : ntr;

    segy_init(nk, hin);
    int* all = sf_intalloc(nk * (int) use);
    sf_intread(all, (size_t) nk * use, hin);

    for (int k = 0; k < nk && k < SF_NKEYS; k++) {
        bool nonzero = false;
        std::vector<int> col((size_t) ntr, 0);
        for (long t = 0; t < use; t++) {
            const int v = all[(size_t) t * nk + k];
            col[(size_t) t] = v;
            if (v) nonzero = true;
        }
        if (nonzero) {
            keynames.push_back(segykeyword(k));
            columns.push_back(col);
        }
    }
    free(all);
    sf_fileclose(hin);
}

/* Replace the tfile columns with values read from the header-copy source. */
static void copy_header_columns(mdio::Dataset& srcds, const char* srcpath,
                                long ntr,
                                std::vector<std::string>& keynames,
                                std::vector<std::vector<int> >& columns)
{
    const std::string headersVar = mdio_headers_variable(srcds, NULL);
    const std::vector<mdio::RangeDescriptor<mdio::Index> > nosl;

    keynames.clear();
    columns.clear();
    for (int k = 0; k < SF_NKEYS; k++) {
        const char* nm = segykeyword(k);
        std::vector<double> vals;
        if (!mdio_read_header_key(srcds, std::string(srcpath), nosl,
                                  headersVar, nm, vals)) continue;
        if ((long) vals.size() != ntr) continue;

        std::vector<int> col((size_t) ntr, 0);
        for (long t = 0; t < ntr; t++)
            col[(size_t) t] = (int) vals[(size_t) t];
        keynames.push_back(nm);
        columns.push_back(col);
    }
}

/* Force delrt to reflect the sample-axis origin (mirrors sfsegywrite). */
static void force_delrt(double o1, long ntr,
                        std::vector<std::string>& keynames,
                        std::vector<std::vector<int> >& columns)
{
    const int dk = segykey("delrt");
    const int v = (int)(1000. * o1);

    for (size_t i = 0; i < keynames.size(); i++)
        if (segykey(keynames[i].c_str()) == dk) {
            for (long t = 0; t < ntr; t++) columns[i][(size_t) t] = v;
            return;
        }

    if (o1 != 0.) {
        keynames.push_back("delrt");
        columns.push_back(std::vector<int>((size_t) ntr, v));
    }
}

/* Per-axis chunk sizes for the data variable.  axes[] is MDIO order, so RSF
   axis ri maps to MDIO axis dim-ri. */
static std::vector<long> resolve_chunks(const std::vector<MdioAxis>& axes)
{
    const int dim = (int) axes.size();

    char* chunkstr = sf_getstring("chunk");
    /* chunk size or named strategy for the data variable: an integer applied to
       every axis (0 or >= axis size means a single full-length chunk), or one of
       auto (default, min(size,128) per axis), full (whole array in one chunk), or
       trace (chunk only the sample axis, full-length trace axes).  Per-axis
       overrides chunk1=,chunk2=,... (RSF axis order; axis 1 = samples) take
       precedence. */
    const std::string strat = (NULL != chunkstr) ? chunkstr : "auto";

    const bool named = (strat == "auto" || strat == "full" || strat == "trace");
    const long fixed = named ? 0 : atol(strat.c_str());

    std::vector<long> chunk((size_t) dim);
    for (int a = 0; a < dim; a++) {
        const long size = axes[(size_t) a].size;
        if (strat == "full")
            chunk[(size_t) a] = size;
        else if (strat == "trace")
            chunk[(size_t) a] = axes[(size_t) a].sample ?
                                ((size < 128) ? size : 128) : size;
        else if (!named)
            chunk[(size_t) a] = (fixed <= 0 || fixed >= size) ? size : fixed;
        else
            chunk[(size_t) a] = (size < 128) ? size : 128;
    }

    for (int ri = 1; ri <= dim; ri++) {
        char key[16];
        int v;
        snprintf(key, sizeof(key), "chunk%d", ri);
        if (sf_getint(key, &v)) {
            const size_t a = (size_t) (dim - ri);
            chunk[a] = (v <= 0 || v >= axes[a].size) ? axes[a].size : v;
        }
    }
    return chunk;
}

/* Text and binary SEG-Y headers into the schema metadata, preferring the
   header-copy source over hfile=/bfile=. */
static void put_segy_headers(nlohmann::json& schema, mdio::Dataset* srcds,
                             const HeaderCopy& cp)
{
    char ahead[SF_EBCBYTES];
    bool have_text = (NULL != srcds && cp.text &&
                      mdio_get_text_header(*srcds, ahead));
    if (!have_text) {
        char* hname = sf_getstring("hfile");
        /* input SEG-Y EBCDIC text header file */
        if (NULL != hname) {
            FILE* fp = fopen(hname, "r");
            if (NULL != fp) {
                memset(ahead, ' ', SF_EBCBYTES);
                fread(ahead, 1, SF_EBCBYTES, fp);
                fclose(fp);
                have_text = true;
            }
        }
    }
    if (have_text) mdio_put_text_header(schema, ahead);

    char bhead[SF_BNYBYTES];
    bool have_bin = (NULL != srcds && cp.binary &&
                     mdio_get_binary_header(*srcds, bhead));
    if (!have_bin) {
        char* bname = sf_getstring("bfile");
        /* input SEG-Y binary header file */
        if (NULL != bname) {
            FILE* fp = fopen(bname, "rb");
            if (NULL != fp) {
                memset(bhead, 0, SF_BNYBYTES);
                fread(bhead, 1, SF_BNYBYTES, fp);
                fclose(fp);
                have_bin = true;
            }
        }
    }
    if (have_bin) mdio_put_binary_header(schema, bhead);
}

/* Uniform dimension coordinates (o + i*d). */
static void write_dim_coords(mdio::Dataset& ds,
                             const std::vector<MdioAxis>& axes)
{
    for (size_t a = 0; a < axes.size(); a++) {
        auto vr = ds.variables.get<mdio::dtypes::float32_t>(axes[a].label);
        if (!vr.status().ok()) continue;
        auto var = vr.value();
        auto vdr = mdio::from_variable<mdio::dtypes::float32_t>(var);
        if (!vdr.status().ok()) continue;
        auto vd = vdr.value();
        const long nc = (long) vd.num_samples();
        auto off = vd.get_flattened_offset();
        float* p = static_cast<float*>(vd.get_data_accessor().data());
        for (long i = 0; i < nc; i++)
            p[off + i] = (float)(axes[a].o + i * axes[a].d);
        /* Await the write: an un-awaited TensorStore future may never flush,
           leaving the coordinate at its NaN fill value. */
        if (!var.Write(vd).status().ok())
            sf_warning("Could not write coordinate \"%s\"",
                       axes[a].label.c_str());
    }
}

/* Build a fresh MDIO from the RSF stream alone.  Lossy: it cannot reproduce
   the parent's schema, so the child will not pass the Axis A/B fidelity
   checks.  Kept for smoke tests behind reduced=y. */
static void write_reduced_child(sf_file in, const char* out_path,
                                const std::string& datavar,
                                const DatasetParams& dsp, bool verb)
{
    sf_warning("************************************************************");
    sf_warning("* reduced=y — NOT a lossless write.                       *");
    sf_warning("* Child will NOT pass mdio_conv Axis A/B fidelity checks. *");
    sf_warning("************************************************************");

    const char* srcpath = dsp.header_source();
    const HeaderCopy cp = get_header_copy();

    off_t n[SF_MAX_DIM];
    const int dim = sf_largefiledims(in, n);
    const long ns = (long) n[0];
    long total = 1;
    for (int i = 0; i < dim; i++) total *= (long) n[i];
    const long ntr = total / ns;

    std::vector<MdioAxis> axes = rsf_axes_from_stream(in, dim, n);
    const double o1 = axes[(size_t) (dim - 1)].o;  /* sample origin -> delrt */

    std::vector<std::string> keynames;
    std::vector<std::vector<int> > columns;
    read_tfile_columns(ntr, keynames, columns);

    mdio::Dataset* srcds = NULL;
    if (NULL != srcpath) {
        auto sf = mdio::Dataset::Open(std::string(srcpath),
                                      mdio::constants::kOpen);
        if (!sf.status().ok())
            sf_error("Cannot open headers source \"%s\"", srcpath);
        srcds = new mdio::Dataset(sf.value());
        if (cp.trace)
            copy_header_columns(*srcds, srcpath, ntr, keynames, columns);
    }

    force_delrt(o1, ntr, keynames, columns);

    char* dname = sf_getstring("name");
    /* dataset name recorded in the MDIO metadata */
    const std::string dsname = (NULL != dname) ? dname : "Madagascar MDIO";

    nlohmann::json schema = mdio_build_schema(dsname, axes, datavar,
                                             keynames, resolve_chunks(axes));
    put_segy_headers(schema, srcds, cp);
    delete srcds;

    if (verb) sf_warning("Creating MDIO \"%s\" (%ld samples x %ld traces)",
                         out_path, ns, ntr);

    auto dsFut = mdio::Dataset::from_json(schema, std::string(out_path),
                                          mdio::constants::kCreate);
    if (!dsFut.status().ok())
        sf_error("Failed to create MDIO \"%s\": %s", out_path,
                 dsFut.status().ToString().c_str());
    mdio::Dataset ds = dsFut.value();

    write_dim_coords(ds, axes);

    /* RSF order matches MDIO row-major order. */
    std::vector<float> buf((size_t) total);
    sf_floatread(buf.data(), total, in);
    if (!write_float_buf(ds, datavar, buf.data(), total))
        sf_error("Failed to write data variable \"%s\"", datavar.c_str());

    for (size_t i = 0; i < keynames.size(); i++)
        if (!write_int_var(ds, keynames[i], columns[i]) && verb)
            sf_warning("Could not write header variable \"%s\"",
                       keynames[i].c_str());

    if (verb) sf_warning("Done");
}

int main(int argc, char* argv[])
{
    sf_init(argc, argv);

    bool verb;
    if (!sf_getbool("verb", &verb)) verb = false;
    /* Verbosity flag */

    sf_file in = sf_input("in");
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    char* path = sf_getstring("mdio");
    /* output MDIO dataset (path or gs://, s3:// URL) */
    if (NULL == path) sf_error("Need mdio=");

    char* dataname = sf_getstring("data");
    /* name of the MDIO data variable (default "seismic") */
    std::string datavar = (NULL != dataname) ? dataname : "seismic";

    DatasetParams dsp = get_dataset_params();

    /* Opt-in lossy reduced-schema path (smoke tests only).  Default is the
       faithful parent-clone finalizer. */
    bool reduced;
    if (!sf_getbool("reduced", &reduced)) reduced = false;
    /* acknowledge reduced-schema write (lossy; not Axis A/B faithful) */

    /* ---- pipe context: require or recover, then verify ---- */
    MdioPipeContext ctx;
    mdio_pipe_read_context(in, ctx);

    if (!ctx.present) {
        const char* source = dsp.recovery_source();
        if (NULL == source) {
            /* No parent MDIO to be faithful to.  The reduced path can still
               originate a fresh (lossy) store from the RSF stream alone; the
               faithful path cannot, so it stays fail-closed. */
            if (reduced) {
                write_reduced_child(in, path, datavar, dsp, verb);
                exit(0);
            }
            sf_error("Missing mdio_* pipe context on input "
                     "(need sfmdioread upstream, or like=/context= recovery; "
                     "or pass reduced=y to originate a new lossy MDIO from "
                     "RSF)");
        }
        recover_pipe_context(source, dataname, ctx, verb);
    }
    validate_pipe_context(ctx);

    /* Prefer the stamped data variable over the CLI default. */
    if (NULL == dataname)
        datavar = ctx.data_variable;
    else if (datavar != ctx.data_variable)
        sf_warning("data=%s overrides mdio_data_variable=%s",
                   datavar.c_str(), ctx.data_variable.c_str());

    mdio::Dataset parent = open_verified_parent(ctx, verb);
    std::vector<MdioAxis> parent_axes = mdio_axes(parent, ctx.data_variable);
    MdioSelection sel = resolve_selection(in, ctx, parent_axes, verb);

    if (!reduced) {
        write_faithful_child(in, path, datavar, ctx, parent, parent_axes,
                             sel, verb);
        exit(0);
    }

    write_reduced_child(in, path, datavar, dsp, verb);
    exit(0);
}
