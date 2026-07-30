/* Shared pipe-context helpers for sfmdioread / sfmdiowrite.
   Stamps and verifies the mdio_* RSF history keys that carry parent
   provenance across a Madagascar filter pipe, plus the finalizer helpers
   (clone/publish, geometry selection, dead-fill, streaming stats, provenance,
   chunk-aligned panels, selective amplitude-skip clone). */

#ifndef _mdio_pipe_hh
#define _mdio_pipe_hh

#include <string>
#include <vector>

#include "mdio2segy.hh"

/* RSF history key names stamped by sfmdioread and required by sfmdiowrite. */
#define MDIO_PIPE_KEY_SOURCE       "mdio_source"
#define MDIO_PIPE_KEY_FINGERPRINT  "mdio_source_fingerprint"
#define MDIO_PIPE_KEY_DATA_VAR     "mdio_data_variable"
#define MDIO_PIPE_KEY_DTYPE        "mdio_data_dtype"
#define MDIO_PIPE_KEY_CONTRACT     "mdio_contract"
#define MDIO_PIPE_KEY_VERSION      "mdio_context_version"
#define MDIO_PIPE_KEY_SELECTION    "mdio_sel"

#define MDIO_PIPE_CONTEXT_VERSION  "1"
#define MDIO_PIPE_CONTRACT_SAME_GEOM "same-geometry"
#define MDIO_PIPE_CONTRACT_TRUNCATE  "truncate-exact"
#define MDIO_PIPE_CONTRACT_SLICE     "slice-exact"
#define MDIO_PIPE_CONTRACT_DECIMATE  "slice-decimated"
#define MDIO_PIPE_CONTRACT_RESAMPLE  "resample"
#define MDIO_PIPE_DTYPE_FLOAT32    "float32"

/* Parsed pipe context carried on an RSF stream (or recovered via like=/context=). */
struct MdioPipeContext {
    std::string source;        /* parent MDIO URI                          */
    std::string fingerprint;   /* "sha256:<hex>" of canonical parent desc  */
    std::string data_variable; /* selected float data variable             */
    std::string dtype;         /* expected "float32" for Phase 1           */
    std::string contract;      /* e.g. "same-geometry"                     */
    std::string version;       /* context schema version ("1")             */
    std::string selection;     /* optional "label:start:step:count;..."    */
    bool        present;       /* true if at least mdio_source was found   */
};

/* Per-axis exact index selection (Phase 4 geometry adapters). */
enum MdioAxisSelKind {
    MDIO_SEL_PRESERVED = 0,
    MDIO_SEL_SUBSET    = 1,  /* contiguous crop, step == 1, shrunk        */
    MDIO_SEL_DECIMATED = 2,  /* integer stride > 1                        */
    MDIO_SEL_RESAMPLED = 3   /* non-integer regular sample grid (o/d)     */
};

struct MdioAxisSel {
    std::string     label;
    long            start;   /* inclusive parent index (unused if RESAMPLED) */
    long            step;    /* >= 1 (unused / 0 if RESAMPLED)               */
    long            count;   /* child length along this axis                 */
    MdioAxisSelKind kind;
    double          o;       /* child origin (always set; stamp + resample)  */
    double          d;       /* child spacing (always set; stamp + resample) */
};

struct MdioSelection {
    std::vector<MdioAxisSel> axes; /* MDIO order, same rank as parent_axes */
    bool same_geometry;            /* every axis preserved                 */
    bool sample_changed;           /* last MDIO axis not preserved         */
    bool spatial_changed;          /* any non-sample axis not preserved    */
    bool exact;                    /* all axes are exact integer selections*/
    std::string contract_name;     /* fidelity contract preset name        */
};

/* Build "sha256:<hex>" over a canonical parent descriptor:
     name\\n datavar\\n dtype\\n
     for each axis: label\\tsize\\to\\td\\tunit\\n
   Digests metadata only — never amplitude bytes.  Reader and writer must
   call this on the same (unwindowed) parent Dataset.  Pass the dtype observed
   via mdio_variable_dtype, not a literal: a constant here would drop out of
   the digest and leave the dtype unbound. */
std::string mdio_pipe_fingerprint(mdio::Dataset& ds,
                                  const std::string& datavar,
                                  const std::string& dtype);

/* Stamp all mdio_* keys onto an RSF output before sf_fileflush.
   selection may be empty (same-geometry / unknown). */
void mdio_pipe_stamp(sf_file out,
                     const std::string& source,
                     const std::string& fingerprint,
                     const std::string& datavar,
                     const std::string& dtype,
                     const std::string& contract = MDIO_PIPE_CONTRACT_SAME_GEOM,
                     const std::string& version  = MDIO_PIPE_CONTEXT_VERSION,
                     const std::string& selection = std::string());

/* Read mdio_* keys from an RSF input.  Sets ctx.present when mdio_source is set.
   Missing individual keys leave the corresponding string empty. */
void mdio_pipe_read_context(sf_file in, MdioPipeContext& ctx);

/* Fill ctx from an already-open parent Dataset (like=/context= recovery).
   Sets source to parent_uri; computes fingerprint for datavar. */
bool mdio_pipe_context_from_parent(mdio::Dataset& ds,
                                   const std::string& parent_uri,
                                   const std::string& datavar,
                                   MdioPipeContext& ctx);

/* ---- Phase 2 same-geometry finalizer helpers ---- */

/* Stats zero-mask tolerance (ingest / mdio_compat STATS_ZERO_ATOL). */
#define MDIO_PIPE_STATS_ZERO_ATOL 1e-8

/* True when uri looks like a remote MDIO path (s3://, gs://, ...). */
bool mdio_pipe_is_remote_uri(const std::string& uri);

/* Recursive copy of a parent Zarr store to child_path (local filesystem
   fast-path, or tensorstore kvstore List/Read/Write for s3:// / gs://).
   Returns false on I/O error. */
bool mdio_clone_store(const std::string& parent_uri,
                      const std::string& child_path,
                      std::string& err);

/* Remove a local path (file or directory tree).  Ignores missing paths. */
void mdio_pipe_remove_path(const std::string& path);

/* Atomically publish tmp_path -> final_path.  Local uses rename (or
   copy+remove).  Remote uses kvstore copy then delete of the temp prefix
   (S3/GCS have no atomic directory rename).  Removes any pre-existing
   final_path.  Returns false on I/O error. */
bool mdio_pipe_publish(const std::string& tmp_path,
                       const std::string& final_path,
                       std::string& err);

/* Float dead-trace fill value (NaN; mdio fill_value_map for float). */
float mdio_pipe_dead_fill(void);

/* Read trace_mask into a flat live[] buffer (row-major, 1 = live, 0 = dead).
   Variable name defaults to "trace_mask".  If absent, fills all-live and
   returns true.  Returns false on a present-but-unreadable mask. */
bool mdio_pipe_read_trace_mask(mdio::Dataset& ds,
                               std::vector<char>& live,
                               const std::string& mask_name = "trace_mask");

/* Apply dead fill to every sample of dead traces in a flat RSF buffer
   (ntr traces x ns samples, sample-fastest within each trace — RSF order).

   first_trace is the buffer's global trace index, so one streamed block can be
   masked without materializing the whole volume.  Traces at or past the end of
   live are left untouched. */
void mdio_pipe_apply_dead_fill(float* buf, long first_trace, long ntr, long ns,
                               const std::vector<char>& live);

/* Recompute SummaryStatistics over finite, |x| > STATS_ZERO_ATOL samples
   and return camelCase JSON (empty histogram).  Matches mdio_compat /
   Axis A.  All-zero / all-fill → empty_stats sentinel. */
std::string mdio_pipe_statsv1_json(const float* buf, long total);

/* Incremental stats accumulator for streamed panels.  Same semantics as
   mdio_pipe_statsv1_json / mdio_compat.compute_summary_statistics_lazy. */
struct MdioStatsAcc {
    long   count;
    double sum;
    double sum_sq;
    double vmin;
    double vmax;
    bool   have;

    MdioStatsAcc()
        : count(0), sum(0), sum_sq(0), vmin(0), vmax(0), have(false) {}

    void update(const float* buf, long n);
    std::string finalize_json() const;
};

/* Stamp statsV1 onto the named data variable of an open child Dataset and
   CommitMetadata().  stats_json is the object (or JSON string) to store. */
bool mdio_pipe_restamp_stats(mdio::Dataset& ds,
                             const std::string& datavar,
                             const std::string& stats_json,
                             std::string& err);

/* Patch dataset-level provenance on a child store (local or remote) via
   the root zarr.json / .zattrs kvstore key: set createdOn (UTC now unless
   created_on given) and optionally append one processing entry under
   attributes.processing. */
bool mdio_pipe_patch_provenance(const std::string& child_path,
                                const nlohmann::json* processing,
                                const char* created_on,
                                std::string& err);

/* Same-geometry preflight: RSF dims / o / d must match parent axes
   (MDIO order = reverse of RSF).  Returns false and fills err on mismatch. */
bool mdio_pipe_assert_same_geometry(sf_file in,
                                    const std::vector<MdioAxis>& parent_axes,
                                    std::string& err);

/* ---- Phase 4 geometry adapters ---- */

/* Encode selection as "label:start:step:count;..." for mdio_sel=. */
std::string mdio_pipe_selection_encode(const MdioSelection& sel);

/* Parse mdio_sel= string into sel (axes filled; flags recomputed).
   parent_axes supplies expected labels/order for validation. */
bool mdio_pipe_selection_parse(const std::string& text,
                               const std::vector<MdioAxis>& parent_axes,
                               MdioSelection& sel,
                               std::string& err);

/* Derive an exact integer selection — or, on the sample axis only, a regular
   in-range interpolating resample grid — by aligning RSF o#/d#/n# to parent
   axes.  Returns false (with err) when any spatial axis is not an exact
   regular subsequence, or when a resampled sample grid is non-positive,
   reversed, or extrapolates past the parent (fail closed). */
bool mdio_pipe_detect_selection(sf_file in,
                                const std::vector<MdioAxis>& parent_axes,
                                MdioSelection& sel,
                                std::string& err);

/* Build selection from reader-side windowing arrays (MDIO order). */
bool mdio_pipe_selection_from_window(const std::vector<MdioAxis>& parent_axes,
                                     const std::vector<long>& f,
                                     const std::vector<long>& jp,
                                     const std::vector<long>& m,
                                     MdioSelection& sel,
                                     std::string& err);

/* Finalize selection flags + fidelity contract_name from per-axis kinds. */
void mdio_pipe_selection_finalize(MdioSelection& sel);

/* Build RangeDescriptors for parent.isel / Variable::slice from sel
   (absolute parent indices, half-open with step). */
std::vector<mdio::RangeDescriptor<mdio::Index> >
mdio_pipe_selection_ranges(const MdioSelection& sel);

/* Clone parent → tmp_path, compact non-amplitude vars onto the origin for
   any start!=0 / step!=1 axes, then trim every variable to sel counts.
   Amplitude shape is trimmed but left unfilled (caller streams it).
   datavar names the amplitude variable to skip compacting. */
bool mdio_pipe_build_child_store(mdio::Dataset& parent,
                                 const std::string& parent_uri,
                                 const MdioSelection& sel,
                                 const std::string& datavar,
                                 const std::string& tmp_path,
                                 std::string& err);

/* After Resize / SelectField / CommitMetadata, repair per-variable zarr.json
   so structured dtypes stay Python-MDIO compatible.  Optionally also cap
   chunk_shape ≤ shape (only safe once TensorStore will not write again). */
bool mdio_pipe_repair_zarr_metadata(const std::string& child_path,
                                    const std::string& parent_uri,
                                    std::string& err,
                                    bool cap_chunks = true);

/* Stamp sample-geometry SEG-Y fields after a sample-axis change:
   samples_per_trace / sample_interval on binary + live traces, and
   delay_recording_time shifted by round(o_c - o_p) ms.  No-op (true)
   when sample axis preserved or headers absent. */
bool mdio_pipe_stamp_sample_geometry(mdio::Dataset& child,
                                     const std::string& child_path,
                                     const std::string& parent_uri,
                                     const std::string& datavar,
                                     const std::vector<MdioAxis>& parent_axes,
                                     const MdioSelection& sel,
                                     std::string& err);

/* Child axis sizes after applying sel (MDIO order). */
std::vector<long> mdio_pipe_child_sizes(const MdioSelection& sel);

/* ---- Phase 5 efficiency helpers ---- */

/* Read the amplitude variable's chunk_shape from store_uri/datavar/zarr.json
   (local or remote).  Returns false when unavailable; chunk left unchanged. */
bool mdio_pipe_chunk_shape(const std::string& store_uri,
                           const std::string& datavar,
                           std::vector<long>& chunk,
                           std::string& err);

/* Choose a panel size (trace count along the streamed RSF-n2 / MDIO
   rank-2 axis) that is a multiple of chunk_along when possible, not
   exceeding panel_cap traces and n_fast.  panel_cap is the user panel=
   upper bound (legacy: max traces).  Falls back to panel_cap when
   chunk_along < 1. */
int mdio_pipe_aligned_panel_traces(int panel_cap, long n_fast,
                                   long chunk_along);

/* Default streaming block budget (MiB) for mdio_pipe_plan_blocks.  One
   chunk row along the pivot is usually already enough concurrency for
   TensorStore, so this stays small on purpose. */
#define MDIO_PIPE_BLOCK_MB 128

/* A chunk-aligned streaming block over the MDIO index space.

   The RSF stream is row-major over the MDIO axes (slowest first), so the
   set of samples with a single index on every axis slower than `pivot`, a
   contiguous run along `pivot`, and the full extent of every faster axis
   is one contiguous run of the stream — readable/writable with a single
   sf_floatread / sf_floatwrite.

   Taking the faster axes in full always covers whole Zarr chunks (a
   truncated final chunk is still covered in full).  So the block covers
   only whole chunks as long as the run along `pivot` is a multiple of the
   chunk extent there and every slower axis has an effective chunk extent
   of 1.  That is what keeps TensorStore from read-modify-writing — and
   re-compressing — a chunk once per index of a coarser-chunked axis. */
struct MdioBlockPlan {
    int  pivot;        /* axis blocked along, in [0, rank-2]            */
    long pivot_block;  /* extent along pivot (multiple of its chunk)    */
    long block_floats; /* pivot_block x product of faster axis extents  */
    bool aligned;      /* false: no whole-chunk plan fits the budget    */
};

/* Plan chunk-aligned streaming blocks over `sizes` given the on-disk
   `chunk` grid and a budget in floats.  Both vectors are in MDIO order
   (slowest first, sample axis last) and must have the same length.

   Picks the slowest viable pivot so each block spans as many whole chunks
   as possible (TensorStore compresses those concurrently).  Returns false
   for rank < 2 or malformed input; sets plan.aligned = false when even one
   chunk along the finest viable pivot blows past 8x the budget, in which
   case the caller should keep its own (unaligned) panel sizing. */
bool mdio_pipe_plan_blocks(const std::vector<long>& sizes,
                           const std::vector<long>& chunk,
                           long budget_floats,
                           MdioBlockPlan& plan);

/* A resolved streaming block loop over an MDIO index space.

   MdioBlockPlan says which axis to block along; this adds the loop bounds
   both mains need to walk it: the extent of the pivot axis, the product of
   the slower axes (visited one index at a time), and the product of the
   faster axes (always taken in full, which is what makes a block one
   contiguous run of the RSF stream). */
struct MdioBlockLoop {
    int  pivot;         /* axis blocked along                             */
    long pivot_block;   /* indices of pivot per block                     */
    long n_pivot;       /* extent of the pivot axis                       */
    long n_outer;       /* product of the axes slower than pivot          */
    long tail_floats;   /* product of the axes faster than pivot          */
    long block_floats;  /* pivot_block x tail_floats: the buffer size     */
    std::vector<long> outer_shape;  /* extents of the axes slower than pivot */
    bool aligned;       /* false: unaligned fallback panels (slow)        */
};

/* Resolve the streaming block loop for `sizes` (MDIO order) against the chunk
   grid of store_uri/datavar.

   panel_cap and blockmb are passed in rather than read here so that panel=
   and blockmb= stay visible to each main's own sfdoc scrape.  Falls back to
   unaligned panels along the fastest spatial axis when no whole-chunk block
   fits the budget, leaving loop.aligned false and a diagnostic in `note`.
   Always yields a usable loop. */
void mdio_pipe_resolve_block_loop(const std::string& store_uri,
                                  const std::string& datavar,
                                  const std::vector<long>& sizes,
                                  int panel_cap, int blockmb,
                                  MdioBlockLoop& loop,
                                  std::string& note);

/* Decode a row-major (last axis fastest) linear index into a multi-index. */
void mdio_pipe_decode_index(long lin, int ndim, const long* shape, long* idx);

/* Same-geometry clone that copies metadata + every non-amplitude object
   but skips amplitude chunk payloads (keeps datavar/zarr.json so the
   chunk grid / codec / fill are inherited).  Amplitude is filled in place
   by the streaming writer. */
bool mdio_clone_store_skip_var(const std::string& parent_uri,
                               const std::string& child_path,
                               const std::string& skip_var,
                               std::string& err);

#endif
