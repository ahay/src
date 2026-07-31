/* Shared helpers for sfmdioread -> sf* -> sfmdiowrite.
   Reader stamps mdio_* keys; writer verifies, clones parent, overwrites data. */

#ifndef _mdio_pipe_hh
#define _mdio_pipe_hh

#include <string>
#include <vector>

#include "mdio2segy.hh"

/* RSF history keys stamped by sfmdioread, required by sfmdiowrite. */
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

/* Pipe context on RSF (or recovered via like=/context=). */
struct MdioPipeContext {
    std::string source;        /* parent MDIO URI                           */
    std::string fingerprint;   /* canonical parent descriptor digest         */
    std::string data_variable;
    std::string dtype;
    std::string contract;
    std::string version;
    std::string selection;     /* optional "label:start:step:count;..."     */
    bool        present;       /* mdio_source was found                     */
};

enum MdioAxisSelKind {
    MDIO_SEL_PRESERVED = 0,
    MDIO_SEL_SUBSET    = 1,  /* contiguous crop, step == 1                  */
    MDIO_SEL_DECIMATED = 2,  /* integer stride > 1                          */
    MDIO_SEL_RESAMPLED = 3   /* non-integer regular sample grid             */
};

struct MdioAxisSel {
    std::string     label;
    long            start;   /* parent index (0 if RESAMPLED)               */
    long            step;    /* >= 1 (0 if RESAMPLED)                       */
    long            count;   /* child axis length                           */
    MdioAxisSelKind kind;
    double          o;       /* child origin                                */
    double          d;       /* child sampling                              */
};

struct MdioSelection {
    std::vector<MdioAxisSel> axes; /* MDIO order, parent rank                 */
    bool same_geometry;            /* every axis preserved                  */
    bool sample_changed;
    bool spatial_changed;
    bool exact;                    /* all axes exact integer selections     */
    std::string contract_name;
};

/* Metadata-only parent digest (name, var, dtype, axes). Use observed dtype. */
std::string mdio_pipe_fingerprint(mdio::Dataset& ds,
                                  const std::string& datavar,
                                  const std::string& dtype);

/* Stamp mdio_* keys before sf_fileflush. */
void mdio_pipe_stamp(sf_file out,
                     const std::string& source,
                     const std::string& fingerprint,
                     const std::string& datavar,
                     const std::string& dtype,
                     const std::string& contract = MDIO_PIPE_CONTRACT_SAME_GEOM,
                     const std::string& version  = MDIO_PIPE_CONTEXT_VERSION,
                     const std::string& selection = std::string());

/* Read mdio_* from RSF; ctx.present if mdio_source set. */
void mdio_pipe_read_context(sf_file in, MdioPipeContext& ctx);

/* Fill ctx from open parent (like=/context= recovery). */
bool mdio_pipe_context_from_parent(mdio::Dataset& ds,
                                   const std::string& parent_uri,
                                   const std::string& datavar,
                                   MdioPipeContext& ctx);

/* Remove local tree or remote prefix; missing ok. */
void mdio_pipe_remove_path(const std::string& path);

/* Publish tmp -> final (rename local; copy+delete remote). */
bool mdio_pipe_publish(const std::string& tmp_path,
                       const std::string& final_path,
                       std::string& err);

/* Clone parent except skip_var chunk payloads (metadata kept). */
bool mdio_clone_store_skip_var(const std::string& parent_uri,
                               const std::string& child_path,
                               const std::string& skip_var,
                               std::string& err);

/* Read trace_mask into live[] (1=live); absent -> all-live. */
bool mdio_pipe_read_trace_mask(mdio::Dataset& ds,
                               std::vector<char>& live,
                               const std::string& mask_name = "trace_mask");

/* NaN-fill dead traces in RSF-order buffer; first_trace is global index. */
void mdio_pipe_apply_dead_fill(float* buf, long first_trace, long ntr, long ns,
                               const std::vector<char>& live);

/* Streamed stats (mdio_compat / Axis A); empty histogram. */
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

/* Stamp statsV1 on data var and CommitMetadata(). */
bool mdio_pipe_restamp_stats(mdio::Dataset& ds,
                             const std::string& datavar,
                             const std::string& stats_json,
                             std::string& err);

/* Patch root zarr.json: createdOn + append attributes.processing. */
bool mdio_pipe_patch_provenance(const std::string& child_path,
                                const nlohmann::json* processing,
                                const char* created_on,
                                std::string& err);

/* Encode selection as "label:start:step:count;...". */
std::string mdio_pipe_selection_encode(const MdioSelection& sel);

/* Parse mdio_sel=; parent_axes supply labels. */
bool mdio_pipe_selection_parse(const std::string& text,
                               const std::vector<MdioAxis>& parent_axes,
                               MdioSelection& sel,
                               std::string& err);

/* Align RSF o#/d#/n# to parent; sample axis may resample (sfremap1). */
bool mdio_pipe_detect_selection(sf_file in,
                                const std::vector<MdioAxis>& parent_axes,
                                MdioSelection& sel,
                                std::string& err);

/* Build selection from reader window arrays (MDIO order). */
bool mdio_pipe_selection_from_window(const std::vector<MdioAxis>& parent_axes,
                                     const std::vector<long>& f,
                                     const std::vector<long>& jp,
                                     const std::vector<long>& m,
                                     MdioSelection& sel,
                                     std::string& err);

/* Child axis sizes after sel (MDIO order). */
std::vector<long> mdio_pipe_child_sizes(const MdioSelection& sel);

/* Clone parent to tmp_path, trim vars to selection; data samples for caller. */
bool mdio_pipe_build_child_store(mdio::Dataset& parent,
                                 const std::string& parent_uri,
                                 const MdioSelection& sel,
                                 const std::string& datavar,
                                 const std::string& tmp_path,
                                 std::string& err);

/* Repair zarr.json for Python MDIO (dtype spelling, cap chunk_shape). */
bool mdio_pipe_repair_zarr_metadata(const std::string& child_path,
                                    const std::string& parent_uri,
                                    std::string& err,
                                    bool cap_chunks = true);

/* Restamp SEG-Y sample geometry after sample-axis change. */
bool mdio_pipe_stamp_sample_geometry(mdio::Dataset& child,
                                     const std::string& child_path,
                                     const std::string& parent_uri,
                                     const std::string& datavar,
                                     const std::vector<MdioAxis>& parent_axes,
                                     const MdioSelection& sel,
                                     std::string& err);

/* Default streaming block budget (MiB). */
#define MDIO_PIPE_BLOCK_MB 128

/* Chunk-aligned block loop: one contiguous RSF read/write per block. */
struct MdioBlockLoop {
    int  pivot;         /* blocked axis                                  */
    long pivot_block;   /* pivot indices per block                       */
    long n_pivot;       /* pivot extent                                  */
    long n_outer;       /* product of slower axes                        */
    long tail_floats;   /* product of faster axes                        */
    long block_floats;  /* pivot_block * tail_floats                     */
    std::vector<long> outer_shape;  /* slower-axis extents                 */
    bool aligned;       /* false: straddles chunks (slow)                */
};

/* Plan block loop for sizes vs chunk grid; note explains fallback. */
void mdio_pipe_resolve_block_loop(const std::string& store_uri,
                                  const std::string& datavar,
                                  const std::vector<long>& sizes,
                                  int blockmb,
                                  MdioBlockLoop& loop,
                                  std::string& note);

/* Row-major linear index -> multi-index. */
void mdio_pipe_decode_index(long lin, int ndim, const long* shape, long* idx);

#endif
