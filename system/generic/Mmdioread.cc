/* Convert an MDIO dataset to RSF.

Reads the MDIO (Zarr) format through the mdio-cpp library and writes the
amplitudes to an RSF file, the per-trace SEG-Y headers to a separate tfile
(compatible with sfsegyread/sfsegywrite), and optionally the SEG-Y EBCDIC text
header (hfile) and 400-byte binary header (bfile) when they are present in the
dataset metadata.

The input may be a local path, or a gs:// or s3:// URL (when mdio-cpp was built
with the corresponding cloud drivers).

The data can be windowed on read with the same parameters as sfwindow
(n#, f#, j#, min#, max#).  Axis 1 is the fast (sample) axis; the remaining
axes index the trace grid, in the order RSF stores them (n2 is the fastest
trace dimension).  Because mdio-cpp slices with unit stride, a contiguous
bounding box is read lazily and any j#>1 decimation is applied afterwards.

Trace headers are read either from one MDIO variable per SEG-Y key, or from a
single structured "headers" variable.
*/

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>

#include "mdio2segy.hh"
#include "mdio_pipe.hh"

/* Read window per axis (MDIO order), from sfwindow n#/f#/j#/... params. */
struct ReadWindow {
    std::vector<long> f;     /* first index on the parent axis           */
    std::vector<long> n;     /* output count along the axis              */
    std::vector<long> j;     /* decimation stride                        */
    std::vector<long> span;  /* parent indices spanned: (n-1)*j + 1      */
};

/* One axis of read window (window.c rules); RSF axis i -> MDIO rank-i. */
static void resolve_axis_window(int i, int a, MdioAxis& ax, ReadWindow& win)
{
    const long n = ax.size;
    const double o = ax.o, d0 = ax.d;
    char  key[8];
    int   iv;
    float av;
    off_t lv;

    /* jump j# (or sampling d#) */
    long j = 1;
    snprintf(key, sizeof(key), "j%d", i);
    if (sf_getint(key, &iv)) {
        j = iv;
    } else {
        snprintf(key, sizeof(key), "d%d", i);
        if (sf_getfloat(key, &av) && d0 != 0.) j = (long)(0.5 + av / d0);
    }
    if (j < 1) j = 1;

    /* start f# (or minimum min#) */
    long ff = 0;
    snprintf(key, sizeof(key), "f%d", i);
    if (sf_getlargeint(key, &lv)) {
        ff = lv;
    } else {
        snprintf(key, sizeof(key), "min%d", i);
        if (sf_getfloat(key, &av) && d0 != 0.) ff = (long)(0.5 + (av - o) / d0);
    }
    if (ff < 0) ff = n + ff;
    if (ff < 0) sf_error("Negative f%d", i);

    const double onew = o + ff * d0;
    const double dnew = d0 * j;

    /* count n# (or maximum max#) */
    long mm;
    snprintf(key, sizeof(key), "n%d", i);
    if (sf_getlargeint(key, &lv)) {
        mm = lv;
    } else {
        snprintf(key, sizeof(key), "max%d", i);
        if (sf_getfloat(key, &av) && dnew != 0.)
            mm = (long)(1.5 + (av - onew) / dnew);
        else
            mm = 1 + (n - 1 - ff) / j;
    }
    if (mm < 1) mm = 1;
    if (ff + (mm - 1) * j > n - 1)
        sf_error("n%d=%ld is too big (axis %s, size %ld)",
                 i, mm, ax.label.c_str(), n);

    win.f[(size_t) a]    = ff;
    win.n[(size_t) a]    = mm;
    win.j[(size_t) a]    = j;
    win.span[(size_t) a] = (mm - 1) * j + 1;

    ax.o = onew;   /* carry windowed geometry to the output */
    ax.d = dnew;
}

static ReadWindow resolve_window(std::vector<MdioAxis>& axes)
{
    const int rank = (int) axes.size();
    ReadWindow win;
    win.f.resize((size_t) rank);
    win.n.resize((size_t) rank);
    win.j.resize((size_t) rank);
    win.span.resize((size_t) rank);
    for (int i = 1; i <= rank; i++)
        resolve_axis_window(i, rank - i, axes[(size_t) (rank - i)], win);
    return win;
}

/* Parent bounding box for mdio-cpp slice (unit stride; j# applied later). */
static std::vector<mdio::RangeDescriptor<mdio::Index> >
window_bounding_box(const std::vector<MdioAxis>& axes, const ReadWindow& win)
{
    std::vector<mdio::RangeDescriptor<mdio::Index> > slices;
    for (size_t a = 0; a < axes.size(); a++) {
        mdio::RangeDescriptor<mdio::Index> r = {
            axes[a].label.c_str(),
            (mdio::Index) win.f[a],
            (mdio::Index)(win.f[a] + win.span[a]),
            1};
        slices.push_back(r);
    }
    return slices;
}

/* Open amplitude out with windowed geometry; RSF axis i = MDIO rank-i. */
static sf_file open_amplitude_output(const std::vector<MdioAxis>& axes,
                                     const ReadWindow& win)
{
    const int rank = (int) axes.size();
    sf_file out = sf_output("out");
    sf_settype(out, SF_FLOAT);
    for (int i = 1; i <= rank; i++) {
        const size_t a = (size_t) (rank - i);
        char key[24];   /* widest is "label" + an int, plus NUL */
        snprintf(key, sizeof(key), "n%d", i);
        sf_putint(out, key, (int) win.n[a]);
        snprintf(key, sizeof(key), "d%d", i);
        sf_putfloat(out, key, (float) axes[a].d);
        snprintf(key, sizeof(key), "o%d", i);
        sf_putfloat(out, key, (float) axes[a].o);
        snprintf(key, sizeof(key), "label%d", i);
        sf_putstring(out, key,
                     axes[a].label.empty() ?
                     (i == 1 ? "Time" : "Trace") : axes[a].label.c_str());
        if (!axes[a].unit.empty()) {
            snprintf(key, sizeof(key), "unit%d", i);
            sf_putstring(out, key, axes[a].unit.c_str());
        }
    }
    return out;
}

static sf_file open_header_output(long ntr, sf_file out)
{
    sf_file hdr = sf_output("tfile");
    sf_putint(hdr, "n1", SF_NKEYS);
    sf_putint(hdr, "n2", (int) ntr);
    sf_settype(hdr, SF_INT);
    segy_init(SF_NKEYS, NULL);
    segy2hist(hdr, SF_NKEYS);

    char* tname = sf_getstring("tfile");
    /* output trace header file */
    if (NULL != out)
        sf_putstring(out, "head", (NULL != tname) ? tname : "tfile");
    return hdr;
}

/* Stamp mdio_* + selection before flush (propagates via sf_fileflush). */
static void stamp_pipe_context(sf_file out, const char* path,
                               const std::string& fingerprint,
                               const std::string& datavar,
                               const std::string& dtype,
                               const std::vector<MdioAxis>& axes,
                               const ReadWindow& win, bool verb)
{
    MdioSelection rsel;
    std::string serr, sel_text;
    std::string contract = MDIO_PIPE_CONTRACT_SAME_GEOM;

    if (mdio_pipe_selection_from_window(axes, win.f, win.j, win.n,
                                        rsel, serr)) {
        sel_text = mdio_pipe_selection_encode(rsel);
        contract = rsel.contract_name;
        if (verb)
            sf_warning("mdio_sel=%s contract=%s",
                       sel_text.c_str(), contract.c_str());
    } else if (verb) {
        sf_warning("Could not encode window selection: %s", serr.c_str());
    }

    mdio_pipe_stamp(out, std::string(path), fingerprint, datavar,
                    dtype, contract, MDIO_PIPE_CONTEXT_VERSION, sel_text);
    sf_fileflush(out, NULL);
}

/* Write optional hfile/bfile from dataset metadata. */
static void write_segy_side_files(mdio::Dataset& ds, bool verb)
{
    char* hname = sf_getstring("hfile");
    /* output SEG-Y EBCDIC text header file */
    if (NULL != hname) {
        char ahead[SF_EBCBYTES];
        if (mdio_get_text_header(ds, ahead)) {
            FILE* fp = fopen(hname, "w");
            if (NULL == fp) sf_error("Cannot open hfile \"%s\"", hname);
            fwrite(ahead, 1, SF_EBCBYTES, fp);
            fclose(fp);
            if (verb) sf_warning("Text header written to \"%s\"", hname);
        } else if (verb) {
            sf_warning("No text header in dataset; hfile not written");
        }
    }

    char* bname = sf_getstring("bfile");
    /* output SEG-Y binary header file */
    if (NULL != bname) {
        char bhead[SF_BNYBYTES];
        if (mdio_get_binary_header(ds, bhead)) {
            FILE* fp = fopen(bname, "wb");
            if (NULL == fp) sf_error("Cannot open bfile \"%s\"", bname);
            fwrite(bhead, 1, SF_BNYBYTES, fp);
            fclose(fp);
            if (verb) sf_warning("Binary header written to \"%s\"", bname);
        } else if (verb) {
            sf_warning("No binary header in dataset; bfile not written");
        }
    }
}

/* Chunk-aligned streaming read over window (contiguous RSF blocks). */
static void stream_amplitude(mdio::Dataset& ds, const std::string& datavar,
                             const char* path,
                             const std::vector<MdioAxis>& axes,
                             const ReadWindow& win, long ntr,
                             sf_file out, bool verb)
{
    const int rank = (int) axes.size();

    int panel_cap;
    if (!sf_getint("panel", &panel_cap)) panel_cap = 64;
    /* max traces per streamed block along the fastest spatial axis; fallback
       sizing only, used when no whole-chunk block fits blockmb= */

    int blockmb;
    if (!sf_getint("blockmb", &blockmb)) blockmb = 128;
    /* block budget (MiB); whole-chunk blocks amortize decompression */
    static_assert(MDIO_PIPE_BLOCK_MB == 128,
                  "sfdoc needs a literal default; keep it equal to the macro");

    MdioBlockLoop loop;
    std::string note;
    mdio_pipe_resolve_block_loop(std::string(path), datavar, win.n,
                                 panel_cap, blockmb, loop, note);
    if (verb && !note.empty()) sf_warning("%s", note.c_str());

    MdioFloat32Var amp_var;
    if (!mdio_get_float32_var(ds, datavar, amp_var))
        sf_error("Cannot open float32 variable \"%s\"", datavar.c_str());

    std::vector<float> pbuf((size_t) loop.block_floats);
    std::vector<long> oidx(loop.outer_shape.size(), 0);

    if (verb)
        sf_warning("Streaming %ld traces in blocks of %ld along \"%s\" "
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

            std::vector<mdio::RangeDescriptor<mdio::Index> > pslices;
            pslices.reserve((size_t) rank);
            for (int a = 0; a < rank; a++) {
                long start, stop;
                if (a < loop.pivot) {
                    start = win.f[(size_t) a] + oidx[(size_t) a];
                    stop  = start + 1;
                } else if (a == loop.pivot) {
                    start = win.f[(size_t) a] + ps;
                    stop  = start + pc;
                } else {
                    start = win.f[(size_t) a];
                    stop  = start + win.n[(size_t) a];
                }
                mdio::RangeDescriptor<mdio::Index> r = {
                    axes[(size_t) a].label.c_str(),
                    (mdio::Index) start,
                    (mdio::Index) stop,
                    1};
                pslices.push_back(r);
            }

            if (!mdio_read_float_block(amp_var, pslices, pbuf.data(), block_n))
                sf_error("Failed to read data block (var \"%s\", %ld floats)",
                         datavar.c_str(), block_n);
            sf_floatwrite(pbuf.data(), block_n, out);
        }
    }
}

/* j#>1: read bbox, decimate in memory (non-same-geometry only). */
static void write_amplitude_decimated(mdio::Dataset& sliced,
                                      const std::string& datavar,
                                      const ReadWindow& win, long total,
                                      sf_file out)
{
    const int rank = (int) win.n.size();
    std::vector<double> draw;
    if (!mdio_read_field(sliced, datavar, draw))
        sf_error("Failed to read data variable \"%s\"", datavar.c_str());

    float* obuf = sf_floatalloc(total);
    std::vector<long> idx((size_t) rank);
    for (long t = 0; t < total; t++) {
        mdio_pipe_decode_index(t, rank, &win.n[0], &idx[0]);
        long s = 0;
        for (int a = 0; a < rank; a++)
            s = s * win.span[(size_t) a] + idx[(size_t) a] * win.j[(size_t) a];
        obuf[t] = (s < (long) draw.size()) ? (float) draw[s] : 0.0f;
    }
    sf_floatwrite(obuf, total, out);
    free(obuf);
}

/* tfile: SF_NKEYS ints per trace on windowed trace grid. */
static void write_trace_headers(
    mdio::Dataset& sliced, const char* path,
    const std::vector<mdio::RangeDescriptor<mdio::Index> >& slices,
    const std::string& headersVar, const ReadWindow& win, long ntr,
    sf_file hdr)
{
    const int tr = (int) win.n.size() - 1;   /* trace-grid axes */

    std::vector<std::vector<int> > keyvals(SF_NKEYS);
    for (int k = 0; k < SF_NKEYS; k++) {
        const char* nm = segykeyword(k);
        std::vector<double> vals;
        if (!mdio_read_header_key(sliced, std::string(path), slices,
                                  headersVar, nm, vals))
            continue;

        std::vector<int> dec((size_t) ntr, 0);
        std::vector<long> idx((size_t) (tr > 0 ? tr : 1));
        for (long t = 0; t < ntr; t++) {
            long s = 0;
            if (tr > 0) {
                mdio_pipe_decode_index(t, tr, &win.n[0], &idx[0]);
                for (int a = 0; a < tr; a++)
                    s = s * win.span[(size_t) a] +
                        idx[(size_t) a] * win.j[(size_t) a];
            }
            dec[(size_t) t] = (s < (long) vals.size()) ? (int) vals[s] : 0;
        }
        keyvals[k] = dec;
    }

    int* itrace = sf_intalloc(SF_NKEYS);
    for (long t = 0; t < ntr; t++) {
        for (int k = 0; k < SF_NKEYS; k++)
            itrace[k] = keyvals[k].empty() ? 0 : keyvals[k][(size_t) t];
        sf_intwrite(itrace, SF_NKEYS, hdr);
    }
    free(itrace);
}

int main(int argc, char* argv[])
{
    sf_init(argc, argv);

    bool verb;
    if (!sf_getbool("verb", &verb)) verb = false;
    /* Verbosity flag */

    char* path = sf_getstring("mdio");
    /* input MDIO dataset (path or gs://, s3:// URL) */
    if (NULL == path) sf_error("Need mdio=");

    char* dataname = sf_getstring("data");
    /* name of the MDIO data variable (default: auto-detect) */
    char* hdrsname = sf_getstring("headers");
    /* name of the MDIO trace-headers variable (default: auto-detect) */

    const char* read = sf_getstring("read");
    /* ---- outputs ---- */
    if (NULL == read) read = "d";

    /* ---- open dataset (lazy) ----
       Remote stores are materialized locally first: Dataset::Open on s3://
       Zarr-V3 drops the bucket from per-variable specs. Stamp still uses the
       original uri so the writer can find the real parent. */
    std::string open_path, mat_err;
    if (!mdio_pipe_materialize_local(std::string(path), open_path, mat_err))
        sf_error("Cannot materialize MDIO \"%s\": %s", path, mat_err.c_str());
    MdioStagedOpen stage_guard;
    if (open_path != std::string(path))
        stage_guard.reset(open_path);

    auto dsFut = mdio::Dataset::Open(open_path, mdio::constants::kOpen);
    if (!dsFut.status().ok())
        sf_error("Cannot open MDIO \"%s\": %s", path,
                 dsFut.status().ToString().c_str());
    mdio::Dataset ds = dsFut.value();

    std::string datavar = mdio_data_variable(ds, dataname);
    if (datavar.empty())
        sf_error("Could not find a data variable in \"%s\"", path);
    std::string headersVar = mdio_headers_variable(ds, hdrsname);

    if (verb) sf_warning("data variable: %s", datavar.c_str());

    /* float32 only — wider dtypes lose precision on SF_FLOAT hop */
    std::string dtype = mdio_variable_dtype(ds, datavar);
    if (dtype != MDIO_PIPE_DTYPE_FLOAT32)
        sf_error("Data variable \"%s\" has dtype \"%s\"; the lossless pipe "
                 "carries float32 only (the RSF stream between sfmdioread and "
                 "sfmdiowrite is SF_FLOAT)",
                 datavar.c_str(),
                 dtype.empty() ? "unknown" : dtype.c_str());

    /* fingerprint unwindowed parent (writer verifies this) */
    std::string fingerprint = mdio_pipe_fingerprint(ds, datavar, dtype);
    if (verb) sf_warning("mdio fingerprint: %s", fingerprint.c_str());

    std::vector<MdioAxis> axes = mdio_axes(ds, datavar);
    const int rank = (int) axes.size();
    if (rank < 1)
        sf_error("Data variable \"%s\" has no dimensions", datavar.c_str());

    ReadWindow win = resolve_window(axes);
    std::vector<mdio::RangeDescriptor<mdio::Index> > slices =
        window_bounding_box(axes, win);

    mdio::Dataset sliced = ds;
    {
        auto r = ds.isel(slices);
        if (!r.status().ok()) sf_error("Failed to slice MDIO dataset");
        sliced = r.value();
    }

    const int sa = rank - 1;             /* MDIO sample axis */
    const long ns = win.n[(size_t) sa];
    long ntr = 1;
    for (int a = 0; a < sa; a++) ntr *= win.n[(size_t) a];
    const long total = ns * ntr;

    if (verb) sf_warning("Reading %ld samples x %ld traces", ns, ntr);

    /* ---- outputs ---- */
    sf_file out = NULL, hdr = NULL;
    if (read[0] != 'h') out = open_amplitude_output(axes, win);
    if (read[0] != 'd') hdr = open_header_output(ntr, out);
    else                segy_init(SF_NKEYS, NULL);

    if (NULL != out)
        stamp_pipe_context(out, path, fingerprint, datavar, dtype,
                           axes, win, verb);

    write_segy_side_files(ds, verb);

    if (NULL != out) {
        bool unit_stride = true;
        for (int a = 0; a < rank; a++)
            if (win.j[(size_t) a] != 1) { unit_stride = false; break; }

        if (unit_stride)
            stream_amplitude(ds, datavar, open_path.c_str(), axes, win, ntr,
                             out, verb);
        else
            write_amplitude_decimated(sliced, datavar, win, total, out);
    }

    if (NULL != hdr)
        write_trace_headers(sliced, open_path.c_str(), slices, headersVar, win,
                            ntr, hdr);

    exit(0);
}
