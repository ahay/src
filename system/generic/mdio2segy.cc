/* MDIO <-> SEG-Y/RSF helpers. See mdio2segy.hh. */

#include <algorithm>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

#include "mdio2segy.hh"

/* SEG-Y key -> MDIO field aliases */
static std::vector<std::string> header_aliases(const std::string& key)
{
    std::vector<std::string> a;
    a.push_back(key);
    if (key == "iline")  { a.push_back("inline"); a.push_back("inline_number");
                           a.push_back("il"); }
    else if (key == "xline") { a.push_back("crossline");
                               a.push_back("crossline_number"); a.push_back("xl"); }
    else if (key == "cdp") { a.push_back("cdp_number"); a.push_back("ensemble"); }
    else if (key == "sx")  { a.push_back("source_coord_x"); a.push_back("source_x");
                            a.push_back("sou_x"); }
    else if (key == "sy")  { a.push_back("source_coord_y"); a.push_back("source_y");
                            a.push_back("sou_y"); }
    else if (key == "gx")  { a.push_back("group_coord_x"); a.push_back("rec_x"); }
    else if (key == "gy")  { a.push_back("group_coord_y"); a.push_back("rec_y"); }
    else if (key == "cdpx"){ a.push_back("cdp_x"); }
    else if (key == "cdpy"){ a.push_back("cdp_y"); }
    return a;
}

/* Read variable as doubles (try common dtypes) */
bool mdio_read_field(mdio::Dataset& ds, const std::string& name,
                     std::vector<double>& out)
{
#define MDIO_TRY_READ(MT)                                                      \
    do {                                                                       \
        auto vr = ds.variables.get<MT>(name);                                  \
        if (vr.status().ok()) {                                                 \
            auto var = vr.value();                                             \
            auto fut = var.Read();                                             \
            if (fut.status().ok()) {                                           \
                auto vd  = fut.value();                                        \
                long n   = (long) vd.num_samples();                           \
                auto off = vd.get_flattened_offset();                          \
                const MT* p =                                                  \
                    static_cast<const MT*>(vd.get_data_accessor().data());     \
                out.resize((size_t) n);                                        \
                for (long i = 0; i < n; i++)                                   \
                    out[(size_t) i] = (double) p[off + i];                     \
                return true;                                                   \
            }                                                                  \
        }                                                                      \
    } while (0)

    MDIO_TRY_READ(mdio::dtypes::int32_t);
    MDIO_TRY_READ(mdio::dtypes::int16_t);
    MDIO_TRY_READ(mdio::dtypes::int8_t);
    MDIO_TRY_READ(mdio::dtypes::int64_t);
    MDIO_TRY_READ(mdio::dtypes::uint32_t);
    MDIO_TRY_READ(mdio::dtypes::uint16_t);
    MDIO_TRY_READ(mdio::dtypes::uint8_t);
    MDIO_TRY_READ(mdio::dtypes::bool_t);
    MDIO_TRY_READ(mdio::dtypes::float32_t);
    MDIO_TRY_READ(mdio::dtypes::float64_t);

#undef MDIO_TRY_READ
    return false;
}

bool mdio_coord(mdio::Dataset& ds, const std::string& label,
                double* o, double* d)
{
    std::vector<double> v;
    if (!mdio_read_field(ds, label, v) || v.empty()) return false;
    *o = v[0];
    *d = (v.size() > 1) ? (v[1] - v[0]) : 1.0;
    return true;
}

bool mdio_get_float32_var(mdio::Dataset& ds, const std::string& datavar,
                          MdioFloat32Var& out)
{
    auto vr = ds.variables.get<mdio::dtypes::float32_t>(datavar);
    if (!vr.status().ok()) return false;
    out = vr.value();
    return true;
}

/* MdioFloat32Var or slice() Variable<> */
template <typename Var>
static bool read_into(Var& var, float* out, long count)
{
    auto fut = var.Read();
    if (!fut.status().ok()) return false;
    auto vd = fut.value();
    auto off = vd.get_flattened_offset();
    const float* p = static_cast<const float*>(vd.get_data_accessor().data());
    const long n = (long) vd.num_samples();
    const long ncopy = (n < count) ? n : count;
    for (long i = 0; i < ncopy; i++) out[i] = p[off + i];
    for (long i = ncopy; i < count; i++) out[i] = 0.0f;
    return true;
}

template <typename Var>
static bool write_from(Var& var, const float* buf, long count)
{
    auto vdr = mdio::from_variable<mdio::dtypes::float32_t>(var);
    if (!vdr.status().ok()) return false;
    auto vd = vdr.value();
    auto off = vd.get_flattened_offset();
    float* p = static_cast<float*>(vd.get_data_accessor().data());
    const long n = (long) vd.num_samples();
    const long ncopy = (n < count) ? n : count;
    for (long i = 0; i < ncopy; i++) p[off + i] = buf[i];
    return var.Write(vd).status().ok();
}

bool mdio_read_float_block(
    MdioFloat32Var& var,
    const std::vector<mdio::RangeDescriptor<mdio::Index> >& slices,
    float* out, long count)
{
    if (!out || count < 0) return false;
    if (slices.empty()) return read_into(var, out, count);

    auto sr = var.slice(slices);
    if (!sr.status().ok()) return false;
    auto svar = sr.value();
    return read_into(svar, out, count);
}

bool mdio_write_float_block(
    MdioFloat32Var& var,
    const std::vector<mdio::RangeDescriptor<mdio::Index> >& slices,
    const float* buf, long count)
{
    if (!buf || count < 0) return false;
    if (slices.empty()) return write_from(var, buf, count);

    /* slice uses absolute coords; view shares kvstore region */
    auto sr = var.slice(slices);
    if (!sr.status().ok()) return false;
    auto svar = sr.value();
    return write_from(svar, buf, count);
}

/* Variable resolution */
static bool has_empty_label(mdio::Dataset& ds, const std::string& name)
{
    auto vr = ds.variables.at(name);
    if (!vr.status().ok()) return false;
    auto dom = vr.value().dimensions();
    int rank = (int) dom.rank();
    auto labels = dom.labels();
    for (int i = 0; i < rank; i++)
        if (std::string(labels[i]).empty()) return true;
    return false;
}

/* 1-D coord var (name == dim label); not per-trace headers */
static bool is_dim_coordinate(mdio::Dataset& ds, const std::string& name)
{
    auto vr = ds.variables.at(name);
    if (!vr.status().ok()) return false;
    auto dom = vr.value().dimensions();
    if ((int) dom.rank() != 1) return false;
    auto labels = dom.labels();
    return std::string(labels[0]) == name;
}

/* Sample dtype candidate; excludes bool/byte/string kinds */
static bool is_sample_dtype(const std::string& dt)
{
    static const char* kSample[] = {
        "int8", "int16", "int32", "int64",
        "uint8", "uint16", "uint32", "uint64",
        "float16", "bfloat16", "float32", "float64",
        "complex64", "complex128", NULL};
    for (int i = 0; kSample[i]; i++)
        if (dt == kSample[i]) return true;
    return false;
}

std::string mdio_variable_dtype(mdio::Dataset& ds, const std::string& name)
{
    auto vr = ds.variables.at(name);
    if (!vr.status().ok()) return "";
    return std::string(vr.value().dtype().name());
}

std::string mdio_data_variable(mdio::Dataset& ds, const char* given)
{
    if (given && *given) return std::string(given);

    auto names = ds.variables.get_iterable_accessor();

    /* prefer seismic/amplitude/data names */
    const char* preferred[] = {"seismic", "amplitude", "data", NULL};
    for (int k = 0; preferred[k]; k++) {
        std::string p = preferred[k];
        if (ds.variables.contains_key(p) && !has_empty_label(ds, p) &&
            is_sample_dtype(mdio_variable_dtype(ds, p)))
            return p;
    }

    /* else largest-rank sample-typed non-structured var */
    std::string best;
    int bestrank = 1;
    for (size_t i = 0; i < names.size(); i++) {
        const std::string& n = names[i];
        auto vr = ds.variables.at(n);
        if (!vr.status().ok()) continue;
        int rank = (int) vr.value().dimensions().rank();
        if (rank < 2 || has_empty_label(ds, n)) continue;
        if (!is_sample_dtype(mdio_variable_dtype(ds, n))) continue;
        if (rank > bestrank) { best = n; bestrank = rank; }
    }
    return best;
}

std::string mdio_headers_variable(mdio::Dataset& ds, const char* given)
{
    if (given && *given) return std::string(given);

    const char* common[] = {"headers", "trace_headers", "TraceHeaders", NULL};
    for (int k = 0; common[k]; k++)
        if (ds.variables.contains_key(common[k])) return std::string(common[k]);

    /* structured var has unlabeled byte dimension */
    auto names = ds.variables.get_iterable_accessor();
    for (size_t i = 0; i < names.size(); i++)
        if (has_empty_label(ds, names[i])) return names[i];

    return "";
}

static std::string header_field_name(mdio::Dataset& ds, const char* key)
{
    std::vector<std::string> cands = header_aliases(key);
    for (size_t i = 0; i < cands.size(); i++)
        if (ds.variables.contains_key(cands[i]) &&
            !is_dim_coordinate(ds, cands[i]))
            return cands[i];
    return "";
}

bool mdio_read_header_key(mdio::Dataset& ds, const std::string& path,
                          const std::vector<mdio::RangeDescriptor<mdio::Index> >& slices,
                          const std::string& headersVar, const char* key,
                          std::vector<double>& out)
{
    /* layout 1: one var per SEG-Y field */
    
    std::string f = header_field_name(ds, key);
    if (!f.empty()) return mdio_read_field(ds, f, out);

    /* layout 2: structured headers + SelectField (reopen per field) */
    if (headersVar.empty()) return false;

    std::vector<std::string> cands = header_aliases(key);
    for (size_t i = 0; i < cands.size(); i++) {
        auto f2 = mdio::Dataset::Open(path, mdio::constants::kOpen);
        if (!f2.status().ok()) continue;
        mdio::Dataset d2 = f2.value();

        mdio::Dataset sl = d2;
        if (!slices.empty()) {
            auto r = d2.isel(slices);
            if (!r.status().ok()) continue;
            sl = r.value();
        }

        auto sf = sl.SelectField(headersVar, cands[i]);
        if (!sf.status().ok()) continue;

        if (mdio_read_field(sl, headersVar, out)) return true;
    }
    return false;
}

/* Axis geometry */
std::vector<MdioAxis> mdio_axes(mdio::Dataset& ds, const std::string& datavar)
{
    std::vector<MdioAxis> axes;
    auto vr = ds.variables.at(datavar);
    if (!vr.status().ok()) return axes;

    auto dom = vr.value().dimensions();
    int rank = (int) dom.rank();
    auto labels = dom.labels();
    auto shape  = dom.shape();

    for (int i = 0; i < rank; i++) {
        MdioAxis a;
        a.label  = std::string(labels[i]);
        a.size   = (long) shape[i];
        a.o      = 0.0;
        a.d      = 1.0;
        a.sample = (i == rank - 1);
        double o, d;
        if (mdio_coord(ds, a.label, &o, &d)) { a.o = o; a.d = d; }
        axes.push_back(a);
    }
    return axes;
}

/* SEG-Y file headers in segy_file_header attrs (mdio-cpp or Python layout) */
static const char* SEGY_FILE_HEADER_VAR = "segy_file_header";

static bool segy_file_header_meta(mdio::Dataset& ds, nlohmann::json& out)
{
    auto vr = ds.variables.at(SEGY_FILE_HEADER_VAR);
    if (!vr.status().ok()) return false;
    out = vr.value().getMetadata();
    return true;
}

static const nlohmann::json* find_node(const nlohmann::json& meta,
                                        const char* name)
{
    /* mdio-cpp: metadata.attributes.name */
    if (meta.contains("metadata") && meta["metadata"].is_object()) {
        const nlohmann::json& m = meta["metadata"];
        if (m.contains("attributes") && m["attributes"].contains(name))
            return &m["attributes"][name];
        if (m.contains(name)) return &m[name];
    }
    /* Python mdio: attributes or top-level */
    if (meta.contains("attributes") && meta["attributes"].contains(name))
        return &meta["attributes"][name];
    if (meta.contains(name)) return &meta[name];
    return NULL;
}

bool mdio_get_text_header(mdio::Dataset& ds, char ahead[SF_EBCBYTES])
{
    nlohmann::json meta;
    if (!segy_file_header_meta(ds, meta)) return false;
    const nlohmann::json* node = find_node(meta, "textHeader");
    if (NULL == node) return false;

    memset(ahead, ' ', SF_EBCBYTES);
    if (node->is_string()) {
        /* string: strip newlines from joined 40 cards */
        std::string raw = node->get<std::string>();
        std::string s;
        s.reserve(raw.size());
        for (size_t i = 0; i < raw.size(); i++)
            if (raw[i] != '\n' && raw[i] != '\r') s.push_back(raw[i]);
        memcpy(ahead, s.data(), std::min((size_t) SF_EBCBYTES, s.size()));
    } else if (node->is_array()) {
        /* array of 80-char cards */
        int off = 0;
        for (size_t i = 0; i < node->size() && off < SF_EBCBYTES; i++) {
            std::string s = (*node)[i].get<std::string>();
            int n = (int) std::min((size_t)(SF_EBCBYTES - off),
                                   std::min((size_t) 80, s.size()));
            memcpy(ahead + off, s.data(), n);
            off += 80;
        }
    } else {
        return false;
    }
    return true;
}

bool mdio_get_binary_header(mdio::Dataset& ds, char bhead[SF_BNYBYTES])
{
    nlohmann::json meta;
    if (!segy_file_header_meta(ds, meta)) return false;
    const nlohmann::json* node = find_node(meta, "binaryHeader");
    if (NULL == node) return false;

    memset(bhead, 0, SF_BNYBYTES);
    if (node->is_array()) {
        for (size_t i = 0; i < node->size() && i < (size_t) SF_BNYBYTES; i++)
            bhead[i] = (char) (int) (*node)[i].get<int>();
    } else if (node->is_string()) {
        std::string s = node->get<std::string>();
        memcpy(bhead, s.data(), std::min((size_t) SF_BNYBYTES, s.size()));
    } else {
        return false;
    }
    return true;
}

/* segy_file_header attrs in schema (create scalar var on first use) */
static nlohmann::json& segy_file_header_attributes(nlohmann::json& schema)
{
    if (!schema.contains("variables") || !schema["variables"].is_array())
        schema["variables"] = nlohmann::json::array();

    nlohmann::json& vars = schema["variables"];
    for (size_t i = 0; i < vars.size(); i++)
        if (vars[i].contains("name") && vars[i]["name"] == SEGY_FILE_HEADER_VAR)
            return vars[i]["metadata"]["attributes"];

    nlohmann::json v;
    v["name"]     = SEGY_FILE_HEADER_VAR;
    v["dataType"] = "int32";
    nlohmann::json dim;
    dim["name"] = SEGY_FILE_HEADER_VAR;
    dim["size"] = 1;
    v["dimensions"] = nlohmann::json::array({dim});
    v["metadata"]["attributes"] = nlohmann::json::object();
    vars.push_back(v);
    return vars[vars.size() - 1]["metadata"]["attributes"];
}

void mdio_put_text_header(nlohmann::json& schema,
                          const char ahead[SF_EBCBYTES])
{
    nlohmann::json& attrs = segy_file_header_attributes(schema);
    attrs["textHeader"] = std::string(ahead, (size_t) SF_EBCBYTES);
}

void mdio_put_binary_header(nlohmann::json& schema,
                            const char bhead[SF_BNYBYTES])
{
    nlohmann::json arr = nlohmann::json::array();
    for (int i = 0; i < SF_BNYBYTES; i++)
        arr.push_back((int) (unsigned char) bhead[i]);
    nlohmann::json& attrs = segy_file_header_attributes(schema);
    attrs["binaryHeader"] = arr;
}

/* Schema construction for reduced=y path */
static std::string iso_now(void)
{
    char buf[32];
    time_t t = time(NULL);
    struct tm* g = gmtime(&t);
    strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S.000000Z", g);
    return std::string(buf);
}

nlohmann::json mdio_build_schema(const std::string& name,
                                 const std::vector<MdioAxis>& axes,
                                 const std::string& datavar,
                                 const std::vector<std::string>& keys,
                                 const std::vector<long>& chunk)
{
    nlohmann::json schema;
    schema["metadata"]["apiVersion"] = "1.0.0";
    schema["metadata"]["name"]       = name;
    schema["metadata"]["createdOn"]  = iso_now();

    nlohmann::json vars = nlohmann::json::array();

    /* uniform coord vars o+i*d */
    for (size_t i = 0; i < axes.size(); i++) {
        nlohmann::json v;
        v["name"]     = axes[i].label;
        v["dataType"] = "float32";
        nlohmann::json dim;
        dim["name"] = axes[i].label;
        dim["size"] = axes[i].size;
        v["dimensions"] = nlohmann::json::array({dim});
        /* omit unitsV1 — RSF unit strings fail MDIO v1 validation */
        vars.push_back(v);
    }

    /* non-sample axis labels */
    std::vector<std::string> tracelabels;
    std::vector<size_t>      traceaxis;   /* parent axis index per label */
    for (size_t i = 0; i < axes.size(); i++)
        if (!axes[i].sample) {
            tracelabels.push_back(axes[i].label);
            traceaxis.push_back(i);
        }

    /* header chunks: full on fast trace axis, data chunks on slow axes */
    std::vector<long> tracechunk(tracelabels.size());
    for (size_t k = 0; k < tracelabels.size(); k++) {
        size_t a = traceaxis[k];
        long size = axes[a].size;
        long c;
        if (k + 1 == tracelabels.size())          /* fast trace axis: full chunk */
            c = size;
        else
            c = (a < chunk.size() && chunk[a] > 0 && chunk[a] <= size)
                ? chunk[a] : size;
        tracechunk[k] = c;
    }

    /* Data variable. */
    {
        nlohmann::json v;
        v["name"]     = datavar;
        v["dataType"] = "float32";
        nlohmann::json dims = nlohmann::json::array();
        for (size_t i = 0; i < axes.size(); i++) dims.push_back(axes[i].label);
        v["dimensions"] = dims;
        if (!chunk.empty()) {
            nlohmann::json cs = nlohmann::json::array();
            for (size_t i = 0; i < chunk.size(); i++) cs.push_back(chunk[i]);
            v["metadata"]["chunkGrid"]["name"] = "regular";
            v["metadata"]["chunkGrid"]["configuration"]["chunkShape"] = cs;
        }
        vars.push_back(v);
    }

    /* One int32 trace-header variable per SEG-Y key. */
    for (size_t k = 0; k < keys.size(); k++) {
        nlohmann::json v;
        v["name"]     = keys[k];
        v["dataType"] = "int32";
        nlohmann::json dims = nlohmann::json::array();
        for (size_t i = 0; i < tracelabels.size(); i++)
            dims.push_back(tracelabels[i]);
        v["dimensions"] = dims;
        /* chunk slow trace axes only (2D default) */
        if (tracelabels.size() >= 2) {
            nlohmann::json cs = nlohmann::json::array();
            for (size_t i = 0; i < tracechunk.size(); i++) cs.push_back(tracechunk[i]);
            v["metadata"]["chunkGrid"]["name"] = "regular";
            v["metadata"]["chunkGrid"]["configuration"]["chunkShape"] = cs;
        }
        vars.push_back(v);
    }

    schema["variables"] = vars;
    return schema;
}
