/* MDIO (mdio-cpp) <-> Madagascar SEG-Y/RSF helpers for sfmdioread/sfmdiowrite. */

#ifndef _mdio2segy_hh
#define _mdio2segy_hh

#include <string>
#include <vector>

#include <mdio/mdio.h>
#include <nlohmann/json.hpp>

extern "C" {
#include <rsf.h>
#include <rsfsegy.h>
}

/* One MDIO dimension (slowest first, sample last). */
struct MdioAxis {
    std::string label;
    long        size;       /* axis extent                                 */
    double      o;          /* origin                                      */
    double      d;          /* sampling                                    */
    std::string unit;
    bool        sample;     /* fastest axis -> RSF n1                      */
};

/* Variable element dtype as TensorStore name, or "". */
std::string mdio_variable_dtype(mdio::Dataset& ds, const std::string& name);

/* Principal data variable; given verbatim else largest sample-typed var. */
std::string mdio_data_variable(mdio::Dataset& ds, const char* given);

/* Headers variable; given else headers/trace_headers or structured var. */
std::string mdio_headers_variable(mdio::Dataset& ds, const char* given);

/* Read SEG-Y key over sliced trace grid (per-field vars or structured headers). */
bool mdio_read_header_key(mdio::Dataset& ds, const std::string& path,
                          const std::vector<mdio::RangeDescriptor<mdio::Index> >& slices,
                          const std::string& headersVar, const char* key,
                          std::vector<double>& out);

/* Data-variable axis geometry (MDIO order). */
std::vector<MdioAxis> mdio_axes(mdio::Dataset& ds, const std::string& datavar);

/* Read whole variable as doubles (common MDIO dtypes). */
bool mdio_read_field(mdio::Dataset& ds, const std::string& name,
                     std::vector<double>& out);

typedef mdio::Variable<mdio::dtypes::float32_t> MdioFloat32Var;

/* Resolve float32 data var once (above streaming loop). */
bool mdio_get_float32_var(mdio::Dataset& ds, const std::string& datavar,
                          MdioFloat32Var& out);

/* Read/write float32 block; optional slice (absolute coords). */
bool mdio_read_float_block(
    MdioFloat32Var& var,
    const std::vector<mdio::RangeDescriptor<mdio::Index> >& slices,
    float* out, long count);

bool mdio_write_float_block(
    MdioFloat32Var& var,
    const std::vector<mdio::RangeDescriptor<mdio::Index> >& slices,
    const float* buf, long count);

/* o/d from coordinate variable. */
bool mdio_coord(mdio::Dataset& ds, const std::string& label,
                double* o, double* d);

/* SEG-Y text/binary from segy_file_header attributes. */
bool mdio_get_text_header(mdio::Dataset& ds, char ahead[SF_EBCBYTES]);
bool mdio_get_binary_header(mdio::Dataset& ds, char bhead[SF_BNYBYTES]);

/* Attach SEG-Y headers to new-dataset schema (segy_file_header attrs). */
void mdio_put_text_header(nlohmann::json& schema,
                          const char ahead[SF_EBCBYTES]);
void mdio_put_binary_header(nlohmann::json& schema,
                            const char bhead[SF_BNYBYTES]);

/* MDIO var name for SEG-Y key (with aliases). */
std::string mdio_header_field(mdio::Dataset& ds, const char* key);

/* Build MDIO v1 schema for new dataset. */
nlohmann::json mdio_build_schema(const std::string& name,
                                 const std::vector<MdioAxis>& axes,
                                 const std::string& datavar,
                                 const std::vector<std::string>& keys,
                                 const std::vector<long>& chunk);

#endif
