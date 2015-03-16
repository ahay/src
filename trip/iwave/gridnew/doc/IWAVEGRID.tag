<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="page">
    <name>index</name>
    <title>IWAVE Grid I/O package</title>
    <filename>index</filename>
  </compound>
  <compound kind="file">
    <name>create_hfile.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>create__hfile_8hh</filename>
    <member kind="function">
      <type>void</type>
      <name>create_hfile</name>
      <anchorfile>create__hfile_8hh.html</anchorfile>
      <anchor>a542618f6a1315204e960f8216427a600</anchor>
      <arglist>(string hfile, int prec=0)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>dot</name>
      <anchorfile>create__hfile_8hh.html</anchorfile>
      <anchor>a9e4ce5fdb1d31a716136cbcfcf96b08d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>doc.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>doc_8h</filename>
  </compound>
  <compound kind="file">
    <name>exchangeinfo.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>exchangeinfo_8h</filename>
    <class kind="struct">s_EXCHANGEINFO</class>
    <member kind="typedef">
      <type>struct s_EXCHANGEINFO</type>
      <name>EXCHANGEINFO</name>
      <anchorfile>exchangeinfo_8h.html</anchorfile>
      <anchor>a42f3e18761e322c6c9bc3d08dee7abab</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ei_setnull</name>
      <anchorfile>exchangeinfo_8h.html</anchorfile>
      <anchor>a268d94641c447c4c66c5ad6e1467c69a</anchor>
      <arglist>(EXCHANGEINFO *einfo)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ei_destroy</name>
      <anchorfile>exchangeinfo_8h.html</anchorfile>
      <anchor>add6ae25b485976ec3e803c655e413506</anchor>
      <arglist>(EXCHANGEINFO *einfo)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>extop.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>extop_8hh</filename>
    <includes id="gridpp_8hh" name="gridpp.hh" local="yes" imported="no">gridpp.hh</includes>
    <class kind="class">TSOpt::ExtOp</class>
    <namespace>TSOpt</namespace>
  </compound>
  <compound kind="file">
    <name>grid.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>grid_8h</filename>
    <class kind="struct">s_axis</class>
    <class kind="struct">grid</class>
    <member kind="define">
      <type>#define</type>
      <name>TOL</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a156b862ebf6d213f5da19b9e3ccb779e</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct s_axis</type>
      <name>axis</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>ae045cbf1e1e13b264909a1fb9ac75e86</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_default_axis</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a92ee2ca8b58522969904c83bf37fbd36</anchor>
      <arglist>(axis *a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_axis</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a27ea20ab0ab6aa1bc2c602b6943a695e</anchor>
      <arglist>(axis *a, size_t n, ireal d, ireal o)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>copy_axis</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a98d94f27eee3a468cefbc35d2492f028</anchor>
      <arglist>(axis *tgt, const axis *src)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_axis</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a688fe60460e05fd246edf01b40c12bc3</anchor>
      <arglist>(axis a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>fprint_axis</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>aafd262df63fb3179387b981672059ad5</anchor>
      <arglist>(FILE *fp, axis a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>fprint_num_axis</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>adcf4bce07588d41e91fdc411785dcae8</anchor>
      <arglist>(FILE *fp, int i, axis a)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>compare_axis</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>ab819660086e2a7c66e662366a860f278</anchor>
      <arglist>(axis a1, axis a2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_default_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a01f8027600522cc9fc783763b2719583</anchor>
      <arglist>(grid *g)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>init_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a05e42af1f6912f163252bec4f5b24820</anchor>
      <arglist>(grid *g, int dim, int gdim)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>copy_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a8d8bc71a6a38757a4109605375eb7878</anchor>
      <arglist>(grid *tgt, const grid *src)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>ad919cc58e1ed9758643d4fa4bc4d0ce6</anchor>
      <arglist>(grid a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>fprint_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a98964590b8470f2843e4706b9c68c159</anchor>
      <arglist>(FILE *fp, grid a)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>compare_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a30040740204e36f554b5cb9e4194915c</anchor>
      <arglist>(const grid g1, const grid g2)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>compatible_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a4649e93ec621db8205bfc8988ad90cef</anchor>
      <arglist>(const grid g1, const grid g2)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_datasize_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a330ac449e15cb27ab637a5d3b087c326</anchor>
      <arglist>(grid g)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>get_global_datasize_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a800020d8b792f2741e09f7ff5a6cb7c6</anchor>
      <arglist>(grid g)</arglist>
    </member>
    <member kind="function">
      <type>ireal</type>
      <name>get_cellvol_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a5b2d131c63cfe4ce3bfdb89bddaf9f18</anchor>
      <arglist>(grid g)</arglist>
    </member>
    <member kind="function">
      <type>ireal</type>
      <name>get_global_cellvol_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>ad222e0634cf85a262c9e7492656c1e5d</anchor>
      <arglist>(grid g)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_panelnum_grid</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>ab39638513bc650980c4c17ff10c09a7d</anchor>
      <arglist>(grid g)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_n</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a43cfda5c13f47f2d1552b7e9d6e0180a</anchor>
      <arglist>(IPNT n, grid g)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_d</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a93a1d18dc7b58df6028354be1c8aa397</anchor>
      <arglist>(_RPNT d, grid g)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_o</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a538d09004bcf4081119b16ab46249dc6</anchor>
      <arglist>(_RPNT o, grid g)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_gs</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a7bc59e3ceb749aa09d1f440ccd9527ef</anchor>
      <arglist>(IPNT gs, grid g)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_ge</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a7e96076542ed6f706be64a3a88b30f66</anchor>
      <arglist>(IPNT ge, grid g)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_ord</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a7757e5d5c23620bbd464d67db033e1ae</anchor>
      <arglist>(IPNT a, grid g)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>grid_union</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>aee000af1b576b6daa074731a97fae271</anchor>
      <arglist>(grid *g, axis const *ax)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>init_step</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a582cd8902115c49ca699161f8f2121dc</anchor>
      <arglist>(grid g, IPNT step, bool fwd)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>next_step</name>
      <anchorfile>grid_8h.html</anchorfile>
      <anchor>a5d8b3905f85f5462499950b7d8c2e676</anchor>
      <arglist>(grid g, IPNT step)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>gridfun.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>gridfun_8hh</filename>
    <includes id="rarray_8h" name="rarray.h" local="yes" imported="no">rarray.h</includes>
    <member kind="function">
      <type>size_t</type>
      <name>getDataSize&lt; RARR &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a272f4d8e7c24e215bc37a04eaca79e59</anchor>
      <arglist>(RARR const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>ireal *</type>
      <name>newData&lt; ireal, RARR &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aa5bcb7547a6e657dccf6c97dea6ec654</anchor>
      <arglist>(RARR &amp;md)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deleteData&lt; ireal, RARR &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aae4078612897269f99db33da7e470f61</anchor>
      <arglist>(ireal **d, RARR **md)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>writeMeta&lt; RARR &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>adb891322716b71a5d50fed1091334dc2</anchor>
      <arglist>(RARR const &amp;md, ostream &amp;s)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>gridio.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>gridio_8h</filename>
    <includes id="grid_8h" name="grid.h" local="yes" imported="no">grid.h</includes>
    <includes id="offsets_8h" name="offsets.h" local="yes" imported="no">offsets.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>IWAVE_USE_FMGR</name>
      <anchorfile>gridio_8h.html</anchorfile>
      <anchor>aa07c024bfde132a4b2f97ea751f6ce6b</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>fseeko</name>
      <anchorfile>gridio_8h.html</anchorfile>
      <anchor>a819f057285a278044b6bc6e66dff2f2f</anchor>
      <arglist>(FILE *stream, off_t offset, int whence)</arglist>
    </member>
    <member kind="function">
      <type>off_t</type>
      <name>ftello</name>
      <anchorfile>gridio_8h.html</anchorfile>
      <anchor>a7a044f39f604e1fc5a83268ca7a4ef45</anchor>
      <arglist>(FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>read_grid</name>
      <anchorfile>gridio_8h.html</anchorfile>
      <anchor>a3175d6ca41b8b8bc6eae1e60b3e1e044</anchor>
      <arglist>(grid *g, const char *fname, FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>par_grid</name>
      <anchorfile>gridio_8h.html</anchorfile>
      <anchor>ad42727c5c50ee68759ccd165f9208a21</anchor>
      <arglist>(grid *g, PARARRAY par, FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>extend_array</name>
      <anchorfile>gridio_8h.html</anchorfile>
      <anchor>a86f3ad9e9fc55e4a856473c5828131ca</anchor>
      <arglist>(ireal *a, const IPNT rags, const IPNT ran, const IPNT gs, const IPNT n, int dim, int ax)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>adj_extend_array</name>
      <anchorfile>gridio_8h.html</anchorfile>
      <anchor>a85c457620e9a73fcb7aae2d70cfc8a64</anchor>
      <arglist>(ireal *a, const IPNT rags, const IPNT ran, const IPNT gs, const IPNT n, int dim, int ax)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rsfread</name>
      <anchorfile>gridio_8h.html</anchorfile>
      <anchor>a0dd522100ec07395eab65c7d0b961cc9</anchor>
      <arglist>(ireal *a, const IPNT gs, const IPNT n, const char *fname, float scale, FILE *fp, int panelindex)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rsfwrite_proto</name>
      <anchorfile>gridio_8h.html</anchorfile>
      <anchor>af4b3f4aae2a6479e147a642ba6686e93</anchor>
      <arglist>(ireal *a, const IPNT rags, const IPNT ran, const char *fname, const char *dname, const char *type, float scale, const char *protohdr, const char *protodata, const grid *protog, int extend, FILE *fp, int panelindex)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rsfwrite</name>
      <anchorfile>gridio_8h.html</anchorfile>
      <anchor>aebc26303b7453981c3c310cadfee492f</anchor>
      <arglist>(ireal *a, const IPNT gs, const IPNT n, const char *fname, float scale, FILE *fp, int panelindex)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>gridops.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>gridops_8hh</filename>
    <includes id="gridpp_8hh" name="gridpp.hh" local="yes" imported="no">gridpp.hh</includes>
    <class kind="class">TSOpt::GridMaskFO</class>
    <class kind="class">TSOpt::GridMaskOp</class>
    <class kind="class">TSOpt::GridWindowFO</class>
    <class kind="class">TSOpt::GridWindowOp</class>
    <class kind="class">TSOpt::GridFwdDerivFO</class>
    <class kind="class">TSOpt::GridAdjDerivFO</class>
    <class kind="class">TSOpt::GridDerivOp</class>
    <class kind="class">TSOpt::GridFwdExtendFO</class>
    <class kind="class">TSOpt::GridAdjExtendFO</class>
    <class kind="class">TSOpt::GridExtendOp</class>
    <class kind="class">TSOpt::HelmFO</class>
    <class kind="class">TSOpt::GridHelmOp</class>
    <namespace>TSOpt</namespace>
    <member kind="typedef">
      <type>GridSpace</type>
      <name>myGridSpace</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>aff5659bac3f9c3b6fb2ab70638f5e715</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>window</name>
      <anchorfile>gridops_8hh.html</anchorfile>
      <anchor>ac9c653cbb18cc5c8a7ffb6324c3b886f</anchor>
      <arglist>(Scalar a, Scalar b, Scalar w, Scalar x)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>gridpp.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>gridpp_8hh</filename>
    <includes id="grid_8h" name="grid.h" local="yes" imported="no">grid.h</includes>
    <includes id="gridio_8h" name="gridio.h" local="yes" imported="no">gridio.h</includes>
    <includes id="rarray_8h" name="rarray.h" local="yes" imported="no">rarray.h</includes>
    <class kind="class">TSOpt::GridDC</class>
    <class kind="class">TSOpt::GridDCF</class>
    <class kind="class">TSOpt::GridSpace</class>
    <namespace>TSOpt</namespace>
    <member kind="define">
      <type>#define</type>
      <name>MAX_NAMELEN</name>
      <anchorfile>gridpp_8hh.html</anchorfile>
      <anchor>aaa0a105941a262d3efc09f89db90ebf0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataSize&lt; RARR &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a272f4d8e7c24e215bc37a04eaca79e59</anchor>
      <arglist>(RARR const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>ireal *</type>
      <name>newData&lt; ireal, RARR &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aa5bcb7547a6e657dccf6c97dea6ec654</anchor>
      <arglist>(RARR &amp;md)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deleteData&lt; ireal, RARR &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aae4078612897269f99db33da7e470f61</anchor>
      <arglist>(ireal **d, RARR **md)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>writeMeta&lt; RARR &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>adb891322716b71a5d50fed1091334dc2</anchor>
      <arglist>(RARR const &amp;md, ostream &amp;s)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>mpigridio.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>mpigridio_8h</filename>
    <includes id="grid_8h" name="grid.h" local="yes" imported="no">grid.h</includes>
    <includes id="gridio_8h" name="gridio.h" local="yes" imported="no">gridio.h</includes>
    <includes id="offsets_8h" name="offsets.h" local="yes" imported="no">offsets.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>N_SLOWEST_AXIS_SLICE</name>
      <anchorfile>mpigridio_8h.html</anchorfile>
      <anchor>ac07833dff29d6a74471b73c3344fd7e4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>mpigridpp.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>mpigridpp_8hh</filename>
    <includes id="gridpp_8hh" name="gridpp.hh" local="yes" imported="no">gridpp.hh</includes>
    <class kind="class">TSOpt::MPIGridSpace</class>
    <namespace>TSOpt</namespace>
  </compound>
  <compound kind="file">
    <name>offsets.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>offsets_8h</filename>
    <member kind="function">
      <type>int</type>
      <name>get_array_offsets</name>
      <anchorfile>offsets_8h.html</anchorfile>
      <anchor>ac5c12689c3f028720cca35711d62c7d8</anchor>
      <arglist>(off_t **offs, size_t *noffs, int dim, const IPNT gs, const IPNT gn, const IPNT ls, const IPNT ln)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>rarray.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>rarray_8h</filename>
    <includes id="exchangeinfo_8h" name="exchangeinfo.h" local="yes" imported="no">exchangeinfo.h</includes>
    <class kind="struct">INFODIM</class>
    <class kind="struct">RARR</class>
    <member kind="typedef">
      <type>int(</type>
      <name>RA_CREATE_FUN</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>acccff3e4b76a45e0059985f479540e40</anchor>
      <arglist>)(RARR *arr, int ndim, IPNT v1, IPNT v2)</arglist>
    </member>
    <member kind="typedef">
      <type>int(</type>
      <name>RA_SET_FUN</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>ae1e4d1ebfe4d88b9bec86bb68dc26006</anchor>
      <arglist>)(RARR *arr, const IPNT v1, const IPNT v2)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_setnull</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>ae9358e31bc43aeabd0c40ef52c63c7e5</anchor>
      <arglist>(RARR *arr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_create_s</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>afd6d8dd6c8a0d78e7d7d1092eb3d6b1d</anchor>
      <arglist>(RARR *arr, int ndim, IPNT gs, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_create_e</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>ad6ed0fdd5aaf83ca117985f702562495</anchor>
      <arglist>(RARR *arr, int ndim, IPNT ge, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_create</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a10345de154b47607cb7e5dad852e9252</anchor>
      <arglist>(RARR *arr, int ndim, IPNT gs, IPNT ge)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_declare_s</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a2f9c54f580f162297f1a3c8c8b8c3a09</anchor>
      <arglist>(RARR *arr, int ndim, IPNT gs, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_declare_e</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a7959da6f927cd695db3894ac04d18ebf</anchor>
      <arglist>(RARR *arr, int ndim, IPNT ge, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_declare</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>afef32079339ce1f5bc8ec533557f86de</anchor>
      <arglist>(RARR *arr, int ndim, IPNT gs, IPNT ge)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_allocate</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>adacb274afb82a85de443d2b88d5d576a</anchor>
      <arglist>(RARR *arr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_destroy</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a5941a679231d2516fc84ef137e65f79f</anchor>
      <arglist>(RARR *arr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_greset_s</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a14f8d30f081684b73ecd32748c865c02</anchor>
      <arglist>(RARR *arr, const IPNT gs, const IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_greset_e</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a2b58cc8b5a719d6cbf2a35dd27546f30</anchor>
      <arglist>(RARR *arr, const IPNT ge, const IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_greset</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>af41a7801279866dde177ab7cfe43a1ac</anchor>
      <arglist>(RARR *arr, const IPNT gs, const IPNT ge)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_offset_s</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a808bbe5aef9dbf6774602eedb10d1c06</anchor>
      <arglist>(RARR *arr, const IPNT os, const IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_offset_e</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a3a67edf135d2149bca1e5fa5f5d3b579</anchor>
      <arglist>(RARR *arr, const IPNT oe, const IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_offset</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a6220433052886a91c2f51936575cb7ff</anchor>
      <arglist>(RARR *arr, const IPNT os, const IPNT oe)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_dump</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a33681cccce4356a126f1077fa5764cf2</anchor>
      <arglist>(const RARR *arr, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_print</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a01cb600917db4bce8cda018ea0f9a9be</anchor>
      <arglist>(RARR *arr, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_fprint</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a9e1774c74faf4a4a133ea7b730054481</anchor>
      <arglist>(RARR *arr, const char *path)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_write</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a9827bd411fc7f4c7941e72ccb0362e20</anchor>
      <arglist>(RARR *arr, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_fwrite</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a6f6ca334c9a1108582792fe611e6db76</anchor>
      <arglist>(RARR *arr, const char *path)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_read</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>ada0caff521cba6c2d3e919d12c33f5f2</anchor>
      <arglist>(RARR *arr, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_fread</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>aaa273469ff87f91357f42dce3c55c9ba</anchor>
      <arglist>(RARR *arr, const char *path)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_printslice</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>af8b708f152c3acb745277453ea7a1286</anchor>
      <arglist>(RARR *arr, FILE *stream, int idim, int islice)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_fprintslice</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>aa402d26e7c45fb19136b2ffe48cd069f</anchor>
      <arglist>(RARR *arr, const char *path, int idim, int islice)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_writeslice</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>afc1d1d883446d1b384056c9b4425bab8</anchor>
      <arglist>(RARR *arr, FILE *stream, int idim, int islice)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_fwriteslice</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a75722fa44a9070325d882d5cfe289bfa</anchor>
      <arglist>(RARR *arr, const char *path, int idim, int islice)</arglist>
    </member>
    <member kind="function">
      <type>ireal</type>
      <name>ra_get</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a31e3bbe47bd7d92664f58569f1f6d16f</anchor>
      <arglist>(const RARR *arr, IPNT li)</arglist>
    </member>
    <member kind="function">
      <type>ireal</type>
      <name>ra_gget</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a1871cfb002d545f29c1399604c8c4d41</anchor>
      <arglist>(const RARR *arr, IPNT gi)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ra_set</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a027872678eb26638bd7184bbd1e43549</anchor>
      <arglist>(RARR *arr, IPNT li, ireal r)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ra_gset</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a3f1399de73c94b39a008387bd324dbd7</anchor>
      <arglist>(RARR *arr, IPNT gi, ireal r)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_size</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a8f6d8840b4c0c7dfc5738c7eb9cbacc2</anchor>
      <arglist>(const RARR *arr, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_datasize</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a8fa89ade73e43e7ff41207f1dfad9999</anchor>
      <arglist>(const RARR *arr, size_t *n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_a_size</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a4af53d2ee18f06191787a630951c941f</anchor>
      <arglist>(const RARR *arr, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_a_datasize</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a991542664857c4355ca52e63eb4810f4</anchor>
      <arglist>(const RARR *arr, size_t *n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_gse</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>aeb6873c12fd215c0ea351d982003d918</anchor>
      <arglist>(const RARR *arr, IPNT gs, IPNT ge)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_se</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>af0cdf3aa7262d6f478d731e20bd0837f</anchor>
      <arglist>(const RARR *arr, IPNT s, IPNT e)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_a_gse</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>aaefe1dd771e1e397ec8e8b6a654043b3</anchor>
      <arglist>(const RARR *arr, IPNT gs, IPNT ge)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_ds</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>adba682204bab03046defd414e8cf9961</anchor>
      <arglist>(const RARR *tgt, const RARR *src, IPNT ds)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_checkbound</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a61999058d81e8019e20ba47c46aa2a15</anchor>
      <arglist>(const RARR *arr, int idim, int li)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_gcheckbound</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a322f2ae6b3de1d565fee34ce8709bf72</anchor>
      <arglist>(const RARR *arr, int idim, int gi)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_ndim</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>ab0f57b314a8d428744c1add04a2c2437</anchor>
      <arglist>(const RARR *arr, int *ndim)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_empty</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>ab8d9148f1a20bdc8426397e1b285a144</anchor>
      <arglist>(RARR *arr, int *empty)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_setempty</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a22bb1e1f8e745a1d4484317c0816fff0</anchor>
      <arglist>(RARR *arr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_setexchangeinfo</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a77ba2ead232b11055b2ee602bb89464c</anchor>
      <arglist>(RARR *arr, EXCHANGEINFO *einfo)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_overlap</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>ada6b45e232e5842b39b68dbce98e60a9</anchor>
      <arglist>(RARR *arr1, RARR *arr2, int *overlap)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_setoverlap</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a289407a8b1e61c1387eaa6d137b03422</anchor>
      <arglist>(RARR *arr1, RARR *arr2)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_zero</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a2420e9343fcf06d0204a7a26a341a804</anchor>
      <arglist>(RARR *arr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_a_zero</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a4d2e09d520f06946df12cd396d10775b</anchor>
      <arglist>(RARR *arr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_deepcopy</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>adf0fd862d7b99d472a30af9797d7ced2</anchor>
      <arglist>(RARR const *src, RARR *tgt)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_copy</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a9103eed66a80eaf8090490b648421239</anchor>
      <arglist>(RARR *arr_des, RARR *arr_src)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_a_copy</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a7ed4fb6bbc677122ca25c1bd6cab7c7d</anchor>
      <arglist>(RARR *arr_des, RARR *arr_src)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_a_scale</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a330f4899b754df9495208ca68a613e49</anchor>
      <arglist>(RARR *tgt, ireal fac)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_a_inner</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a1aed7545554ac8b71e411c7930622ca6</anchor>
      <arglist>(RARR const *arr1, RARR const *arr2, ireal *ip)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_axpy</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a9deb5eec97763638b75a9f0b0116998d</anchor>
      <arglist>(RARR *arry, RARR const *arrx, ireal a)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ra_compare_meta</name>
      <anchorfile>rarray_8h.html</anchorfile>
      <anchor>a6ea5b0f13e0797053eeca3ecd130ea18</anchor>
      <arglist>(const RARR *a, const RARR *b)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>wcreate_hfile.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/grid/include/</path>
    <filename>wcreate__hfile_8hh</filename>
    <member kind="function">
      <type>void</type>
      <name>wcreate_hfile</name>
      <anchorfile>wcreate__hfile_8hh.html</anchorfile>
      <anchor>a48a50fea61d81b4107759cd4c1e1bc86</anchor>
      <arglist>(string hfile, int prec=0)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>grid</name>
    <filename>structgrid.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>dim</name>
      <anchorfile>structgrid.html</anchorfile>
      <anchor>a8abc5fe4c70542c90707a15b7bd4b51e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>gdim</name>
      <anchorfile>structgrid.html</anchorfile>
      <anchor>adb8316308fd30e88553b8ee1e56a708c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>axis</type>
      <name>axes</name>
      <anchorfile>structgrid.html</anchorfile>
      <anchor>a9903bf52100a2e39c4f0fce78eccc4b1</anchor>
      <arglist>[RARR_MAX_NDIM]</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>INFODIM</name>
    <filename>structINFODIM.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>n</name>
      <anchorfile>structINFODIM.html</anchorfile>
      <anchor>add91debacfa3b22b128f91e7eb765322</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>gs</name>
      <anchorfile>structINFODIM.html</anchorfile>
      <anchor>a28a92929fc8f51cf42faf396989a52b7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ge</name>
      <anchorfile>structINFODIM.html</anchorfile>
      <anchor>a4607db83bcbce2e794a3b474a83012ef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>n0</name>
      <anchorfile>structINFODIM.html</anchorfile>
      <anchor>a3579267126b7b81343a5bef7ba09be6a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>gs0</name>
      <anchorfile>structINFODIM.html</anchorfile>
      <anchor>ac2a7ca83864c366dcd2770004ac032ef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ge0</name>
      <anchorfile>structINFODIM.html</anchorfile>
      <anchor>a1b426c7907da1c925a30a9ad72644245</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>RARR</name>
    <filename>structRARR.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>ndim</name>
      <anchorfile>structRARR.html</anchorfile>
      <anchor>a1a6b3a01bdbbe809c3ca8bc48e2cf311</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ireal *</type>
      <name>_s</name>
      <anchorfile>structRARR.html</anchorfile>
      <anchor>a5a04d39a8c59ab7177b9f8978aaac253</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ireal *</type>
      <name>_s0</name>
      <anchorfile>structRARR.html</anchorfile>
      <anchor>a7a0541b1980002751bd879470c26a766</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>INFODIM</type>
      <name>_dims</name>
      <anchorfile>structRARR.html</anchorfile>
      <anchor>adf8b8cd98b735b6b0864d4bb1fbcf9d4</anchor>
      <arglist>[RARR_MAX_NDIM]</arglist>
    </member>
    <member kind="variable">
      <type>ireal *restrict</type>
      <name>_s1</name>
      <anchorfile>structRARR.html</anchorfile>
      <anchor>adc8f401df62fdbd158962c8b1b1fb2b7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>s_axis</name>
    <filename>structs__axis.html</filename>
    <member kind="variable">
      <type>size_t</type>
      <name>n</name>
      <anchorfile>structs__axis.html</anchorfile>
      <anchor>aa463154814a60c6086a9f50f469cd2bd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ireal</type>
      <name>d</name>
      <anchorfile>structs__axis.html</anchorfile>
      <anchor>a2b93b145f4bf786f4d5d3f9ce8cd05f2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ireal</type>
      <name>o</name>
      <anchorfile>structs__axis.html</anchorfile>
      <anchor>a31a4517f8d49e693cc711b2aad5ef4a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>id</name>
      <anchorfile>structs__axis.html</anchorfile>
      <anchor>a2943ccf3c7be3fd89b217dcff8b1ef6a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>s_EXCHANGEINFO</name>
    <filename>structs__EXCHANGEINFO.html</filename>
    <member kind="variable">
      <type>void *</type>
      <name>buf</name>
      <anchorfile>structs__EXCHANGEINFO.html</anchorfile>
      <anchor>a2868838c146ce18c572f8bfd5592c9d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>MPI_Datatype</type>
      <name>type</name>
      <anchorfile>structs__EXCHANGEINFO.html</anchorfile>
      <anchor>a726493c8aedba6c755941f106575c934</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>MPI_Datatype</type>
      <name>type2</name>
      <anchorfile>structs__EXCHANGEINFO.html</anchorfile>
      <anchor>abf93a80b0598b567a31ec47d1bebd5c3</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>TSOpt</name>
    <filename>namespaceTSOpt.html</filename>
    <class kind="class">TSOpt::ExtOp</class>
    <class kind="class">TSOpt::GridMaskFO</class>
    <class kind="class">TSOpt::GridMaskOp</class>
    <class kind="class">TSOpt::GridWindowFO</class>
    <class kind="class">TSOpt::GridWindowOp</class>
    <class kind="class">TSOpt::GridFwdDerivFO</class>
    <class kind="class">TSOpt::GridAdjDerivFO</class>
    <class kind="class">TSOpt::GridDerivOp</class>
    <class kind="class">TSOpt::GridFwdExtendFO</class>
    <class kind="class">TSOpt::GridAdjExtendFO</class>
    <class kind="class">TSOpt::GridExtendOp</class>
    <class kind="class">TSOpt::HelmFO</class>
    <class kind="class">TSOpt::GridHelmOp</class>
    <class kind="class">TSOpt::GridDC</class>
    <class kind="class">TSOpt::GridDCF</class>
    <class kind="class">TSOpt::GridSpace</class>
    <class kind="class">TSOpt::MPIGridSpace</class>
    <member kind="typedef">
      <type>GridSpace</type>
      <name>myGridSpace</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>aff5659bac3f9c3b6fb2ab70638f5e715</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::ExtOp</name>
    <filename>classTSOpt_1_1ExtOp.html</filename>
    <templarg></templarg>
    <member kind="function">
      <type></type>
      <name>ExtOp</name>
      <anchorfile>classTSOpt_1_1ExtOp.html</anchorfile>
      <anchor>a9ca3328e86b3015d6ccfdaef2c56a8e9</anchor>
      <arglist>(ExtOp const &amp;A)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ExtOp</name>
      <anchorfile>classTSOpt_1_1ExtOp.html</anchorfile>
      <anchor>aff6a35d1dcb01411a33c28dc728b3d86</anchor>
      <arglist>(GridSpace const &amp;_dom, GridSpace const &amp;_rng, string dtype, FILE *_str)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ExtOp</name>
      <anchorfile>classTSOpt_1_1ExtOp.html</anchorfile>
      <anchor>ac8ee92ef4668b958e7482da6ac0f47e4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classTSOpt_1_1ExtOp.html</anchorfile>
      <anchor>afb203e70202860793690f7a6563b52d3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; float &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classTSOpt_1_1ExtOp.html</anchorfile>
      <anchor>aee6550d60a873958f3c216e75c99c1ac</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; float &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classTSOpt_1_1ExtOp.html</anchorfile>
      <anchor>a2f602356ccda898dfd2a711211932591</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1ExtOp.html</anchorfile>
      <anchor>aa727821c853e33dc4400cad77d109e3b</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classTSOpt_1_1ExtOp.html</anchorfile>
      <anchor>a3d3c621d3b3f0de405f9b1ee6883cbcf</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classTSOpt_1_1ExtOp.html</anchorfile>
      <anchor>a1c78f679205f735df3ad4a011d1e60bd</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridMaskFO</name>
    <filename>classTSOpt_1_1GridMaskFO.html</filename>
    <member kind="function">
      <type></type>
      <name>GridMaskFO</name>
      <anchorfile>classTSOpt_1_1GridMaskFO.html</anchorfile>
      <anchor>a33b942066138bcabb979cc0f174778dc</anchor>
      <arglist>(IPNT const &amp;_siw, IPNT const &amp;_eiw, bool _bias=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridMaskFO</name>
      <anchorfile>classTSOpt_1_1GridMaskFO.html</anchorfile>
      <anchor>ad740d14bd7df8751a7977206bf487ca2</anchor>
      <arglist>(GridMaskFO const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classTSOpt_1_1GridMaskFO.html</anchorfile>
      <anchor>a12507e3e881cec9466c825e36bbe788c</anchor>
      <arglist>(LocalDataContainer&lt; ireal &gt; &amp;x, LocalDataContainer&lt; ireal &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classTSOpt_1_1GridMaskFO.html</anchorfile>
      <anchor>aa54795db559bfed3118fd7f482c1bb86</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridMaskOp</name>
    <filename>classTSOpt_1_1GridMaskOp.html</filename>
    <member kind="function">
      <type></type>
      <name>GridMaskOp</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>ab8aca2bf1a40cedafe91747351055cc7</anchor>
      <arglist>(Space&lt; ireal &gt; const &amp;_dom, Vector&lt; ireal &gt; const &amp;_bg, RPNT const &amp;sw=RPNT_0, RPNT const &amp;ew=RPNT_0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridMaskOp</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>abf79cd139e6c834105d80072277c7325</anchor>
      <arglist>(GridMaskOp const &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GridMaskOp</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>af11fec159a1454d52659cea46d8a3414</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; ireal &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>a6398e06ef715ea716c43c219e53d15a2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; ireal &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>a8317c0b855d84097f3ba720b6f0afc61</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>a95fa014ca62476c44ce596a40ca0b3ee</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>abe2705c16641f26e2ec7ced1d4a9cd72</anchor>
      <arglist>(Vector&lt; ireal &gt; const &amp;x, Vector&lt; ireal &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>a58698b879c8ade4244012e7ff1db0869</anchor>
      <arglist>(Vector&lt; ireal &gt; const &amp;x, Vector&lt; ireal &gt; const &amp;dx, Vector&lt; ireal &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>abfac5ffcda857aea6b469e545a6580ef</anchor>
      <arglist>(Vector&lt; ireal &gt; const &amp;x, Vector&lt; ireal &gt; const &amp;dy, Vector&lt; ireal &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv2</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>ada2ff5cf3db1385bd8c030e8e9511733</anchor>
      <arglist>(const Vector&lt; ireal &gt; &amp;x, const Vector&lt; ireal &gt; &amp;dx0, const Vector&lt; ireal &gt; &amp;dx1, Vector&lt; ireal &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv2</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>a2c42be44760978ee5c1d8dba7c725557</anchor>
      <arglist>(const Vector&lt; ireal &gt; &amp;x, const Vector&lt; ireal &gt; &amp;dx0, const Vector&lt; ireal &gt; &amp;dy, Vector&lt; ireal &gt; &amp;dx1) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; ireal &gt; *</type>
      <name>clone</name>
      <anchorfile>classTSOpt_1_1GridMaskOp.html</anchorfile>
      <anchor>ae245d6bc67eccc51c4ed6d0bf6639b05</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridWindowFO</name>
    <filename>classTSOpt_1_1GridWindowFO.html</filename>
    <member kind="function">
      <type></type>
      <name>GridWindowFO</name>
      <anchorfile>classTSOpt_1_1GridWindowFO.html</anchorfile>
      <anchor>a89e7a6eeb78dc95b3b66710941061556</anchor>
      <arglist>(IPNT const &amp;_iw, bool _bias=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridWindowFO</name>
      <anchorfile>classTSOpt_1_1GridWindowFO.html</anchorfile>
      <anchor>aba48ade173f4c2ad0792fd7bdc831566</anchor>
      <arglist>(GridWindowFO const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classTSOpt_1_1GridWindowFO.html</anchorfile>
      <anchor>abd59a2f3c204f263695579bf4f69a849</anchor>
      <arglist>(LocalDataContainer&lt; ireal &gt; &amp;x, LocalDataContainer&lt; ireal &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classTSOpt_1_1GridWindowFO.html</anchorfile>
      <anchor>a8905d5aa32bb033cb5523863e2c2ec71</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridWindowOp</name>
    <filename>classTSOpt_1_1GridWindowOp.html</filename>
    <member kind="function">
      <type></type>
      <name>GridWindowOp</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>ab86e647b1bc37e8fbb2291060e157689</anchor>
      <arglist>(Space&lt; ireal &gt; const &amp;_dom, Vector&lt; ireal &gt; const &amp;_bg, RPNT const &amp;sw=RPNT_0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridWindowOp</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>a20769a045a9f2bb517e2abcb0618014c</anchor>
      <arglist>(GridWindowOp const &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GridWindowOp</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>a106212bd68b2127e97bb2adbdc019a71</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; ireal &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>af9bb69b2c80e5442ce5e3472bf6ac1ed</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; ireal &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>a950d679d60830495540f14695b9ad7cc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>a4ebd89abb0df79e58a066df7d0a4e04e</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>a6147e77ed74d60f6ef2a37a44d8a83cd</anchor>
      <arglist>(Vector&lt; ireal &gt; const &amp;x, Vector&lt; ireal &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>a13448d3bc295b47363a19d49ba3e3723</anchor>
      <arglist>(Vector&lt; ireal &gt; const &amp;x, Vector&lt; ireal &gt; const &amp;dx, Vector&lt; ireal &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>ae2a45cb7f9d6244cb5709809c376d5a2</anchor>
      <arglist>(Vector&lt; ireal &gt; const &amp;x, Vector&lt; ireal &gt; const &amp;dy, Vector&lt; ireal &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv2</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>aaa71615daf1e9510d09c5766d1ea3d26</anchor>
      <arglist>(const Vector&lt; ireal &gt; &amp;x, const Vector&lt; ireal &gt; &amp;dx0, const Vector&lt; ireal &gt; &amp;dx1, Vector&lt; ireal &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv2</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>afd4919bef0b87aacdb20301a78e2d881</anchor>
      <arglist>(const Vector&lt; ireal &gt; &amp;x, const Vector&lt; ireal &gt; &amp;dx0, const Vector&lt; ireal &gt; &amp;dy, Vector&lt; ireal &gt; &amp;dx1) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; ireal &gt; *</type>
      <name>clone</name>
      <anchorfile>classTSOpt_1_1GridWindowOp.html</anchorfile>
      <anchor>ad77029d57b7aee6aa2c4f7a6593410a4</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridFwdDerivFO</name>
    <filename>classTSOpt_1_1GridFwdDerivFO.html</filename>
    <member kind="function">
      <type></type>
      <name>GridFwdDerivFO</name>
      <anchorfile>classTSOpt_1_1GridFwdDerivFO.html</anchorfile>
      <anchor>a62e6acc3bd8fd97fd2b91e24b4823f16</anchor>
      <arglist>(int _dir, ireal _fac)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridFwdDerivFO</name>
      <anchorfile>classTSOpt_1_1GridFwdDerivFO.html</anchorfile>
      <anchor>a9716e44d8a8d0349e87a78c2b81d4c2c</anchor>
      <arglist>(GridFwdDerivFO const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classTSOpt_1_1GridFwdDerivFO.html</anchorfile>
      <anchor>a447681ef2d6d2e088ab4d65e03753bc2</anchor>
      <arglist>(LocalDataContainer&lt; ireal &gt; &amp;x, LocalDataContainer&lt; ireal &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classTSOpt_1_1GridFwdDerivFO.html</anchorfile>
      <anchor>ada127ff191b076a10c541434fec0ce9b</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridAdjDerivFO</name>
    <filename>classTSOpt_1_1GridAdjDerivFO.html</filename>
    <member kind="function">
      <type></type>
      <name>GridAdjDerivFO</name>
      <anchorfile>classTSOpt_1_1GridAdjDerivFO.html</anchorfile>
      <anchor>a632eb4153dd2282609d9dcecc6df2c82</anchor>
      <arglist>(int _dir, ireal _fac)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridAdjDerivFO</name>
      <anchorfile>classTSOpt_1_1GridAdjDerivFO.html</anchorfile>
      <anchor>a53b79af575acfa00f348d3fb4c753549</anchor>
      <arglist>(GridAdjDerivFO const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classTSOpt_1_1GridAdjDerivFO.html</anchorfile>
      <anchor>a21f910807cccec30f0339dc3149f1815</anchor>
      <arglist>(LocalDataContainer&lt; ireal &gt; &amp;x, LocalDataContainer&lt; ireal &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classTSOpt_1_1GridAdjDerivFO.html</anchorfile>
      <anchor>aae3a59f84a022a8c5991653ce0f5ea19</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridDerivOp</name>
    <filename>classTSOpt_1_1GridDerivOp.html</filename>
    <member kind="function">
      <type></type>
      <name>GridDerivOp</name>
      <anchorfile>classTSOpt_1_1GridDerivOp.html</anchorfile>
      <anchor>a8a6012ada9120092285f1d76422d4264</anchor>
      <arglist>(Space&lt; ireal &gt; const &amp;_dom, int _dir, ireal scale=REAL_ONE)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridDerivOp</name>
      <anchorfile>classTSOpt_1_1GridDerivOp.html</anchorfile>
      <anchor>ad4efec0c5c8d248234cff4cc7aea4f34</anchor>
      <arglist>(GridDerivOp const &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GridDerivOp</name>
      <anchorfile>classTSOpt_1_1GridDerivOp.html</anchorfile>
      <anchor>ab97026eb49cc7100e006b956e5bdb48f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; ireal &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classTSOpt_1_1GridDerivOp.html</anchorfile>
      <anchor>af756e704d05da35878332543f5eef1ff</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; ireal &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classTSOpt_1_1GridDerivOp.html</anchorfile>
      <anchor>a6fa31fc3dd5e15d2122343a2330249cb</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1GridDerivOp.html</anchorfile>
      <anchor>a1ed9d9aac83b82afe0b0a5541f0e5b6e</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classTSOpt_1_1GridDerivOp.html</anchorfile>
      <anchor>a75542a369accecf2a0e53907f8c84e73</anchor>
      <arglist>(Vector&lt; ireal &gt; const &amp;x, Vector&lt; ireal &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classTSOpt_1_1GridDerivOp.html</anchorfile>
      <anchor>a322a1fb7b51da1688fda21322851fbd4</anchor>
      <arglist>(Vector&lt; ireal &gt; const &amp;x, Vector&lt; ireal &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; ireal &gt; *</type>
      <name>clone</name>
      <anchorfile>classTSOpt_1_1GridDerivOp.html</anchorfile>
      <anchor>a5e60557401bb472d93bfc8bf2422ac46</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridFwdExtendFO</name>
    <filename>classTSOpt_1_1GridFwdExtendFO.html</filename>
    <member kind="function">
      <type></type>
      <name>GridFwdExtendFO</name>
      <anchorfile>classTSOpt_1_1GridFwdExtendFO.html</anchorfile>
      <anchor>a173c4bd1177470d87179b3b7981a7dd7</anchor>
      <arglist>(int _n_ext, bool _ext, ireal _fac)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridFwdExtendFO</name>
      <anchorfile>classTSOpt_1_1GridFwdExtendFO.html</anchorfile>
      <anchor>a28a3a415594fccb733bd5cd831dc0e3f</anchor>
      <arglist>(GridFwdExtendFO const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classTSOpt_1_1GridFwdExtendFO.html</anchorfile>
      <anchor>a8f6c347cebede6c86dc183bcafb85c71</anchor>
      <arglist>(LocalDataContainer&lt; ireal &gt; &amp;x, LocalDataContainer&lt; ireal &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classTSOpt_1_1GridFwdExtendFO.html</anchorfile>
      <anchor>a7d26b0b86d27c9d3fbde318f3e09b7cb</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridAdjExtendFO</name>
    <filename>classTSOpt_1_1GridAdjExtendFO.html</filename>
    <member kind="function">
      <type></type>
      <name>GridAdjExtendFO</name>
      <anchorfile>classTSOpt_1_1GridAdjExtendFO.html</anchorfile>
      <anchor>afdac699a524218cad9395c9d47fa4e91</anchor>
      <arglist>(int _n_ext, bool _ext, ireal _fac)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridAdjExtendFO</name>
      <anchorfile>classTSOpt_1_1GridAdjExtendFO.html</anchorfile>
      <anchor>ae22a3190048fca50ea9f8dd5831f23f1</anchor>
      <arglist>(GridAdjExtendFO const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classTSOpt_1_1GridAdjExtendFO.html</anchorfile>
      <anchor>ae65bb751df20a9623f52741e964f4147</anchor>
      <arglist>(LocalDataContainer&lt; ireal &gt; &amp;x, LocalDataContainer&lt; ireal &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classTSOpt_1_1GridAdjExtendFO.html</anchorfile>
      <anchor>a7d32ea94bb903c6650d22b8284162445</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridExtendOp</name>
    <filename>classTSOpt_1_1GridExtendOp.html</filename>
    <member kind="function">
      <type></type>
      <name>GridExtendOp</name>
      <anchorfile>classTSOpt_1_1GridExtendOp.html</anchorfile>
      <anchor>a54b80e69941b3f1a5111e0194132b780</anchor>
      <arglist>(Space&lt; ireal &gt; const &amp;_dom, Space&lt; ireal &gt; const &amp;_rng)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridExtendOp</name>
      <anchorfile>classTSOpt_1_1GridExtendOp.html</anchorfile>
      <anchor>ada56b25e3d6d888b99d328fb1ac48793</anchor>
      <arglist>(GridExtendOp const &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GridExtendOp</name>
      <anchorfile>classTSOpt_1_1GridExtendOp.html</anchorfile>
      <anchor>a54b341846dc7eb3dee70bc5a5c7f84c3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; ireal &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classTSOpt_1_1GridExtendOp.html</anchorfile>
      <anchor>a3d807711d26e931159226457ff9be4c0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; ireal &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classTSOpt_1_1GridExtendOp.html</anchorfile>
      <anchor>ac2ba12941fe462f94a559c1a8f900653</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1GridExtendOp.html</anchorfile>
      <anchor>a52c57d3aac16506f113d3f387a1e0447</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classTSOpt_1_1GridExtendOp.html</anchorfile>
      <anchor>a5c50d6c95f398fba9c8eae3a3b609bdf</anchor>
      <arglist>(Vector&lt; ireal &gt; const &amp;x, Vector&lt; ireal &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classTSOpt_1_1GridExtendOp.html</anchorfile>
      <anchor>a545da29b53012d0756173ad40bad3ed5</anchor>
      <arglist>(Vector&lt; ireal &gt; const &amp;x, Vector&lt; ireal &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; ireal &gt; *</type>
      <name>clone</name>
      <anchorfile>classTSOpt_1_1GridExtendOp.html</anchorfile>
      <anchor>ab72a52ccf95720753a1e1a3f7d75ea39</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::HelmFO</name>
    <filename>classTSOpt_1_1HelmFO.html</filename>
    <member kind="function">
      <type></type>
      <name>HelmFO</name>
      <anchorfile>classTSOpt_1_1HelmFO.html</anchorfile>
      <anchor>a13fabe673ae1cd8e1f045905311b0921</anchor>
      <arglist>(IPNT const &amp;_narr, RPNT const &amp;_darr, ireal _scale1=1.0f, ireal _scale2=1.0f, ireal _power=0.0f, ireal _datum=0.0f, int _DirichletSides=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>HelmFO</name>
      <anchorfile>classTSOpt_1_1HelmFO.html</anchorfile>
      <anchor>ad82215635282f2b7197758d36ad14a50</anchor>
      <arglist>(HelmFO const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classTSOpt_1_1HelmFO.html</anchorfile>
      <anchor>adf0db7ce3f1baac4373ae42c5bbdb276</anchor>
      <arglist>(LocalDataContainer&lt; ireal &gt; &amp;x, LocalDataContainer&lt; ireal &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classTSOpt_1_1HelmFO.html</anchorfile>
      <anchor>ae3a5c3431ae8453cfcf417e5919f594f</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridHelmOp</name>
    <filename>classTSOpt_1_1GridHelmOp.html</filename>
    <member kind="function">
      <type></type>
      <name>GridHelmOp</name>
      <anchorfile>classTSOpt_1_1GridHelmOp.html</anchorfile>
      <anchor>ab8b6e3be0e55e2fe1f2201a92c9d43dc</anchor>
      <arglist>(GridHelmOp const &amp;A)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridHelmOp</name>
      <anchorfile>classTSOpt_1_1GridHelmOp.html</anchorfile>
      <anchor>a00b689f39e71e05708b984d74ee0562f</anchor>
      <arglist>(Space&lt; float &gt; const &amp;_dom, RPNT _weights, float _power=0.0f, float _datum=0.0f, int _DirichletSides=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GridHelmOp</name>
      <anchorfile>classTSOpt_1_1GridHelmOp.html</anchorfile>
      <anchor>a036dd27486bb3ffe5e6dedcb5a39d865</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>LinearOp&lt; float &gt; *</type>
      <name>clone</name>
      <anchorfile>classTSOpt_1_1GridHelmOp.html</anchorfile>
      <anchor>abe90b23194ce48edbec4e0004027eb8e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; float &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classTSOpt_1_1GridHelmOp.html</anchorfile>
      <anchor>a36c1e03fcb5b3c12b4cbe8efc4d65e27</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; float &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classTSOpt_1_1GridHelmOp.html</anchorfile>
      <anchor>a5140c8fe1f1bfa8b7272ce49a9bed2f5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1GridHelmOp.html</anchorfile>
      <anchor>a6ce62a99524e323df3f4be87e8e3528e</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classTSOpt_1_1GridHelmOp.html</anchorfile>
      <anchor>a1062cf3eafc1a7a373144dd4d96644e9</anchor>
      <arglist>(const Vector&lt; float &gt; &amp;x, Vector&lt; float &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classTSOpt_1_1GridHelmOp.html</anchorfile>
      <anchor>a1236fe01dcba9e8375b198fe4eaae4fb</anchor>
      <arglist>(const Vector&lt; float &gt; &amp;x, Vector&lt; float &gt; &amp;y) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridDC</name>
    <filename>classTSOpt_1_1GridDC.html</filename>
    <member kind="function">
      <type></type>
      <name>GridDC</name>
      <anchorfile>classTSOpt_1_1GridDC.html</anchorfile>
      <anchor>a2811088d69f67ad6b7149bf0bff108c1</anchor>
      <arglist>(string const &amp;_protohdr, string const &amp;_protodata, grid const &amp;_protog, string const &amp;_data_format, string data_type, bool incore=false, ostream &amp;_outfile=cerr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GridDC</name>
      <anchorfile>classTSOpt_1_1GridDC.html</anchorfile>
      <anchor>ae01fe46567291407671e39cda7a6b78a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isInCore</name>
      <anchorfile>classTSOpt_1_1GridDC.html</anchorfile>
      <anchor>a23a7d597477aab0aa50797d2ef1ba63a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getProtohdr</name>
      <anchorfile>classTSOpt_1_1GridDC.html</anchorfile>
      <anchor>a135c8fc66f31bee35a5c9f18194fc1d1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1GridDC.html</anchorfile>
      <anchor>a20db70291795f0b51ca2d9b578e256d4</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>open_p</name>
      <anchorfile>classTSOpt_1_1GridDC.html</anchorfile>
      <anchor>a65d7694914384284774e62921a67c415</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>ContentPackage&lt; ireal, RARR &gt; &amp;</type>
      <name>get</name>
      <anchorfile>classTSOpt_1_1GridDC.html</anchorfile>
      <anchor>a40a2c0fc409042c427f89e1db8090e14</anchor>
      <arglist>(bool &amp;more)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>ContentPackage&lt; ireal, RARR &gt; const &amp;</type>
      <name>get</name>
      <anchorfile>classTSOpt_1_1GridDC.html</anchorfile>
      <anchor>a13dbe5a1ab5400540f1edaa09f8c1683</anchor>
      <arglist>(bool &amp;more) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>put</name>
      <anchorfile>classTSOpt_1_1GridDC.html</anchorfile>
      <anchor>a3b0821fa5ad5b66b6bbd94eb32bafc80</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>reset</name>
      <anchorfile>classTSOpt_1_1GridDC.html</anchorfile>
      <anchor>a5ae0790388cec7ef9656be6137909578</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridDCF</name>
    <filename>classTSOpt_1_1GridDCF.html</filename>
    <member kind="function">
      <type></type>
      <name>GridDCF</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>a46fbf53fc89dff7aba8e7e0441c6d469</anchor>
      <arglist>(string _protohdr, bool _incore=false, ostream &amp;_outfile=cerr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridDCF</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>a97dbd7175a449b57ab9f5c0ff8c5b5bb</anchor>
      <arglist>(GridDCF const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GridDCF</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>a58c05b6a5a92af5d2dc850c9e23a1e75</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>PackageContainerFactory&lt; ireal, RARR &gt; *</type>
      <name>clone</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>ab92e660321bda7e1212cca3fff27a654</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>grid const &amp;</type>
      <name>getGrid</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>a18ba8c47134a16ff64bc35e0b43f0a71</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ireal</type>
      <name>getCellVol</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>a32b4a65419ff8dcf280dbd90e86e3234</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ireal</type>
      <name>getScaleFactor</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>aea98bd03addd6a0246d2de01ac5af03f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getFilename</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>a05af9f6326fdd016e64efdd8a46d11a6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>compare</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>a84c8378ddea1612c027668dcd74f888a</anchor>
      <arglist>(DataContainerFactory const &amp;dcf) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isCompatible</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>a453ba61567aa1ad506d5a781adfb8ac5</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isIncore</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>ac86c8082edbc9beb19f6130d0243d255</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>a0960fb7c7f4a371db315790e2785fa86</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>GridDC *</type>
      <name>buildGridDC</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>a1459a0631d43487be6ed39a7075e1211</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>PackageContainer&lt; ireal, RARR &gt; *</type>
      <name>buildPC</name>
      <anchorfile>classTSOpt_1_1GridDCF.html</anchorfile>
      <anchor>a69d3e8ba45b176789a7352f5735c7582</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::GridSpace</name>
    <filename>classTSOpt_1_1GridSpace.html</filename>
    <member kind="function">
      <type></type>
      <name>GridSpace</name>
      <anchorfile>classTSOpt_1_1GridSpace.html</anchorfile>
      <anchor>a14cc44c8287237d253cc77aafacc8f14</anchor>
      <arglist>(string hdr, string dtype=&quot;notype&quot;, bool incore=false, ostream &amp;outfile=cerr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridSpace</name>
      <anchorfile>classTSOpt_1_1GridSpace.html</anchorfile>
      <anchor>a690b7537d1c6473b04535cb158067af4</anchor>
      <arglist>(grid const &amp;g, std::string tag, std::string thdr=&quot;&quot;, std::string fmt=&quot;native_ireal&quot;, bool incore=false, ostream &amp;outfile=cerr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GridSpace</name>
      <anchorfile>classTSOpt_1_1GridSpace.html</anchorfile>
      <anchor>a37bdafb3f02be13b27fffa0dc7f40b92</anchor>
      <arglist>(GridSpace const &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GridSpace</name>
      <anchorfile>classTSOpt_1_1GridSpace.html</anchorfile>
      <anchor>a654c9f6e7e27cec5b0e5c07ed159c600</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>LinearAlgebraPackage&lt; ireal &gt; const &amp;</type>
      <name>getLAP</name>
      <anchorfile>classTSOpt_1_1GridSpace.html</anchorfile>
      <anchor>a52ccefb6b4c55facf7b92dee863e5e8a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DataContainerFactory const &amp;</type>
      <name>getDCF</name>
      <anchorfile>classTSOpt_1_1GridSpace.html</anchorfile>
      <anchor>a1c29846b564a4887979dbb95183b029b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>grid const &amp;</type>
      <name>getGrid</name>
      <anchorfile>classTSOpt_1_1GridSpace.html</anchorfile>
      <anchor>a3f8b892ddfd2451bea49cd40d88f4c8b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isIncore</name>
      <anchorfile>classTSOpt_1_1GridSpace.html</anchorfile>
      <anchor>ae83bcfa5341448d0e3361ff74da1c1cc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1GridSpace.html</anchorfile>
      <anchor>ad8763d2f7748fa3b1044762ae309c5ef</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::MPIGridSpace</name>
    <filename>classTSOpt_1_1MPIGridSpace.html</filename>
    <member kind="function">
      <type></type>
      <name>MPIGridSpace</name>
      <anchorfile>classTSOpt_1_1MPIGridSpace.html</anchorfile>
      <anchor>aa53b65948a845e8298a61b01f286d8ab</anchor>
      <arglist>(string hdr, string dtype=&quot;notype&quot;, bool _incore=false, ostream &amp;_outfile=cerr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MPIGridSpace</name>
      <anchorfile>classTSOpt_1_1MPIGridSpace.html</anchorfile>
      <anchor>ae2fc89b4bfadb5588c2580d29638bde2</anchor>
      <arglist>(MPIGridSpace const &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~MPIGridSpace</name>
      <anchorfile>classTSOpt_1_1MPIGridSpace.html</anchorfile>
      <anchor>aacebc809e82cf4f2a703a633ea882554</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>LinearAlgebraPackage&lt; ireal &gt; const &amp;</type>
      <name>getLAP</name>
      <anchorfile>classTSOpt_1_1MPIGridSpace.html</anchorfile>
      <anchor>a6990bf8bcb619e226e7b140a09208cff</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DataContainerFactory const &amp;</type>
      <name>getDCF</name>
      <anchorfile>classTSOpt_1_1MPIGridSpace.html</anchorfile>
      <anchor>af2bc12819998dcf57be05c03cde210d1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>grid const &amp;</type>
      <name>getGrid</name>
      <anchorfile>classTSOpt_1_1MPIGridSpace.html</anchorfile>
      <anchor>abeae2edd5ac255bdbea2a1996fc46624</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isIncore</name>
      <anchorfile>classTSOpt_1_1MPIGridSpace.html</anchorfile>
      <anchor>a94aff9f9bbcd85867079f638f01a60f5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1MPIGridSpace.html</anchorfile>
      <anchor>acd6d17070433110dedd903f33670e248</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
</tagfile>
