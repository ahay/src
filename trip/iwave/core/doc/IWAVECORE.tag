<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="page">
    <name>index</name>
    <title>IWAVE Core Classes</title>
    <filename>index</filename>
  </compound>
  <compound kind="file">
    <name>bcast_params.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>bcast__params_8h</filename>
    <member kind="function">
      <type>int</type>
      <name>bcast_params</name>
      <anchorfile>bcast__params_8h.html</anchorfile>
      <anchor>a1cb021c13b1b0b78c8183edcee2ad4b1</anchor>
      <arglist>(PARARRAY *par, FILE *stream)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>create_sten_2k.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>create__sten__2k_8h</filename>
    <includes id="fd_8h" name="fd.h" local="yes" imported="no">fd.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>DEP_F</name>
      <anchorfile>create__sten__2k_8h.html</anchorfile>
      <anchor>a168538719b22838c100ed53c64eb9b63</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DEP_DFDZ</name>
      <anchorfile>create__sten__2k_8h.html</anchorfile>
      <anchor>adfc8b7d36e1cd3856b6a9312611b3c1e</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DEP_DFDX</name>
      <anchorfile>create__sten__2k_8h.html</anchorfile>
      <anchor>ab1ca0ba740a48da4c4bc8f74b55c265a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DEP_DFDY</name>
      <anchorfile>create__sten__2k_8h.html</anchorfile>
      <anchor>a9fab69499ec4f20994aeccc48db84dd9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DEP_LAPF</name>
      <anchorfile>create__sten__2k_8h.html</anchorfile>
      <anchor>a33edc9bdd420b82b70216d46551f062e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>create_sten2_2k</name>
      <anchorfile>create__sten__2k_8h.html</anchorfile>
      <anchor>aad2802a72f9b807c10368b4f0bdf99eb</anchor>
      <arglist>(FILE *stream, IWaveInfo const &amp;ic, int k, int ndim, IPNT gtype[RDOM_MAX_NARR], int stendm[RDOM_MAX_NARR][RDOM_MAX_NARR], STENCIL *sten)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>doc.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>doc_8h</filename>
  </compound>
  <compound kind="file">
    <name>doc_notes.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>doc__notes_8h</filename>
  </compound>
  <compound kind="file">
    <name>doc_parallel.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>doc__parallel_8h</filename>
  </compound>
  <compound kind="file">
    <name>exchange.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>exchange_8h</filename>
    <includes id="stencil_8h" name="stencil.h" local="yes" imported="no">stencil.h</includes>
    <includes id="rdomain_8h" name="rdomain.h" local="yes" imported="no">rdomain.h</includes>
    <member kind="function">
      <type>int</type>
      <name>ex_compute</name>
      <anchorfile>exchange_8h.html</anchorfile>
      <anchor>a7eeb5c5a93bf851e102ba9b889c7f5f4</anchor>
      <arglist>(int iarr, STENCIL *sten, RDOM *dom, IPNT a_gs, IPNT a_ge, IPNT r_gs[], IPNT r_ge[], int frcvempty[])</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>fd.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>fd_8h</filename>
    <includes id="rdomain_8h" name="rdomain.h" local="yes" imported="no">rdomain.h</includes>
    <includes id="exchange_8h" name="exchange.h" local="yes" imported="no">exchange.h</includes>
    <includes id="model_8h" name="model.h" local="yes" imported="no">model.h</includes>
    <includes id="iwinfo_8hh" name="iwinfo.hh" local="yes" imported="no">iwinfo.hh</includes>
    <member kind="define">
      <type>#define</type>
      <name>INCLUDE_BOUNDARY_PNTS</name>
      <anchorfile>fd_8h.html</anchorfile>
      <anchor>ab32b969abf783dfe30e94f152a549c88</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DUAL_GRID</name>
      <anchorfile>fd_8h.html</anchorfile>
      <anchor>ace5122de560069a28377cfc51568b66a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>PRIMAL_GRID</name>
      <anchorfile>fd_8h.html</anchorfile>
      <anchor>a036fb3e982456bbe0d4fd98f5763a5b2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>fd_update</name>
      <anchorfile>fd_8h.html</anchorfile>
      <anchor>a20321e89ae4de1a7dc31d70a24a42ccf</anchor>
      <arglist>(int ia, int iv, IWaveInfo const &amp;ic)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>fd_isdyn</name>
      <anchorfile>fd_8h.html</anchorfile>
      <anchor>aea7b96bfbe74eac692face4d89626377</anchor>
      <arglist>(int i, IWaveInfo const &amp;ic)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>fd_numsubsteps</name>
      <anchorfile>fd_8h.html</anchorfile>
      <anchor>a61722c532597b1abe3e6df4cce9243f5</anchor>
      <arglist>(IWaveInfo const &amp;ic)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>fd_readgrid</name>
      <anchorfile>fd_8h.html</anchorfile>
      <anchor>a204aebdb61214f9a4cb54f3b5ae49d7b</anchor>
      <arglist>(PARARRAY *par, FILE *stream, IMODEL *mdl, std::string gridkey)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>fd_modelcrea</name>
      <anchorfile>fd_8h.html</anchorfile>
      <anchor>ac5b2c5ef6d8dc738fbc0150866cc752f</anchor>
      <arglist>(IPNT cdims, IPNT crank, PARARRAY *par, FILE *stream, IMODEL *model, IWaveInfo const &amp;ic)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>istate.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>istate_8hh</filename>
    <includes id="iwinfo_8hh" name="iwinfo.hh" local="yes" imported="no">iwinfo.hh</includes>
    <includes id="iwave_8h" name="iwave.h" local="yes" imported="no">iwave.h</includes>
    <includes id="revolve_8h" name="revolve.h" local="yes" imported="no">revolve.h</includes>
    <class kind="class">TSOpt::IWaveSampler</class>
    <class kind="class">TSOpt::IWaveTree</class>
    <class kind="class">TSOpt::IWaveSim</class>
    <namespace>TSOpt</namespace>
    <member kind="function">
      <type>void</type>
      <name>IWaveEnvironment</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>a6ef420026e4d9483441a90f284f9c344</anchor>
      <arglist>(int argc, char **argv, int ts, PARARRAY **par, FILE **stream)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>IWaveApply</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>a4ec9561eac48757e078bcf7c66041d63</anchor>
      <arglist>(int argc, char **argv)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>iwave.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>iwave_8h</filename>
    <includes id="fd_8h" name="fd.h" local="yes" imported="no">fd.h</includes>
    <class kind="struct">s_PARALLELINFO</class>
    <class kind="struct">s_IWAVE</class>
    <member kind="typedef">
      <type>struct s_PARALLELINFO</type>
      <name>PARALLELINFO</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a109ecef7256427cf1fdd950e5e00fb7e</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct s_IWAVE</type>
      <name>IWAVE</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>aa92c3c9ad002003c549491d2c0174bdf</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>initparallel</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a4322a646407565091ac5d3041ad8ec0e</anchor>
      <arglist>(int ts)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initparallel_global</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>ad2ea5d7a20cb033ad75db0d6685f6154</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>initparallel_local</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a20c2e240e69db1c295bdea4ea31dc810</anchor>
      <arglist>(PARARRAY, FILE *)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>initpinfo</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>ab58f3736dc43109bfd7062f425899035</anchor>
      <arglist>(PARALLELINFO *pinfo, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>initexch</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a2bcf8c837b8d46fb69353312dd54fd0a</anchor>
      <arglist>(PARALLELINFO *pinfo, int ndim, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>destroypinfo</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a12e51cee9fdb6e5bdc5cbac72aed67f7</anchor>
      <arglist>(PARALLELINFO *pinfo, int fabort)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>iwave_construct</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>ae7f6ca0e69846e172c1b09432c86655b</anchor>
      <arglist>(IWAVE *state, PARARRAY *pars, FILE *stream, IWaveInfo const &amp;ic)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>iwave_printf</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a7ecca1b03428c60243b5f32341604a07</anchor>
      <arglist>(IWAVE *state, PARARRAY *pars, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>iwave_dynamic_init</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a9c4b35558629cce0291219c33acebd05</anchor>
      <arglist>(IWAVE *state, int it, IWaveInfo const &amp;ic)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>iwave_destroy</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>aae0df9f1ffc8e0054baf67176095c4fd</anchor>
      <arglist>(IWAVE *state, FD_MODELDEST d)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>quietexit</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a437e5233c92a52eca341ff5bc349727a</anchor>
      <arglist>(PARARRAY *pars, FILE **stream)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>abortexit</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a98fe4eb96342ca4fd0320a3cc515b413</anchor>
      <arglist>(int err, PARARRAY *pars, FILE **stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>initoutstream</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a14d372fb6a8836b98eebbc72e7d44e94</anchor>
      <arglist>(FILE **stream, int rk, int sz)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>readinput</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a4eda167ade0487f6ddd20fd645dd176a</anchor>
      <arglist>(PARARRAY **pars, FILE *stream, int argc, char **argv)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>readparpars</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a70752404f8530e755c518d44f61cb216</anchor>
      <arglist>(PARARRAY *pars, int *stats, int *nopts, int *printact, FILE *out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeparallel</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>a211e11d04814771ec7a5c4829823f4b4</anchor>
      <arglist>(PARALLELINFO *pinfo)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>setrecvexchange</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>ad5524a902569943bc13ac670870167d6</anchor>
      <arglist>(IMODEL *model, PARALLELINFO *pinfo, FILE *stream, IWaveInfo const &amp;ic)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>dump_pi</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>ae0b071fa81122ba9dc1c6e578db4cc4a</anchor>
      <arglist>(PARARRAY *pars, PARALLELINFO *pinfo, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dump_ac</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>aab0507fb4bfa1313d4dbbcf5f981ddc0</anchor>
      <arglist>(PARARRAY *pars, IMODEL *model, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dump_rs</name>
      <anchorfile>iwave_8h.html</anchorfile>
      <anchor>ad6d7cb725884e09bf7ce6cf6647f7ae3</anchor>
      <arglist>(PARARRAY *pars, IMODEL *model, FILE *stream, int sends)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>iwinfo.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>iwinfo_8hh</filename>
    <includes id="model_8h" name="model.h" local="yes" imported="no">model.h</includes>
    <class kind="struct">s_field</class>
    <class kind="struct">s_iokeys</class>
    <class kind="class">IWaveInfo</class>
    <class kind="struct">TSOpt::s_task_reln</class>
    <namespace>TSOpt</namespace>
    <member kind="define">
      <type>#define</type>
      <name>IWAVEMAXDATA</name>
      <anchorfile>iwinfo_8hh.html</anchorfile>
      <anchor>a3d8024458190d8c61f1b671937bd63d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>void(*</type>
      <name>FD_TIMESTEP</name>
      <anchorfile>iwinfo_8hh.html</anchorfile>
      <anchor>af86746c1676d4f2f4edb19e4f4462432</anchor>
      <arglist>)(std::vector&lt; RDOM * &gt; dom, bool fwd, int iv, void *fdpars)</arglist>
    </member>
    <member kind="typedef">
      <type>int(*</type>
      <name>FD_MODELINIT</name>
      <anchorfile>iwinfo_8hh.html</anchorfile>
      <anchor>a1a65e1ee5dcae43795427dec5d2b11a5</anchor>
      <arglist>)(PARARRAY *pars, FILE *stream, grid const &amp;g, ireal dt, void **specs)</arglist>
    </member>
    <member kind="typedef">
      <type>void(*</type>
      <name>FD_MODELDEST</name>
      <anchorfile>iwinfo_8hh.html</anchorfile>
      <anchor>ac4dcc79c60cc0245ad62564b27055543</anchor>
      <arglist>)(void **specs)</arglist>
    </member>
    <member kind="typedef">
      <type>int(*</type>
      <name>FD_TIMEGRID</name>
      <anchorfile>iwinfo_8hh.html</anchorfile>
      <anchor>afc8872ebc93b0a13c2f878d94aa714e6</anchor>
      <arglist>)(PARARRAY *pars, FILE *stream, grid const &amp;g, ireal &amp;dt)</arglist>
    </member>
    <member kind="typedef">
      <type>int(*</type>
      <name>FD_STENCIL</name>
      <anchorfile>iwinfo_8hh.html</anchorfile>
      <anchor>a27b9833ec8e1f9e42d128fe37fad6017</anchor>
      <arglist>)(void *specs, FILE *stream, int ndim, IPNT gtype[RDOM_MAX_NARR], STENCIL *sten)</arglist>
    </member>
    <member kind="typedef">
      <type>void(*</type>
      <name>FD_CHECK</name>
      <anchorfile>iwinfo_8hh.html</anchorfile>
      <anchor>adbb9133d8b5a78e76e9ec53f577c7215</anchor>
      <arglist>)(RDOM *dom, void *specs, FILE *stream)</arglist>
    </member>
    <member kind="typedef">
      <type>struct s_field</type>
      <name>FIELD</name>
      <anchorfile>iwinfo_8hh.html</anchorfile>
      <anchor>a1a24ba90e54929ea11199f782af4539f</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct s_iokeys</type>
      <name>IOKEY</name>
      <anchorfile>iwinfo_8hh.html</anchorfile>
      <anchor>aca7d8644dd669745ffd1d8ba57c165f0</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct TSOpt::s_task_reln</type>
      <name>TASK_RELN</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>ac7455791dcb5d897fc66c1ce818ea2f7</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>pow2</name>
      <anchorfile>iwinfo_8hh.html</anchorfile>
      <anchor>a37a984317bf282800383303959feffa9</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>IOTask</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>aed328f8ca196779d8b5c946f428b5a81</anchor>
      <arglist>(std::vector&lt; TASK_RELN * &gt; &amp;tr, int order, bool fwd, IWaveInfo const &amp;ic)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>IOTaskWriter</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>ab1833032111ef5b389f132a11ca5a5c9</anchor>
      <arglist>(std::vector&lt; TASK_RELN * &gt; const &amp;tr, ostream &amp;str)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>iwop.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>iwop_8hh</filename>
    <includes id="istate_8hh" name="istate.hh" local="yes" imported="no">istate.hh</includes>
    <class kind="class">TSOpt::IWaveSpace</class>
    <class kind="class">TSOpt::IWaveOp</class>
    <namespace>TSOpt</namespace>
    <member kind="define">
      <type>#define</type>
      <name>DEFAULT_SNAPS</name>
      <anchorfile>iwop_8hh.html</anchorfile>
      <anchor>a3e46a0090537c1b8126280e74bd10d65</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>model.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>model_8h</filename>
    <includes id="stencil_8h" name="stencil.h" local="yes" imported="no">stencil.h</includes>
    <includes id="rdomain_8h" name="rdomain.h" local="yes" imported="no">rdomain.h</includes>
    <class kind="struct">TIMESTEPINDEX</class>
    <class kind="struct">IMODEL</class>
    <member kind="typedef">
      <type>struct IMODEL</type>
      <name>IMODEL</name>
      <anchorfile>model_8h.html</anchorfile>
      <anchor>a197eafe0f299cd5c272fd8cfcadb2751</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>int(*</type>
      <name>IMODELINIT_FUN</name>
      <anchorfile>model_8h.html</anchorfile>
      <anchor>a9d0c8b90bddc03dc75953595076b580f</anchor>
      <arglist>)(PARARRAY *pars, FILE *stream, IMODEL *model)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>next_step</name>
      <anchorfile>model_8h.html</anchorfile>
      <anchor>a6364dbbf3b7a7c5324c2186a585bcec3</anchor>
      <arglist>(TIMESTEPINDEX *ts)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>less_than</name>
      <anchorfile>model_8h.html</anchorfile>
      <anchor>a037a39274b35a6d968654a5d63e52802</anchor>
      <arglist>(TIMESTEPINDEX t1, TIMESTEPINDEX t2)</arglist>
    </member>
    <member kind="function">
      <type>ireal</type>
      <name>get_time</name>
      <anchorfile>model_8h.html</anchorfile>
      <anchor>a920acb07cd90c1d2e7cab9a869283478</anchor>
      <arglist>(TIMESTEPINDEX ts)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>im_construct</name>
      <anchorfile>model_8h.html</anchorfile>
      <anchor>a9e7f76939ab99ad4dd09d3f4233e469d</anchor>
      <arglist>(IMODEL *model)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>im_destroy</name>
      <anchorfile>model_8h.html</anchorfile>
      <anchor>a6bf86c4826874b199714f4ff9aeed9c3</anchor>
      <arglist>(IMODEL *model, void(*destr)(void **))</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>im_setndim</name>
      <anchorfile>model_8h.html</anchorfile>
      <anchor>a4fbff36e81831162e188def1d837ac21</anchor>
      <arglist>(IMODEL *model)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>rdomain.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>rdomain_8h</filename>
    <class kind="struct">RDOM</class>
    <member kind="function">
      <type>int</type>
      <name>rd_a_setnull</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>aaa6e945d05bbda943ae34ef332010fe1</anchor>
      <arglist>(RDOM *dom)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_setnull</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a647d93160acbf7bfe51ca53751cfeb79</anchor>
      <arglist>(RDOM *dom, int iarr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_create_s</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a554189fbf7117b2e5be211e9789e6f49</anchor>
      <arglist>(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dn[])</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_create_e</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a319329b01c1753e7aed6b19c7bfb87b5</anchor>
      <arglist>(RDOM *dom, int narr, int ndim, IPNT dge[], IPNT dn[])</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_create</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>aba4b93b6bc62769206c5c6f0829fc70b</anchor>
      <arglist>(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dge[])</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_create_s</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>adaf90f02491c1a07df3b092d2fc92d6f</anchor>
      <arglist>(RDOM *dom, int ndim, IPNT gs, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_create_e</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a6cf51b424bd7e169375936ce14a14f07</anchor>
      <arglist>(RDOM *dom, int ndim, IPNT ge, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_create</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a980743fc97e25460caef482ac30f6893</anchor>
      <arglist>(RDOM *dom, int ndim, IPNT gs, IPNT ge)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_declare_s</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ade0acfb67bdc56f10b6e62d9f4a1c7b9</anchor>
      <arglist>(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dn[])</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_declare_e</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>aacaf74bfbc325bc00a35f4a9549be784</anchor>
      <arglist>(RDOM *dom, int narr, int ndim, IPNT dge[], IPNT dn[])</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_declare</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ac35a57be1da9c076f96d550985a11c48</anchor>
      <arglist>(RDOM *dom, int narr, int ndim, IPNT dgs[], IPNT dge[])</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_declare_s</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a7e1f318443f99b6ec651d54362d02238</anchor>
      <arglist>(RDOM *dom, int ndim, IPNT gs, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_declare_e</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a06034b8a175e2f41981ddb29d52c8077</anchor>
      <arglist>(RDOM *dom, int ndim, IPNT ge, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_declare</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a8325bd22558cdf7c9f2be39488d674e4</anchor>
      <arglist>(RDOM *dom, int ndim, IPNT gs, IPNT ge)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_allocate</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ad4f9ea93f5df5562bf74c7ebc8e6755c</anchor>
      <arglist>(RDOM *dom)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_allocate</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a58f56ac027f8724baa0a209b27451ab3</anchor>
      <arglist>(RDOM *dom, int iarr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_destroy</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>afd71f4c200b4060dc0825acc14d4b0c3</anchor>
      <arglist>(RDOM *dom)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_greset</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>abcf648513a1c321019171892c4149773</anchor>
      <arglist>(RDOM *dom, IPNT dgs[], IPNT dge[])</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_greset_s</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ab562a24b169e999dd572e20ac4770ed3</anchor>
      <arglist>(RDOM *dom, int iarr, IPNT gs, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_greset_e</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a28e3bc656f4277948cb3b0b301d886eb</anchor>
      <arglist>(RDOM *dom, int iarr, IPNT ge, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_greset</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ae67acd7550e554c4a1c98f7db5905290</anchor>
      <arglist>(RDOM *dom, int iarr, IPNT gs, IPNT ge)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_offset_s</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a62d992fcf404188716d4aa5854c9b126</anchor>
      <arglist>(RDOM *dom, int iarr, IPNT os, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_offset_e</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ab7a3ad4aef4c69a12cf14ee3feaade47</anchor>
      <arglist>(RDOM *dom, int iarr, IPNT oe, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_offset</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a6825501e7de35b58178da5e068601dea</anchor>
      <arglist>(RDOM *dom, int iarr, IPNT os, IPNT oe)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_dump</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a28ac9c69197e69413a240319718156fb</anchor>
      <arglist>(const RDOM *dom, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_dump</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a131e5b2b27b438be35d210b0f62996ae</anchor>
      <arglist>(const RDOM *dom, int iarr, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_print</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a8f0ddac987325c79b8d3d39242dc299e</anchor>
      <arglist>(RDOM *dom, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_fprint</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ad323921327625ab7b141e2db986a94ae</anchor>
      <arglist>(RDOM *dom, const char *path)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_fsprint</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a3fadcc57aa91863cf95fac9f9c1e6b7a</anchor>
      <arglist>(RDOM *dom, const char *path)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_print</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a88b520bfbfb9d3550fa9b2e0cd9ee52f</anchor>
      <arglist>(RDOM *dom, int iarr, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_fprint</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a0e1758e58656d99e77f07f7b1ff457c9</anchor>
      <arglist>(RDOM *dom, int iarr, const char *path)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_write</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a012f9d985e1fc91e94b884ab651a20f1</anchor>
      <arglist>(RDOM *dom, int iarr, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_fwrite</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ad33a8d6fbf1d7297d84d79fde69cd944</anchor>
      <arglist>(RDOM *dom, int iarr, const char *path)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_printslice</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a4a54cc16d135ac4514fa1e687aa62d0d</anchor>
      <arglist>(RDOM *dom, int iarr, FILE *stream, int idim, int islice)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_fprintslice</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a561de6499f20d1960502ef2976bf2cfe</anchor>
      <arglist>(RDOM *dom, int iarr, const char *path, int idim, int islice)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_writeslice</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ae0bb860dfa5abc7dfb9ad2abb68c986e</anchor>
      <arglist>(RDOM *dom, int iarr, FILE *stream, int idim, int islice)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_fwriteslice</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>abe5dcea3791be00b37e7f283577d7127</anchor>
      <arglist>(RDOM *dom, int iarr, const char *path, int idim, int islice)</arglist>
    </member>
    <member kind="function">
      <type>ireal</type>
      <name>rd_get</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a8e48e748996fd38349bdd8036af39c66</anchor>
      <arglist>(const RDOM *dom, int iarr, IPNT li)</arglist>
    </member>
    <member kind="function">
      <type>ireal</type>
      <name>rd_gget</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a196cb4fc951f958935fddacbc5c6b908</anchor>
      <arglist>(const RDOM *dom, int iarr, IPNT gi)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rd_set</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a095f506f28589b458dfc3f4787884bb6</anchor>
      <arglist>(RDOM *dom, int iarr, IPNT li, ireal r)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rd_gset</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ae45e48e1f607cc91d492f509f6aa3bfd</anchor>
      <arglist>(RDOM *dom, int iarr, IPNT gi, ireal r)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_size</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a6517371f4347390e9d860af93da038e7</anchor>
      <arglist>(RDOM *dom, int iarr, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_size</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ad711b8e2f014c7230b3c540cc128f7a4</anchor>
      <arglist>(RDOM *dom, int iarr, IPNT n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_gse</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>aee60a686de83ae8fddc64f1c238515d4</anchor>
      <arglist>(const RDOM *dom, int iarr, IPNT gs, IPNT ge)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_gse</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a8e2375f39d14d1fd3bf53660e65f74c2</anchor>
      <arglist>(const RDOM *dom, int iarr, IPNT gs, IPNT ge)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_ndim</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ae3d1c11e82cb480d02719f6413d66d22</anchor>
      <arglist>(const RDOM *dom, int iarr, int *ndim)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_narr</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a5ffd516a3f1b9551fafbfc2b6dc6b179</anchor>
      <arglist>(const RDOM *dom, int *narr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_setempty</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a15eb8d587d575704bbbac4828d0285fb</anchor>
      <arglist>(RDOM *dom, int iarr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_empty</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a179e6ad8529ffa8e315187b5c5a40019</anchor>
      <arglist>(RDOM *dom, int iarr, int *empty)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_setexchangeinfo</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>ae3451698a67b8e2d40bd8ee1aa11a7d3</anchor>
      <arglist>(RDOM *dom, int iarr, EXCHANGEINFO *einfo)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_overlap</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a1db06b649505d6e0334107ed5a98527f</anchor>
      <arglist>(RDOM *dom1, int iarr1, RDOM *dom2, int iarr2, int *overlap)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_setoverlap</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a1a92372de3b53efd9b5ff83fd30b0138</anchor>
      <arglist>(RDOM *dom1, int iarr1, RDOM *dom2, int iarr2)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_inner</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>a04fb6456da8e953fdf64afe8aa53de8b</anchor>
      <arglist>(RDOM const *dom1, RDOM const *dom2, ireal *ip)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_zero</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>abf6bcd397cc55c4ce4145a8de7bbb416</anchor>
      <arglist>(RDOM *dom)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>rd_a_scale</name>
      <anchorfile>rdomain_8h.html</anchorfile>
      <anchor>aaeb8810c60c55de893d1d0377c672faf</anchor>
      <arglist>(RDOM *dom, int iarr, ireal fac)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>revolve.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>revolve_8h</filename>
    <class kind="class">Checkpoint</class>
    <class kind="class">Schedule</class>
    <class kind="class">Online</class>
    <class kind="class">Online_r2</class>
    <class kind="class">Online_r3</class>
    <class kind="class">Arevolve</class>
    <class kind="class">Moin</class>
    <class kind="class">Offline</class>
    <class kind="class">Revolve</class>
    <namespace>ACTION</namespace>
    <member kind="define">
      <type>#define</type>
      <name>checkup</name>
      <anchorfile>revolve_8h.html</anchorfile>
      <anchor>a4d9bb7b76904058421e65974cc9afe40</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>repsup</name>
      <anchorfile>revolve_8h.html</anchorfile>
      <anchor>a72be6312c1deba95502160c5f28d0395</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MAXINT</name>
      <anchorfile>revolve_8h.html</anchorfile>
      <anchor>ada488b1a153e29f9ae0b098de0d6912d</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <name>action</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6f</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>advance</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fa5caa012e1a233bb4ceac641c540e3c6a</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>takeshot</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fa1194129978cb5215feb60a0a1991f817</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>restore</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fabfae685cfc5f1c90c1fca18580138db5</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>firsturn</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fa01684db13bde5f0a8ae7ade293978897</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>youturn</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fab281176f2370fd21fa55650a828e534d</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>terminate</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fad3102ce1cd56c9da8d72ee340a550e49</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>error</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fab6437f24a8d3da111127c1c25618e244</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>stencil.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>stencil_8h</filename>
    <class kind="struct">STENCIL_MASK</class>
    <class kind="struct">STENCIL</class>
    <member kind="function">
      <type>int</type>
      <name>mask_setnull</name>
      <anchorfile>stencil_8h.html</anchorfile>
      <anchor>a766d3bceff1c97997b3cc3f402e6a29b</anchor>
      <arglist>(STENCIL_MASK *mask)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>mask_create</name>
      <anchorfile>stencil_8h.html</anchorfile>
      <anchor>a6f2aea73fc339ab6e8c5a6b583964c5b</anchor>
      <arglist>(STENCIL_MASK *mask, int ip, int ir, int n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>mask_destroy</name>
      <anchorfile>stencil_8h.html</anchorfile>
      <anchor>a3df589e38e79f4e23d141782544b20e4</anchor>
      <arglist>(STENCIL_MASK *mask)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>mask_set</name>
      <anchorfile>stencil_8h.html</anchorfile>
      <anchor>a311a6473c8a6a5d3855741deff53bc72</anchor>
      <arglist>(STENCIL_MASK *mask, int i, const IPNT ind)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>mask_get</name>
      <anchorfile>stencil_8h.html</anchorfile>
      <anchor>a65bcaec639834f22f8cc8feb18c531b8</anchor>
      <arglist>(STENCIL_MASK *mask, int i, IPNT ind)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sten_setnull</name>
      <anchorfile>stencil_8h.html</anchorfile>
      <anchor>ad63ba423129fc29ae69a988892de898d</anchor>
      <arglist>(STENCIL *sten)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sten_create</name>
      <anchorfile>stencil_8h.html</anchorfile>
      <anchor>a10a69d12fd6639f93dd6d4fe590b123f</anchor>
      <arglist>(STENCIL *sten, int nmask)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sten_destroy</name>
      <anchorfile>stencil_8h.html</anchorfile>
      <anchor>a4df5897398a94be7875731ef144355f3</anchor>
      <arglist>(STENCIL *sten)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sten_set</name>
      <anchorfile>stencil_8h.html</anchorfile>
      <anchor>a76c4a670a2ddab2b74e5d99ca229c76a</anchor>
      <arglist>(STENCIL *sten, int imask, STENCIL_MASK *mask)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sten_get</name>
      <anchorfile>stencil_8h.html</anchorfile>
      <anchor>a04ac328c2f5bb6003867b7dfb2159847</anchor>
      <arglist>(STENCIL *sten, int imask, STENCIL_MASK *mask)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sten_out</name>
      <anchorfile>stencil_8h.html</anchorfile>
      <anchor>a66dad6c28a7a366534d1a2104e6d232a</anchor>
      <arglist>(STENCIL *sten, FILE *stream, const char *ind2str_fun(int))</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>usage_selfdoc.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/core/include/</path>
    <filename>usage__selfdoc_8h</filename>
    <member kind="function">
      <type>Typical parameter list May be and used for include parameters on command</type>
      <name>line</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>a3a322f8374010e0d578889b3149e30f6</anchor>
      <arglist>(for example in Flow)</arglist>
    </member>
    <member kind="variable">
      <type>Typical parameter list May be</type>
      <name>copied</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>a42c2135b42864057d1a356680edfaa82</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Typical parameter list May be</type>
      <name>edited</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>a4792c07677d965156b1a6113f554df59</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Typical parameter list May be and used for</type>
      <name>input</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>a7beac068c17687fb99c931200916f3ae</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Typical parameter list May be and used for include parameters on command or</type>
      <name>place</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>aa14742f388f275fab86e14cd99cd2b42</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Typical parameter list May be and used for include parameters on command or in file&lt; foo &gt; and include</type>
      <name>par</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>aaac2dc988f9a64da4fc9ccb32c85b66a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Typical parameter list May be and used for include parameters on command or in file&lt; foo &gt; and include included explicitly in command line overrides parameter with same</type>
      <name>key</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>a2127a947894a3f1910fc2e56ca1c69d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Typical parameter list May be and used for include parameters on command or in file&lt; foo &gt; and include included explicitly in command line overrides parameter with same in par</type>
      <name>file</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>a8c45d4c09018e137e2f648db5a994cae</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Typical parameter list May be and used for include parameters on command or in file&lt; foo &gt; and include included explicitly in command line overrides parameter with same in par Invoke single threaded execution of xxxx</type>
      <name>by</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>ace03d2e053b535b6f1bf92e979385727</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Typical parameter list May be and used for include parameters on command or in file&lt; foo &gt; and include included explicitly in command line overrides parameter with same in par Invoke single threaded execution of xxxx</type>
      <name>sfxxxx</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>ab6a52cbf7a3150b83c4846f81e337af9</anchor>
      <arglist>[parameters][Madagascar install]</arglist>
    </member>
    <member kind="variable">
      <type>Typical parameter list May be and used for include parameters on command or in file&lt; foo &gt; and include included explicitly in command line overrides parameter with same in par Invoke single threaded execution of xxxx</type>
      <name>or</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>a3fd0023df7f56212df75a3bbdfa85cd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Typical parameter list May be and used for include parameters on command or in file&lt; foo &gt; and include included explicitly in command line overrides parameter with same in par Invoke single threaded execution of xxxx xxxx</type>
      <name>x</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>a3e02e3b72a0ef684d25e3cae20eb6442</anchor>
      <arglist>[parameters][standalone install]</arglist>
    </member>
    <member kind="variable">
      <type>non optional values indicated corner</type>
      <name>brackets</name>
      <anchorfile>usage__selfdoc_8h.html</anchorfile>
      <anchor>a84cc1b6a7e179790032b922f8143ebc2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>notes</name>
    <title>Notes on Modeling with IWAVE</title>
    <filename>notes</filename>
  </compound>
  <compound kind="page">
    <name>parallel</name>
    <title>Parallel Simulations with IWAVE</title>
    <filename>parallel</filename>
  </compound>
  <compound kind="class">
    <name>Arevolve</name>
    <filename>classArevolve.html</filename>
    <base>Online</base>
    <member kind="function">
      <type></type>
      <name>Arevolve</name>
      <anchorfile>classArevolve.html</anchorfile>
      <anchor>a477e1eb0c7b52030c1728c24b4532e8a</anchor>
      <arglist>(int sn, Checkpoint *c, ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Arevolve</name>
      <anchorfile>classArevolve.html</anchorfile>
      <anchor>a35cc05ccd0a10251178b3727bcc43fd4</anchor>
      <arglist>(Arevolve &amp;o)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>tmin</name>
      <anchorfile>classArevolve.html</anchorfile>
      <anchor>a01644f4e8118a83d2da58e01fe76e573</anchor>
      <arglist>(int steps, int snaps, ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sumtmin</name>
      <anchorfile>classArevolve.html</anchorfile>
      <anchor>ac671251a7935449d6790c5d99022221b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>mintmin</name>
      <anchorfile>classArevolve.html</anchorfile>
      <anchor>ab3be74edd7aa104a89d7750e3a244a3e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_fine</name>
      <anchorfile>classArevolve.html</anchorfile>
      <anchor>af5e27f9d5e6a844fee4da34ef376b194</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>akt_cp</name>
      <anchorfile>classArevolve.html</anchorfile>
      <anchor>a2d0357ef54ec8fd8d532a8e750177a2c</anchor>
      <arglist>(int cp)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_fine</name>
      <anchorfile>classArevolve.html</anchorfile>
      <anchor>a087160017858759a29eaa6311b1b8e27</anchor>
      <arglist>(int f)</arglist>
    </member>
    <member kind="function">
      <type>enum ACTION::action</type>
      <name>revolve</name>
      <anchorfile>classArevolve.html</anchorfile>
      <anchor>a78c1877c2aff0c647aa3ee9dc3e52943</anchor>
      <arglist>(ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Arevolve</name>
      <anchorfile>classArevolve.html</anchorfile>
      <anchor>a7902d38d27dc9aff739c11b972c7c967</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Checkpoint</name>
    <filename>classCheckpoint.html</filename>
    <member kind="function">
      <type></type>
      <name>Checkpoint</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>a191a2cd6845e05a7e080e38410c73bd5</anchor>
      <arglist>(int s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_ch</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>aaa175d4ab2a01fbf78ca87f6f90cbf0e</anchor>
      <arglist>(ostream &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_ord_ch</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>ad9571036e895f50012f3f90b4f64101f</anchor>
      <arglist>(ostream &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_ord_ch</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>a3f166eb2e9b481bfd2e4caaddce17d40</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Checkpoint</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>a7fff1a4a6a71c2150a969fbe40d01c34</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>vector&lt; int &gt;</type>
      <name>ch</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>a92775bd341f8060a9ae9c153c54dabad</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>vector&lt; int &gt;</type>
      <name>ord_ch</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>af84ad611370aefd4bb4300a6056bbc78</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>vector&lt; int &gt;</type>
      <name>number_of_writes</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>af76f3a084deb3ec3ba271e85c2271112</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>vector&lt; int &gt;</type>
      <name>number_of_reads</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>a15313ebb370d0c7e63091de2d98cc4ab</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>advances</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>a79ae34d7f31d3dc90c197fcd58cd36ec</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>takeshots</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>a5faf6a4d9111b8d20c5e4263dbc70a68</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>commands</name>
      <anchorfile>classCheckpoint.html</anchorfile>
      <anchor>ad02592bc6fbb1088ffbee0ee17cfadc7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>IMODEL</name>
    <filename>structIMODEL.html</filename>
    <member kind="variable">
      <type>void *</type>
      <name>specs</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>ac1f42e1045d1c00262c45eff1f81999c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>grid</type>
      <name>g</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>a7464eadb75bc3fc0c5a8d380efef9597</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>grid</type>
      <name>gl</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>a699a83219b2349679539adc0a96290ac</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>IPNT</type>
      <name>nls</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>aa15a9dd5adc3a8d8ff0aad15c311e350</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>IPNT</type>
      <name>nrs</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>a52c8857911aacf68d54a8313593750c8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nnei</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>a48102bb60dff6db6f4d447d5b6d247c7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>TIMESTEPINDEX</type>
      <name>tsind</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>a6d545f4821423088d8f1f44acf4e517d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>RDOM</type>
      <name>ld_a</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>a0383bd4f1385129e7b8c632b81a7e5f9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>RDOM</type>
      <name>ld_c</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>a205a86e4f1593c0ac777ec3b67c24cd7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>RDOM</type>
      <name>ld_p</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>a66d52d0a4592f44c898bca00b089192a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>RDOM *</type>
      <name>ld_s</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>a4c7ad708da78259279f3cf5b5bb0a27a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>RDOM *</type>
      <name>ld_r</name>
      <anchorfile>structIMODEL.html</anchorfile>
      <anchor>a401ab9c4e2d896f2b39453ca314d5dbf</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IWaveInfo</name>
    <filename>classIWaveInfo.html</filename>
    <member kind="function">
      <type>std::string</type>
      <name>get_iwave_model</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a0dc8bbb2ddc0356b3898e95d68deeb60</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>IOKEY *</type>
      <name>get_iwave_iokeys</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>aa7ac5b7fb68c9b8f4657b7ddfcbb3083</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>FIELD *</type>
      <name>get_iwave_fields</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a1ec04f8334954afadca530522fdbd781</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>FD_MODELINIT</type>
      <name>get_minit</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a44299654fc2a4e4d50d5dca837ca4e6d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>FD_MODELDEST</type>
      <name>get_mdest</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a1cd966666f4fc2885a48e8d4364d26ca</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>FD_TIMEGRID</type>
      <name>get_timegrid</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a83a78bfb5c1e16e876842aa122892d5a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>FD_TIMESTEP</type>
      <name>get_timestep</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a3819e43d802e588d553fa79e56017d62</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>FD_STENCIL</type>
      <name>get_stencil</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>af1aad529f4c5f656d034636e74ae726f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>FD_CHECK</type>
      <name>get_check</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>aefb4760490d0d2f148ab3e2f3ac37ab0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_num_fields</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a87375a9cf94c1d76aa60defae529fe70</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_num_iokeys</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a05c0dc5441f65e19daa7cc1b7a4b16e4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write_iwave_fields</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a07c654c050f77676e368f98de628df5a</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write_iwave_iokeys</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>abbdc13fd56bb37dcadf41518d9208b06</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static std::string</type>
      <name>iwave_model</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>aa01576db6f2dc9e89727a900450cd6bf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static FIELD</type>
      <name>iwave_fields</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a2234e215ad3f4994d43da57375c9bb7d</anchor>
      <arglist>[]</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static IOKEY</type>
      <name>iwave_iokeys</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a1510eeffad9bab211b3bbceda5e3104a</anchor>
      <arglist>[]</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static FD_MODELINIT</type>
      <name>minit</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a8afc38e0207e95ec7197e70f3db8e25e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static FD_MODELDEST</type>
      <name>mdest</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a7886bb46086fc846af0f8c188db80a08</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static FD_TIMESTEP</type>
      <name>timestep</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a2b93f6fd7240609da7ab0052f85e9fd1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static FD_TIMEGRID</type>
      <name>timegrid</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>a015dc5971fdf5a6cc9472afa7841a672</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static FD_STENCIL</type>
      <name>createstencil</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>ab5719dc54b74e848909e37065573c099</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static FD_CHECK</type>
      <name>check</name>
      <anchorfile>classIWaveInfo.html</anchorfile>
      <anchor>aa16e1a9c8d9b0b5cae3834d2fe8ef02f</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Moin</name>
    <filename>classMoin.html</filename>
    <base>Online</base>
    <member kind="function">
      <type></type>
      <name>Moin</name>
      <anchorfile>classMoin.html</anchorfile>
      <anchor>aece1b001bc3d1f3bc25a707095830600</anchor>
      <arglist>(int sn, Checkpoint *c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Moin</name>
      <anchorfile>classMoin.html</anchorfile>
      <anchor>ad98afe8dee372b23ddb85b47d3a76e62</anchor>
      <arglist>(Moin &amp;o)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>is_dispensable</name>
      <anchorfile>classMoin.html</anchorfile>
      <anchor>a0094d8052af06cba5b89a27f7e44e993</anchor>
      <arglist>(int *index)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_lmin</name>
      <anchorfile>classMoin.html</anchorfile>
      <anchor>a749f9ae82fe8fa9fbe4dc61a7ab52625</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>adjust_cp</name>
      <anchorfile>classMoin.html</anchorfile>
      <anchor>a797a067a394d96f5338ae034177e61cf</anchor>
      <arglist>(int index)</arglist>
    </member>
    <member kind="function">
      <type>enum ACTION::action</type>
      <name>revolve</name>
      <anchorfile>classMoin.html</anchorfile>
      <anchor>a989803314de2df388d7122ba806422f7</anchor>
      <arglist>(ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_fine</name>
      <anchorfile>classMoin.html</anchorfile>
      <anchor>afebed87d7fdf86304b4622e378f46e4d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_fine</name>
      <anchorfile>classMoin.html</anchorfile>
      <anchor>a70d1df7542b35bd8320e8ed7e27c95c6</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Moin</name>
      <anchorfile>classMoin.html</anchorfile>
      <anchor>a4f6d381f0b6178abc1f625d871871c04</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Offline</name>
    <filename>classOffline.html</filename>
    <base>Schedule</base>
    <member kind="function">
      <type></type>
      <name>Offline</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a2d3220ad0108852ca04590958772cade</anchor>
      <arglist>(int st, int sn, Checkpoint *c, ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Offline</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>ac6dcc290191fd5b7eaf428fe379d7eec</anchor>
      <arglist>(int sn, Checkpoint *c, Online *o, int f, ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Offline</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>ac10c88a579bd89e8871a3774c7ff5f1c</anchor>
      <arglist>(Schedule *o, ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Offline</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a55fa5bbb9bf7954d69863299ffb9a791</anchor>
      <arglist>(Offline &amp;o)</arglist>
    </member>
    <member kind="function">
      <type>ACTION::action</type>
      <name>revolve</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a43d667d403c268cc3f35f32fa98ccb4b</anchor>
      <arglist>(ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_check</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>acb18087b52f127a815a3ce01a043d69e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_capo</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a27b1066f1e64d0ce54169b2ffa9e6f36</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_fine</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>ad9843c8e705f5047d3b905730a767bff</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_snaps</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a79e30655f87258823ce0ff79d7e5907f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_commands</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a24f36eaa5bc00303b1c6e777e3859f8a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_steps</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a64706455d3d8470b4bb0a85c8398c70b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_online</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>aa3fe55aec75663f09c108599b5b6cb10</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>vector&lt; int &gt;</type>
      <name>get_num_ch</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a91207fd7bec1622c93eed6378b093e2d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_num_ch</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a858310565a07d8221ea38d5ef6e381c1</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_fine</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a3eee76f6e9ba53a016fdd9801f4c5942</anchor>
      <arglist>(int f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_capo</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a598ffc03a36fd61828af6bac202fbaf1</anchor>
      <arglist>(int c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Offline</name>
      <anchorfile>classOffline.html</anchorfile>
      <anchor>a747c2340b6fa12ff100926a2d678ae34</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Online</name>
    <filename>classOnline.html</filename>
    <base>Schedule</base>
    <member kind="function">
      <type></type>
      <name>Online</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>aa85380fe158831898712eae84ae23aa9</anchor>
      <arglist>(int sn, Checkpoint *c, bool o)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Online</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>a300e883fb0aa2a43d5bb78de464b98d1</anchor>
      <arglist>(Online &amp;o)</arglist>
    </member>
    <member kind="function">
      <type>ACTION::action</type>
      <name>revolve</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>a194103a871ed7dcbc6e205b9410c1dc0</anchor>
      <arglist>(ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_check</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>aba3c22a23c819b0259bffaee62d46702</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_capo</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>aa6b49ad15e5152a464ab3832336c82ba</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_fine</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>a7ba4a274293b2573a11fefd25241ff5c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_fine</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>a47946740dbe9015401dd27ca790214e4</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_output</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>af2561d4bcacd775bba43c8a647cdce66</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_capo</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>a6e93856629c20696300ae23d58d6ccb4</anchor>
      <arglist>(int c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Online</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>abaf1ef34225654b1857ed73f33198f73</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>check</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>aa19e682cc36f837dd4211417cf3cbbb8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>output</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>a0bb41828a79b2ffc241bd37708e790c9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>capo</name>
      <anchorfile>classOnline.html</anchorfile>
      <anchor>adbc4abb7ec5a84dfc0a914067ed36e64</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Online_r2</name>
    <filename>classOnline__r2.html</filename>
    <base>Online</base>
    <member kind="function">
      <type></type>
      <name>Online_r2</name>
      <anchorfile>classOnline__r2.html</anchorfile>
      <anchor>abf9728525ac82f194e132f00d1166f06</anchor>
      <arglist>(int sn, Checkpoint *c, bool o)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Online_r2</name>
      <anchorfile>classOnline__r2.html</anchorfile>
      <anchor>a9076f08b155003826777dd8618fa1adf</anchor>
      <arglist>(Online &amp;o)</arglist>
    </member>
    <member kind="function">
      <type>ACTION::action</type>
      <name>revolve</name>
      <anchorfile>classOnline__r2.html</anchorfile>
      <anchor>aa65d5b7cab0b2dbb15b60984bf17e5ba</anchor>
      <arglist>(ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_check</name>
      <anchorfile>classOnline__r2.html</anchorfile>
      <anchor>a2dc94fe7329b0356e034210a2b25347c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_capo</name>
      <anchorfile>classOnline__r2.html</anchorfile>
      <anchor>a74526c8021afd07df962d74be7b3f1e3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_fine</name>
      <anchorfile>classOnline__r2.html</anchorfile>
      <anchor>a82178028d287557e8e1c2caf2cafbbaa</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_output</name>
      <anchorfile>classOnline__r2.html</anchorfile>
      <anchor>aeea825f571d32a4db7929a757a79f525</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_fine</name>
      <anchorfile>classOnline__r2.html</anchorfile>
      <anchor>a1ee10797a0e25fbd126fc93d40d57432</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Online_r2</name>
      <anchorfile>classOnline__r2.html</anchorfile>
      <anchor>a95f43d4c759ba82ecc9b3aa7311031dc</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Online_r3</name>
    <filename>classOnline__r3.html</filename>
    <base>Online</base>
    <member kind="function">
      <type></type>
      <name>Online_r3</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a8bf256f26b804082210b779fa8c669f0</anchor>
      <arglist>(int sn, Checkpoint *c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Online_r3</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a18ba9fcb6f185a7591fd125832b43e17</anchor>
      <arglist>(Online_r3 &amp;o)</arglist>
    </member>
    <member kind="function">
      <type>ACTION::action</type>
      <name>revolve</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>aec8e7d497eb947ff711238d4dd4c093e</anchor>
      <arglist>(ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_capo</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a65bdc9eb08d4bbc84e090f7971337d55</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_fine</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>ab19a499276e8cff80efe910243d62e23</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_fine</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a96698b5c1dcc621a0b835ad432e789d4</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>choose_cp</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a4a5d61992e5b6b2506c7b9db42b410e8</anchor>
      <arglist>(int number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>tdiff_akt</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a36b8ea405cb18e2a32fbc27279a1054e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>akt_cp</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a9b9f0e2505b8087ac17caecc2f8311c2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_tdiff_end</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>ae0d7eb6b2a47b0e43cc275c5a5256f30</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Online_r3</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>af0b5ef7bd381c301d6dd108c28b2860f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>forward</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a5d7c16f0bd32baf5a455c557424dd578</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>ind_now</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a32e3d65519603ef3a74f869fc198092e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>cp</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a3781e074dae72c0c4c71a8f5f2d47984</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>vector&lt; int &gt;</type>
      <name>ch3</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a4b032f5c3cf1a685b5555e09c373b6ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>vector&lt; int &gt;</type>
      <name>tdiff</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a26cdf6fb77e077379f9695fdd1a8199d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>vector&lt; int &gt;</type>
      <name>tdiff_end</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a5ed68a57d68fd460bc9c53d187479d0f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>vector&lt; bool &gt;</type>
      <name>cp_fest</name>
      <anchorfile>classOnline__r3.html</anchorfile>
      <anchor>a8588a4ec60b1eacbc672712c2e7f73b3</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>RDOM</name>
    <filename>structRDOM.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>narr</name>
      <anchorfile>structRDOM.html</anchorfile>
      <anchor>adfe905ce2b87fc5f633d0d189acfb622</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>RARR</type>
      <name>_s</name>
      <anchorfile>structRDOM.html</anchorfile>
      <anchor>a93311b63dd53fa8e794d7c33d5e0cd04</anchor>
      <arglist>[RDOM_MAX_NARR]</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Revolve</name>
    <filename>classRevolve.html</filename>
    <member kind="function">
      <type></type>
      <name>Revolve</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a7521f7c848016cc51265b8b53362da34</anchor>
      <arglist>(int st, int sn, ostream &amp;_str)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Revolve</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a722d87bc06aeb550902060797a8f9d60</anchor>
      <arglist>(int st, int sn, int sn_ram, ostream &amp;_str)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Revolve</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a1a0ffadefa1c5177a2c83677a21e8939</anchor>
      <arglist>(int sn, ostream &amp;_str)</arglist>
    </member>
    <member kind="function">
      <type>ACTION::action</type>
      <name>revolve</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a87d59c03ea625d5d7f051b2f38130f26</anchor>
      <arglist>(int *check, int *capo, int *fine, int snaps, int *info, bool *where_to_put)</arglist>
    </member>
    <member kind="function">
      <type>ACTION::action</type>
      <name>revolve</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a0f234600584d255d1d4c16fe1fbe7a32</anchor>
      <arglist>(int *check, int *capo, int *fine, int snaps, int *info)</arglist>
    </member>
    <member kind="function">
      <type>ACTION::action</type>
      <name>revolve</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a2503649ea0668d754b0d85c82b88af3c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>adjust</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>afafc52fed672e84f4423afe8c6182e9e</anchor>
      <arglist>(int steps, ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>maxrange</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>aed7c5aa144830c44fcb60062e9a82d1f</anchor>
      <arglist>(int ss, int tt, ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>expense</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>aa9579b4055e16e7c82206fa7dccae30a</anchor>
      <arglist>(int steps, int snaps, ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>numforw</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>af71474fcf9754c5d7ed53bcaa01769d4</anchor>
      <arglist>(int steps, int snaps, ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>turn</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>ad6cdf16e5a0e0e3612c52956ea618601</anchor>
      <arglist>(int fine)</arglist>
    </member>
    <member kind="function">
      <type>vector&lt; int &gt;</type>
      <name>get_write_and_read_counts</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a0448b4edb967ad3008c023df09ae46c5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_number_of_writes_i</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a2717864d626195b3ce25496dc478a511</anchor>
      <arglist>(int l, int c, int i)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_number_of_reads_i</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a77ae8e11d002c12827c5e00af70cf6d8</anchor>
      <arglist>(int l, int c, int i)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getadvances</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>ad6e0078c7b8c57fe5847f3e34f0529db</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getcheck</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a79755fad6a8ea9f96a76927884b34fe0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getcheckram</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a18119ed7d8be5ff14e62e92813dfe999</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getcheckrom</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>aa1017bac3402abfee5eba61fab960685</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getcapo</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>abf0a689bc1b44de48dee4207ab10d4bf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getfine</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a2c44d3d893bf9d63468e19ec8d79607b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getinfo</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>ae4638eaab059644f1776c28e12599074</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getoldcapo</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>afda59c654fba048752eb28e7d20dae34</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>getwhere</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a426b9a14bcb03332dca641d9e49180b2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_info</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a28a7ddba6061c54eb90af0b6acc0122a</anchor>
      <arglist>(int inf)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>number_of_writes</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a037d1b43b465eaa97cbcc6c133597c0c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>number_of_reads</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a58c34759d9b57ba91414e865acbdad4a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_r</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a81125a4086242df648a4b633cbf23541</anchor>
      <arglist>(int steps, int snaps, ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_r</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>a20b0cce5ba2ce2ef0e4955d551609a82</anchor>
      <arglist>(ostream &amp;str)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Revolve</name>
      <anchorfile>classRevolve.html</anchorfile>
      <anchor>ab5730b2aaf4b02a54fa94188bb151496</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>s_field</name>
    <filename>structs__field.html</filename>
    <member kind="variable">
      <type>std::string</type>
      <name>field</name>
      <anchorfile>structs__field.html</anchorfile>
      <anchor>a8a413ed0f966da5b0f987c094022782f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>dynamic</name>
      <anchorfile>structs__field.html</anchorfile>
      <anchor>a5bcb03d4e62601bef798721a10032171</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>substep</name>
      <anchorfile>structs__field.html</anchorfile>
      <anchor>a66c337a7bf855093dd9a736454a69204</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>IPNT</type>
      <name>gtype</name>
      <anchorfile>structs__field.html</anchorfile>
      <anchor>ad84ee8a64be4454b96bc9d07dca60013</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>s_iokeys</name>
    <filename>structs__iokeys.html</filename>
    <member kind="variable">
      <type>std::string</type>
      <name>keyword</name>
      <anchorfile>structs__iokeys.html</anchorfile>
      <anchor>ab4e8da04a2266056f68a543f5e10b1f1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>rarrindex</name>
      <anchorfile>structs__iokeys.html</anchorfile>
      <anchor>a7086f841855a79ac5f3d1ad6e26f3ea8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>input</name>
      <anchorfile>structs__iokeys.html</anchorfile>
      <anchor>ab41e237c00072ad295fbca1de654dca6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>active</name>
      <anchorfile>structs__iokeys.html</anchorfile>
      <anchor>ac6f4043c036839630f4f2797143e6c4e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>s_IWAVE</name>
    <filename>structs__IWAVE.html</filename>
    <member kind="variable">
      <type>PARALLELINFO</type>
      <name>pinfo</name>
      <anchorfile>structs__IWAVE.html</anchorfile>
      <anchor>a3a441e16df9384493649437df1093b43</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>IMODEL</type>
      <name>model</name>
      <anchorfile>structs__IWAVE.html</anchorfile>
      <anchor>ab1fababb36b37134f166a62f32f96b0a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>stats</name>
      <anchorfile>structs__IWAVE.html</anchorfile>
      <anchor>afdf5a8ca6d48ed50b4c2c305d0cfabdc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nopts</name>
      <anchorfile>structs__IWAVE.html</anchorfile>
      <anchor>aec1d349a33811bcdfd71b585f1c080a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>printact</name>
      <anchorfile>structs__IWAVE.html</anchorfile>
      <anchor>ab5fe8f250dee8bd0bd517b234aa46d5f</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>s_PARALLELINFO</name>
    <filename>structs__PARALLELINFO.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>ndim</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>a10bceab0c02d25485d0f91d990a7b45d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nnei</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>a7f9a5a67cc9de1901ee397947e3e74a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>wsize</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>abc576da7cc4f593c0fb228454ffc117b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>wrank</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>a944c75d0962c0b5020fb611539faccf0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>threadsupp</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>a71d8326d063f171937c737a008e7c5c8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>MPI_Comm</type>
      <name>ccomm</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>ae5281254eb98b068a3e128421497507c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>IPNT</type>
      <name>cdims</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>a782a490a2fa5436f1286de0b34956d56</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>IPNT</type>
      <name>crank</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>a4cd216b3df4c00b29ee59ab7465c2738</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>lrank</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>a4ada4217f03699fc7229e02edf3b5588</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>sranks</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>a2c4b2379f51b60a206abc0c841ee3db0</anchor>
      <arglist>[IWAVE_NNEI]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>rranks</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>a7be9c6a3f56676f2858e9da4c2f23fa2</anchor>
      <arglist>[IWAVE_NNEI]</arglist>
    </member>
    <member kind="variable">
      <type>EXCHANGEINFO</type>
      <name>seinfo</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>a7f1ce95a9fd738dc55df4ad7212c7940</anchor>
      <arglist>[RDOM_MAX_NARR][IWAVE_NNEI]</arglist>
    </member>
    <member kind="variable">
      <type>EXCHANGEINFO</type>
      <name>reinfo</name>
      <anchorfile>structs__PARALLELINFO.html</anchorfile>
      <anchor>a3136f651443c9c7a8a4b930f9105f724</anchor>
      <arglist>[RDOM_MAX_NARR][IWAVE_NNEI]</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Schedule</name>
    <filename>classSchedule.html</filename>
    <member kind="function">
      <type></type>
      <name>Schedule</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>a06305170f48c526ed3ece80db9c3ad72</anchor>
      <arglist>(int sn, Checkpoint *c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Schedule</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>ae072912e8434b080be40d004ac339ec5</anchor>
      <arglist>(int sn)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual ACTION::action</type>
      <name>revolve</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>aad6aac9f950f7e7b9bab1a045538178d</anchor>
      <arglist>(ostream &amp;str)=0</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>numforw</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>a6370f7184fe14314c11cb6a24f69a766</anchor>
      <arglist>(int steps, int snaps, ostream &amp;str)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>get_capo</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>a1d01cf0c9fae6a62b2c4b01ca8ccb032</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>get_fine</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>af1432856204a8d673d875bd4d59a96e7</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>get_check</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>a5e0a388cc9887905901e15d07dc137e1</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>set_fine</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>ae3cbfb9f4dcb8f5f0b70cd040ef4b9b3</anchor>
      <arglist>(int f)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>set_capo</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>aba03899153eb3f8dd4f8c72852ba2666</anchor>
      <arglist>(int c)=0</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_snaps</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>a8ec7e587ae4feac2178eaf9d020d76f9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Checkpoint *</type>
      <name>get_CP</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>a453dd1eba6ef53ad9ee7e933566cb246</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_advances</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>a8460f05e520123b31d01667c8c0d3db2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_shots</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>a31bf1bc41282605ffd7b4ebedb353f2c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_commands</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>aa7370d4ee9f51fe1bc9dc508f0dac8cf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>get_info</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>a264566e4ed2108c5a26fb552798c205a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Schedule</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>a4806b985197d35c00b9e707c0ed87998</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Checkpoint *</type>
      <name>checkpoint</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>ae3cf24eb75c1bfb66e9c346885a58f2d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>snaps</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>ae96f87a9e9878e840400e26ae3a0a438</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>info</name>
      <anchorfile>classSchedule.html</anchorfile>
      <anchor>a52605218e72aafa2d6b238ac151e1e50</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>STENCIL</name>
    <filename>structSTENCIL.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>nmask</name>
      <anchorfile>structSTENCIL.html</anchorfile>
      <anchor>a6b652b9509f12968b120ced2a838a3c1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>STENCIL_MASK *</type>
      <name>masks</name>
      <anchorfile>structSTENCIL.html</anchorfile>
      <anchor>a851735f78c97e94ebf78f027e69a638f</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>STENCIL_MASK</name>
    <filename>structSTENCIL__MASK.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>ip</name>
      <anchorfile>structSTENCIL__MASK.html</anchorfile>
      <anchor>a08ac7d7709299d02cf48c36f28f1ec24</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ir</name>
      <anchorfile>structSTENCIL__MASK.html</anchorfile>
      <anchor>a9931a3b54fc34e50d9bc8668d41330b8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>n</name>
      <anchorfile>structSTENCIL__MASK.html</anchorfile>
      <anchor>a990abb38439ee05e233cc8e50197ff8e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>IPNT *</type>
      <name>s</name>
      <anchorfile>structSTENCIL__MASK.html</anchorfile>
      <anchor>aeb7a3fc3550451b6a18f97a3cefb965f</anchor>
      <arglist></arglist>
    </member>
    <docanchor file="structSTENCIL__MASK">stencil</docanchor>
  </compound>
  <compound kind="struct">
    <name>TIMESTEPINDEX</name>
    <filename>structTIMESTEPINDEX.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>it</name>
      <anchorfile>structTIMESTEPINDEX.html</anchorfile>
      <anchor>a1c4ef504d5fb763465921b278a93d089</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>iv</name>
      <anchorfile>structTIMESTEPINDEX.html</anchorfile>
      <anchor>ae5fc9bc82c2767d2056d9fca07140dfe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>niv</name>
      <anchorfile>structTIMESTEPINDEX.html</anchorfile>
      <anchor>a274756be0d82c09ba541a158c056ff44</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ireal</type>
      <name>dt</name>
      <anchorfile>structTIMESTEPINDEX.html</anchorfile>
      <anchor>a09615d896959e957022ef9e8b15bf588</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>ACTION</name>
    <filename>namespaceACTION.html</filename>
    <member kind="enumeration">
      <name>action</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6f</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>advance</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fa5caa012e1a233bb4ceac641c540e3c6a</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>takeshot</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fa1194129978cb5215feb60a0a1991f817</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>restore</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fabfae685cfc5f1c90c1fca18580138db5</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>firsturn</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fa01684db13bde5f0a8ae7ade293978897</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>youturn</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fab281176f2370fd21fa55650a828e534d</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>terminate</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fad3102ce1cd56c9da8d72ee340a550e49</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>error</name>
      <anchorfile>namespaceACTION.html</anchorfile>
      <anchor>ac0ae89c835896822a9ef372f0e5e9b6fab6437f24a8d3da111127c1c25618e244</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>TSOpt</name>
    <filename>namespaceTSOpt.html</filename>
    <class kind="class">TSOpt::IWaveSampler</class>
    <class kind="class">TSOpt::IWaveTree</class>
    <class kind="class">TSOpt::IWaveSim</class>
    <class kind="struct">TSOpt::s_task_reln</class>
    <class kind="class">TSOpt::IWaveSpace</class>
    <class kind="class">TSOpt::IWaveOp</class>
    <member kind="typedef">
      <type>struct TSOpt::s_task_reln</type>
      <name>TASK_RELN</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>ac7455791dcb5d897fc66c1ce818ea2f7</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>IWaveEnvironment</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>a6ef420026e4d9483441a90f284f9c344</anchor>
      <arglist>(int argc, char **argv, int ts, PARARRAY **par, FILE **stream)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>IWaveApply</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>a4ec9561eac48757e078bcf7c66041d63</anchor>
      <arglist>(int argc, char **argv)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>IOTask</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>aed328f8ca196779d8b5c946f428b5a81</anchor>
      <arglist>(std::vector&lt; TASK_RELN * &gt; &amp;tr, int order, bool fwd, IWaveInfo const &amp;ic)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>IOTaskWriter</name>
      <anchorfile>namespaceTSOpt.html</anchorfile>
      <anchor>ab1833032111ef5b389f132a11ca5a5c9</anchor>
      <arglist>(std::vector&lt; TASK_RELN * &gt; const &amp;tr, ostream &amp;str)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::IWaveSampler</name>
    <filename>classTSOpt_1_1IWaveSampler.html</filename>
    <member kind="function">
      <type></type>
      <name>IWaveSampler</name>
      <anchorfile>classTSOpt_1_1IWaveSampler.html</anchorfile>
      <anchor>acaa36090e763736ce3ec3813c421062b</anchor>
      <arglist>(IWAVE *state, string key, PARARRAY &amp;pars, FILE *stream)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~IWaveSampler</name>
      <anchorfile>classTSOpt_1_1IWaveSampler.html</anchorfile>
      <anchor>ab0696169cd017827b9d7e4036dd71b09</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumAxes</name>
      <anchorfile>classTSOpt_1_1IWaveSampler.html</anchorfile>
      <anchor>a9fa2714c9bfdd4e1ae0fa8502c8057fe</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>axis const &amp;</type>
      <name>getAxis</name>
      <anchorfile>classTSOpt_1_1IWaveSampler.html</anchorfile>
      <anchor>a639eaebb237a962707f9f471fe59b816</anchor>
      <arglist>(int i) const </arglist>
    </member>
    <member kind="function">
      <type>ireal</type>
      <name>getCellVol</name>
      <anchorfile>classTSOpt_1_1IWaveSampler.html</anchorfile>
      <anchor>aeb653422170ad0d363be6493c15e36c0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ireal</type>
      <name>getRecipCellVol</name>
      <anchorfile>classTSOpt_1_1IWaveSampler.html</anchorfile>
      <anchor>a8503484361f49367584f740e1a09be38</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sample</name>
      <anchorfile>classTSOpt_1_1IWaveSampler.html</anchorfile>
      <anchor>a224f83e0a1d4a01d9092548da74fbfc1</anchor>
      <arglist>(grid g, IPNT step, bool fwd, bool input, IWAVE *state, int ridx, int iwdx, FILE *stream, bool dryrun=false, ostream &amp;drystr=cerr)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::IWaveTree</name>
    <filename>classTSOpt_1_1IWaveTree.html</filename>
    <member kind="function">
      <type></type>
      <name>IWaveTree</name>
      <anchorfile>classTSOpt_1_1IWaveTree.html</anchorfile>
      <anchor>a9c39d8db74eeb896207c4b2dc30cc686</anchor>
      <arglist>(PARARRAY &amp;_pars, FILE *_stream, IWaveInfo const &amp;_ic, int order=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IWaveTree</name>
      <anchorfile>classTSOpt_1_1IWaveTree.html</anchorfile>
      <anchor>ad18e9504a5f96b27060e1729c5fa9ee8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; IWAVE * &gt; &amp;</type>
      <name>getStateArray</name>
      <anchorfile>classTSOpt_1_1IWaveTree.html</anchorfile>
      <anchor>a69a2a231a76d48735ec32c1043d7294e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; IWAVE * &gt; const &amp;</type>
      <name>getStateArray</name>
      <anchorfile>classTSOpt_1_1IWaveTree.html</anchorfile>
      <anchor>a59eb8252a9cbabb5dda2d8b5fd836032</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; IWAVE * &gt; &amp;</type>
      <name>getRefStateArray</name>
      <anchorfile>classTSOpt_1_1IWaveTree.html</anchorfile>
      <anchor>a4c4d143dc6443c983eb48d243816fd54</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; IWAVE * &gt; const &amp;</type>
      <name>getRefStateArray</name>
      <anchorfile>classTSOpt_1_1IWaveTree.html</anchorfile>
      <anchor>af619637075ef7dba681db2d3bd1d7f27</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; RDOM * &gt; const &amp;</type>
      <name>getRDOMArray</name>
      <anchorfile>classTSOpt_1_1IWaveTree.html</anchorfile>
      <anchor>a28321778ef3ff9acb6a87f0fbc5d97b0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; RDOM * &gt; const &amp;</type>
      <name>getRefRDOMArray</name>
      <anchorfile>classTSOpt_1_1IWaveTree.html</anchorfile>
      <anchor>af865bc95c977386f35117a7f0f690540</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1IWaveTree.html</anchorfile>
      <anchor>a44b93339f4c672e21b67a2d04c071758</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::IWaveSim</name>
    <filename>classTSOpt_1_1IWaveSim.html</filename>
    <member kind="function">
      <type></type>
      <name>IWaveSim</name>
      <anchorfile>classTSOpt_1_1IWaveSim.html</anchorfile>
      <anchor>af47104e9a90a038718c79ca3f989903c</anchor>
      <arglist>(int order, bool fwd, PARARRAY &amp;par, FILE *stream, IWaveInfo const &amp;_ic, int printact=0, int snaps=0, bool dryrun=false, ostream &amp;drystr=cerr, ostream &amp;announce=cerr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IWaveSim</name>
      <anchorfile>classTSOpt_1_1IWaveSim.html</anchorfile>
      <anchor>a70d50958551f8362534ef6cd8af71cf0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>run</name>
      <anchorfile>classTSOpt_1_1IWaveSim.html</anchorfile>
      <anchor>a3ce00315fcbf4c5f568cb3a24743097f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; IWAVE * &gt; const &amp;</type>
      <name>getStateArray</name>
      <anchorfile>classTSOpt_1_1IWaveSim.html</anchorfile>
      <anchor>a796d865bcf1a61a7898c399fdc2e7a28</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; RDOM * &gt; const &amp;</type>
      <name>getRDOMArray</name>
      <anchorfile>classTSOpt_1_1IWaveSim.html</anchorfile>
      <anchor>ac01baa9729961686adf3eaf7cb460830</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>printgrid</name>
      <anchorfile>classTSOpt_1_1IWaveSim.html</anchorfile>
      <anchor>a1ab22fb117948fb6f049a60152bdd69a</anchor>
      <arglist>(FILE *fp) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1IWaveSim.html</anchorfile>
      <anchor>a6fb59483796d4289f26291c91fe5e8e5</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>TSOpt::s_task_reln</name>
    <filename>structTSOpt_1_1s__task__reln.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>iwaveindex</name>
      <anchorfile>structTSOpt_1_1s__task__reln.html</anchorfile>
      <anchor>aa2f2b88ba21d4e8e95b726f81d84b7f1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>keyword</name>
      <anchorfile>structTSOpt_1_1s__task__reln.html</anchorfile>
      <anchor>a928c1707f63054a764a5bf477dec029d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>rarrindex</name>
      <anchorfile>structTSOpt_1_1s__task__reln.html</anchorfile>
      <anchor>a33f613302f2c6baa5b8b01494f90d7bd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>input</name>
      <anchorfile>structTSOpt_1_1s__task__reln.html</anchorfile>
      <anchor>ae00096b32e740f7ea8fa4132ee59f7d9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::IWaveSpace</name>
    <filename>classTSOpt_1_1IWaveSpace.html</filename>
    <member kind="function">
      <type></type>
      <name>IWaveSpace</name>
      <anchorfile>classTSOpt_1_1IWaveSpace.html</anchorfile>
      <anchor>aeb58acaf4009fc1ecc1b5f00a726310d</anchor>
      <arglist>(PARARRAY const &amp;par, IWaveInfo const &amp;ic, bool input, ostream &amp;outfile=cerr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IWaveSpace</name>
      <anchorfile>classTSOpt_1_1IWaveSpace.html</anchorfile>
      <anchor>ab6b265fcabd5c65c973c0a0cd90b8c76</anchor>
      <arglist>(IWaveSpace const &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IWaveSpace</name>
      <anchorfile>classTSOpt_1_1IWaveSpace.html</anchorfile>
      <anchor>ac2e670482384456b469ee00a1e5ed8da</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DataContainer *</type>
      <name>buildDataContainer</name>
      <anchorfile>classTSOpt_1_1IWaveSpace.html</anchorfile>
      <anchor>aa114888de1fd3d4cfa0136b90719dd67</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSize</name>
      <anchorfile>classTSOpt_1_1IWaveSpace.html</anchorfile>
      <anchor>a8430efd6c0c47c50e007ed88849ddfdb</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; ireal &gt; const &amp;</type>
      <name>operator[]</name>
      <anchorfile>classTSOpt_1_1IWaveSpace.html</anchorfile>
      <anchor>acbf30501cb6555c089f19b6bd95f6294</anchor>
      <arglist>(size_t i) const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::string &gt;</type>
      <name>getKeys</name>
      <anchorfile>classTSOpt_1_1IWaveSpace.html</anchorfile>
      <anchor>ae2fb8ffad882c5b5e914a3b5db1c29b4</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TSOpt::IWaveOp</name>
    <filename>classTSOpt_1_1IWaveOp.html</filename>
    <member kind="function">
      <type></type>
      <name>IWaveOp</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a446d257d3268683b9495209cdbf5ea54</anchor>
      <arglist>(PARARRAY _pars, FILE *_stream, bool _dryrun=false, ostream &amp;_drystr=cerr, ostream &amp;_announce=cerr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IWaveOp</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a91c94ce33f8459c6abd0aee780cc5808</anchor>
      <arglist>(IWaveOp const &amp;x)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IWaveOp</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a786522fa186215f1d51ee15bc71e0763</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const IWaveSpace &amp;</type>
      <name>getIWaveDomain</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a4b55c9ad5ff00b47cd5eb3ad1b2d1dae</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const IWaveSpace &amp;</type>
      <name>getIWaveRange</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a89d4c5f2d5f7e5ec69a64672ad62b469</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; ireal &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>ab8e053cd6b9618a9d8f591b3d1961faa</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; ireal &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a3dd7e77143429cd02100fe2589ef7cc5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>PARARRAY &amp;</type>
      <name>getPar</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a44521e5f6cfdc740e367ec777a23698f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>PARARRAY const &amp;</type>
      <name>getPar</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a3f1630a931e9f31499acb37a1ac97f21</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>aa1d42958480405d61a5e6da742025ac2</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a6c5b70537a26793f836a4a17aba876e8</anchor>
      <arglist>(const Vector&lt; ireal &gt; &amp;x, Vector&lt; ireal &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a32a71829e7903519a7e9156d4098cbec</anchor>
      <arglist>(const Vector&lt; ireal &gt; &amp;x, const Vector&lt; ireal &gt; &amp;dx, Vector&lt; ireal &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>ab760a64c516832aaff03dbcca28a7efe</anchor>
      <arglist>(const Vector&lt; ireal &gt; &amp;x, const Vector&lt; ireal &gt; &amp;dy, Vector&lt; ireal &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv2</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a3de56630620ce22b469a421b19f45951</anchor>
      <arglist>(const Vector&lt; ireal &gt; &amp;x, const Vector&lt; ireal &gt; &amp;dx0, const Vector&lt; ireal &gt; &amp;dx1, Vector&lt; ireal &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv2</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>ac32cb098d7ed634dedf38073a9861e5c</anchor>
      <arglist>(const Vector&lt; ireal &gt; &amp;x, const Vector&lt; ireal &gt; &amp;dx0, const Vector&lt; ireal &gt; &amp;dy, Vector&lt; ireal &gt; &amp;dx1) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; ireal &gt; *</type>
      <name>clone</name>
      <anchorfile>classTSOpt_1_1IWaveOp.html</anchorfile>
      <anchor>a3bceac2e7c84bd3b03de78df7658a26e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
</tagfile>
