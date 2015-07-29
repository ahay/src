<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="page">
    <name>index</name>
    <title>IWAVE Basic Utilities Package</title>
    <filename>index</filename>
  </compound>
  <compound kind="file">
    <name>cstd.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/base/include/</path>
    <filename>cstd_8h</filename>
  </compound>
  <compound kind="file">
    <name>doc.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/base/include/</path>
    <filename>doc_8h</filename>
  </compound>
  <compound kind="file">
    <name>iwave_fopen.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/base/include/</path>
    <filename>iwave__fopen_8h</filename>
    <includes id="utils_8h" name="utils.h" local="yes" imported="no">utils.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>UNLINK_TMPS</name>
      <anchorfile>iwave__fopen_8h.html</anchorfile>
      <anchor>ac2bfc19430003fbd7e4db6b369aa379b</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>FILE *</type>
      <name>iwave_fopen</name>
      <anchorfile>iwave__fopen_8h.html</anchorfile>
      <anchor>a8f146c5713bad5d63f990c73d3729e2e</anchor>
      <arglist>(char **name, const char *mode, const char *proto, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>FILE *</type>
      <name>iwave_const_fopen</name>
      <anchorfile>iwave__fopen_8h.html</anchorfile>
      <anchor>ae1d3cc1a9ba0ca8a2e6aabdb97751b3c</anchor>
      <arglist>(const char *name, const char *mode, const char *proto, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>iwave_fclose</name>
      <anchorfile>iwave__fopen_8h.html</anchorfile>
      <anchor>abab5e5b217ac5a9b0186250f8557de77</anchor>
      <arglist>(FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>iwave_fdestroy</name>
      <anchorfile>iwave__fopen_8h.html</anchorfile>
      <anchor>ac8c88db3cc017c7697c3d74d652e5a3c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>iwave_fprintall</name>
      <anchorfile>iwave__fopen_8h.html</anchorfile>
      <anchor>a605d28ef6feaead82d779f6d627afba0</anchor>
      <arglist>(FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>const char *</type>
      <name>iwave_getproto</name>
      <anchorfile>iwave__fopen_8h.html</anchorfile>
      <anchor>a89b8a90d07bc39403dd4603ee0aa638b</anchor>
      <arglist>(const char *name)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>parser.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/base/include/</path>
    <filename>parser_8h</filename>
    <includes id="utils_8h" name="utils.h" local="yes" imported="no">utils.h</includes>
    <includes id="iwave__fopen_8h" name="iwave_fopen.h" local="yes" imported="no">iwave_fopen.h</includes>
    <class kind="struct">s_WORD</class>
    <class kind="struct">s_KEYVAL</class>
    <class kind="struct">s_PSLINK</class>
    <class kind="struct">s_PARARRAY</class>
    <namespace>RVL</namespace>
    <member kind="typedef">
      <type>struct s_WORD</type>
      <name>WORD</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a87900bbfb8b58df6c578788019c1273c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct s_KEYVAL</type>
      <name>KEYVAL</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>aa5fde07346236a20295bd6627cc6e69e</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct s_PSLINK</type>
      <name>PSLINK</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a2eee17f7254ce6ac97c378725eccd878</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct s_PARARRAY</type>
      <name>PARARRAY</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a46b6c542c3daf8f458e93bf556bb3d72</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>WORD *</type>
      <name>word_new</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>ad314cd6934ec82ddae91cd9cd3b0106a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>word_delete</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a60aa52de796242d8617032078ddd2caa</anchor>
      <arglist>(WORD **w)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>word_reset</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a943c53679a8f1e19337baf8c8c9631f1</anchor>
      <arglist>(WORD *w)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>word_assign</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a8a80dfe7ae836810e7652e1172ae1e19</anchor>
      <arglist>(WORD *w, const char *str, int len)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>word_whitechar</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>ade76196c336bca6eda5cd59978fc97d9</anchor>
      <arglist>(char c)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>word_copy</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a3ec0967bcefa5709fb0183179d17de32</anchor>
      <arglist>(WORD *tgt, WORD src)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>word_read</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a31b9492f09129836f1123be1cda34553</anchor>
      <arglist>(WORD *w, char **src)</arglist>
    </member>
    <member kind="function">
      <type>KEYVAL *</type>
      <name>kv_new</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>ae47ac8099cbff8460c92f7b04c298041</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>kv_delete</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>aeb9ca2260972c2a7b7443cfbce2ac641</anchor>
      <arglist>(KEYVAL **pair)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>kv_reset</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>afb2fede6c63343d37d8ff3df501e1e54</anchor>
      <arglist>(KEYVAL *pair)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>kv_check</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a712ecf61cd7af15d640c9d044e9cc277</anchor>
      <arglist>(KEYVAL src)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>kv_copy</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a081a54565b0d07115345d81ef64bd439</anchor>
      <arglist>(KEYVAL *tgt, KEYVAL src)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>kv_read</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>afc8e9e18651fbff935db63acaedad101</anchor>
      <arglist>(KEYVAL *kv, char **src)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>kv_print</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a0a9b2716adc5cceb8765d506a4546725</anchor>
      <arglist>(KEYVAL kv)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>kv_fprint</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a7b0606cb7749171c48a827e78cd3d24c</anchor>
      <arglist>(KEYVAL kv, FILE *f)</arglist>
    </member>
    <member kind="function">
      <type>PSLINK *</type>
      <name>pslink_new</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a527a991cbf0f765b3e2f44f15056205f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pslink_delete</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>ad260f4dbde17377184638f5ab888f1ff</anchor>
      <arglist>(PSLINK **p)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pslink_setnull</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a5b19a3a4ba6fafcaaae7940bcaaeb6e7</anchor>
      <arglist>(PSLINK **p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>pslink_front</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a5daf07c36a877015a64d499d1b507296</anchor>
      <arglist>(PSLINK **p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>pslink_back</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a7aa9d681b54892f95d37b17f7b483b17</anchor>
      <arglist>(PSLINK **p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>pslink_read</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>ae74c6b6150939065016a4a78aada8e10</anchor>
      <arglist>(PSLINK **p, char **str)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>pslink_findfirst</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>a02a1d0a3c306273b4232dddea6579ec5</anchor>
      <arglist>(PSLINK *par, WORD skey, WORD *sval)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>pslink_findlast</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>aecf14000ea87713af0583fa956b73a9f</anchor>
      <arglist>(PSLINK *par, WORD skey, WORD *sval)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>pslink_setfirst</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>aa1b304297a2ed79e751b7f85a653cc80</anchor>
      <arglist>(PSLINK **par, WORD skey, WORD sval)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>pslink_setlast</name>
      <anchorfile>parser_8h.html</anchorfile>
      <anchor>aecb9237b415ba71db08cb1a1b2c16109</anchor>
      <arglist>(PSLINK **par, WORD skey, WORD sval)</arglist>
    </member>
    <member kind="function">
      <type>PARARRAY *</type>
      <name>ps_new</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga22610c013d9baddbff98547d73b959f3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ps_delete</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gae7588685e9fcf8ea7c70bc421c896f43</anchor>
      <arglist>(PARARRAY **p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_setnull</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga403dda40e717a2a6fde0ce23f0bd7972</anchor>
      <arglist>(PARARRAY *parr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_printall</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>ga64ce217021050379834f6f38c5987101</anchor>
      <arglist>(PARARRAY parr, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_createfp</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>ga7aaefb64870ef48bbbaae5326bd39e73</anchor>
      <arglist>(PARARRAY *parr, FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_createfile_fproto</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>ga1a231cf9b18762847f072bf86d4e3dfe</anchor>
      <arglist>(PARARRAY *parr, FILE **stream, const char *proto, const char *filename)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_createfile</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>gafd1f9e4d4112cd2a5a1a358f536747b7</anchor>
      <arglist>(PARARRAY *parr, const char *filename)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_createargs</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>gaed3ad5b588d5854b59fc2ddc4fd43b62</anchor>
      <arglist>(PARARRAY *parr, int argc, char **argv)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_copy</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>gad44b69c03bec322c9e87f992cb4f4ea9</anchor>
      <arglist>(PARARRAY **tgt, PARARRAY src)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffcstring</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>gaf826510b46fad583252705525f124672</anchor>
      <arglist>(PARARRAY parr, const char *key, char **p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffchar</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga0ae7a683ffb52fa1e178d0d177e1cd7e</anchor>
      <arglist>(PARARRAY parr, const char *key, char *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffshort</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga0bfd6c5d1c6c1e14fdc6d47982df927a</anchor>
      <arglist>(PARARRAY parr, const char *key, short *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffint</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga4a5876ce56fbb1b103a524c863d7792d</anchor>
      <arglist>(PARARRAY parr, const char *key, int *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_fflong</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>gae7cbab0c40ca721c6bdc6aa5b872747e</anchor>
      <arglist>(PARARRAY parr, const char *key, long *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffushort</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga615269e94a42c65015a0e87ef2054582</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned short *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffuint</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga1d1b56b5c891de96dd19ea5ba9b3c14d</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned int *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffulong</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga5cc39d8804c55f4ef86e4bca41c314ac</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned long *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_fffloat</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>gacbc0af2c3b9aa11c95665344fb1ac0e7</anchor>
      <arglist>(PARARRAY parr, const char *key, float *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffdouble</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga17e18df3535d2e580af71e6693510560</anchor>
      <arglist>(PARARRAY parr, const char *key, double *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffreal</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga7e403f404da8e23ec947139f83113aac</anchor>
      <arglist>(PARARRAY parr, const char *key, ireal *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flcstring</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>gad600512aa48c23bca168afbd56360759</anchor>
      <arglist>(PARARRAY parr, const char *key, char **p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flchar</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga5e96dc7ee7d414859d6d91c925e2de00</anchor>
      <arglist>(PARARRAY parr, const char *key, char *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flshort</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>gabf116354c182aa9b1afd22208eba244d</anchor>
      <arglist>(PARARRAY parr, const char *key, short *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flint</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga5c39554369dfd58ba1c78f7511a3dea8</anchor>
      <arglist>(PARARRAY parr, const char *key, int *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_fllong</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga666f5e165f499ef37b9af543998dcf16</anchor>
      <arglist>(PARARRAY parr, const char *key, long *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flushort</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga23e235903067cc52b6b02f43cd8fb70a</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned short *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_fluint</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga57e5e102e3cc0687b55c0214f676ec1c</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned int *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flulong</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga07df22608785d6d6fc961db32ed34bbd</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned long *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flfloat</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>gaf47a79ee9eeb11f468154c0ff2c3a0fb</anchor>
      <arglist>(PARARRAY parr, const char *key, float *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_fldouble</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga5c692ae40d51f5c1300e05e0b9357a52</anchor>
      <arglist>(PARARRAY parr, const char *key, double *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flreal</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>gac460bbe480477a492b8f49eb48476e6c</anchor>
      <arglist>(PARARRAY parr, const char *key, ireal *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfcstring</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga8e8a1042527967504e4fca64f71d9bef</anchor>
      <arglist>(PARARRAY parr, const char *key, const char *val)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfchar</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gae97418a9c1f5515024d8e29ea39c6424</anchor>
      <arglist>(PARARRAY parr, const char *key, char p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfshort</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga74571071be0305000f3ce6a43c1101e1</anchor>
      <arglist>(PARARRAY parr, const char *key, short p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfint</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga30d103654f68511eca81ced1ffbf0e69</anchor>
      <arglist>(PARARRAY parr, const char *key, int p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sflong</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gab06af338eeaffc655e817a01dfec951f</anchor>
      <arglist>(PARARRAY parr, const char *key, long p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfushort</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga4278e02d6089c43797d40fd4ba1e1a9d</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned short p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfuint</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga3f6da9922f01b68cfac2a3f398c857bf</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned int p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfulong</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga4fd84b864795ef403a8575cfc0ac3a15</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned long p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sffloat</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gaf4f62c7b33f05bd632b6252058a3352a</anchor>
      <arglist>(PARARRAY parr, const char *key, float p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfdouble</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga04b714f7fd7aeaab2afb1c5f45978cdc</anchor>
      <arglist>(PARARRAY parr, const char *key, double p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfreal</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga105a6eb9a20b1847f5cfcf56edc9dc7a</anchor>
      <arglist>(PARARRAY parr, const char *key, ireal p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slcstring</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga47ed5c5bb39d5579e2c0fe47ef147dc9</anchor>
      <arglist>(PARARRAY parr, const char *key, const char *val)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slchar</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga833aab3c778e10765b133974a7470102</anchor>
      <arglist>(PARARRAY parr, const char *key, char p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slshort</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga66a745782546f06eadb3b627a7e706e5</anchor>
      <arglist>(PARARRAY parr, const char *key, short p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slint</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gad9d8015cab8bc775997486106645a184</anchor>
      <arglist>(PARARRAY parr, const char *key, int p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sllong</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga0e8695938688e7c8ca9b043c6491add9</anchor>
      <arglist>(PARARRAY parr, const char *key, long p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slushort</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gab16848ce5d857bfc249685abab98490d</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned short p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sluint</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga37784810b2472c56419b515d51dcb766</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned int p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slulong</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gac279ebee6e5509e3687a15bc36050486</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned long p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slfloat</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gac797f58997f7b291c34702fd39cbc9d8</anchor>
      <arglist>(PARARRAY parr, const char *key, float p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sldouble</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gabfe9137bdf2fafb0e68c3cb6ebda6673</anchor>
      <arglist>(PARARRAY parr, const char *key, double p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slreal</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga30aa66866932e4e69cfbb35f0feed844</anchor>
      <arglist>(PARARRAY parr, const char *key, ireal p)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a3c0789d4f7c051d3ce0ebcb99f0e9143</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, string &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gac8ed76b18b430feb617bd2abe41ceda2</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, char &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga3900c5e5488602f104b1e1d4fad95d83</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, short &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga94ae60450fd62a060e9657df7f64d18a</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, int &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga4b052205a4bae64c413a435afdeb68e6</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, long &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga21407734aa4593653ee79ad686f49fe6</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, unsigned short &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga942bd570653f63468d6ad27db005c54f</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, unsigned int &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga408d0b569973fed89d3532c1b71ebcc3</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, unsigned long &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gac72841c30b2e09f7997527e9896f7168</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, float &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gafc15e3471a2b1ab7e13cbcc8947e65b6</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, double &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga9c3dbb432ad3b6448de1858a6f64e39b</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, bool &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>valparse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gaccf70be59e5711495e5933b4e092d95c</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, T def)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>valparse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gaae0f6b6aad214dc4aea9c12af739f833</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>usempi.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/base/include/</path>
    <filename>usempi_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>MPICH_SKIP_MPICXX</name>
      <anchorfile>usempi_8h.html</anchorfile>
      <anchor>a077ffa213bd6af61ce94d0b0e1d4942a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MPI_Comm</name>
      <anchorfile>usempi_8h.html</anchorfile>
      <anchor>ae0a6553e8d5ccdb943ae1b816e4ac5ee</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MPI_PROC_NULL</name>
      <anchorfile>usempi_8h.html</anchorfile>
      <anchor>a62f1a2971d72ca736c9997d98568ff73</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MPI_Datatype</name>
      <anchorfile>usempi_8h.html</anchorfile>
      <anchor>ada6230d3c1696537923dc5a9446d7521</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MPI_DATATYPE_NULL</name>
      <anchorfile>usempi_8h.html</anchorfile>
      <anchor>a248e55efe62ab7ba021822dbe63e2cf2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>utils.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/iwave/base/include/</path>
    <filename>utils_8h</filename>
    <includes id="cstd_8h" name="cstd.h" local="yes" imported="no">cstd.h</includes>
    <includes id="usempi_8h" name="usempi.h" local="yes" imported="no">usempi.h</includes>
    <class kind="struct">s_SIZEDSTRING</class>
    <member kind="define">
      <type>#define</type>
      <name>inline</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a00d24c7231be28dbaf71f5408f30e44c</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>restrict</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a080abdcb9c02438f1cd2bb707af25af8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DT_CSTRING</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a1ffdafef8e5ad392c64f9ce46f71eb26</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DT_CHAR</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a5c6c994e380c763b5e2a5cab960e1eb4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DT_SHORT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a6858f84fb4436265ba11ab6726bbe0e3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DT_INT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a4e79b9fbe659e8934ee69b723cea4ee1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DT_LONG</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a3e996ab40386ca1965152452a38df738</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DT_FLOAT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>af1b7782f3dbb8e4b80a33f0d8cc4dd23</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DT_DOUBLE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a9ed82f0ea030c565c495ba7884d91a09</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DT_USHORT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a741be5df831b796ebb0968d26eeec39d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DT_UINT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a07fd3b7713c0d75657763a38db785d40</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DT_ULONG</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>aff6072d0efebcfc749f43acc48c8e2c4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>DT_REAL</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a7e39c895d11815f969806841e2d81ece</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>ireal</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a4442f038871a1a9a574595e979c6dd75</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>REAL_NAN</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ae32666337a62b7c6ae17bab939752dbd</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>REAL_ZERO</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ab8551802efb9fcdd01a8a07aec737c2e</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>REAL_ONE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>add4e84af378b21ff6824d144f0350c49</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>REAL_EPS</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>aebb5e6716e06431296af4d1a71744dec</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>REAL_MIN</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ae19c49417e2720027ed571edc915b5f2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>REAL_MAX</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a70fc4e60483cc1a5cb39dd935640cadc</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IWAVE_MPI_REAL</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ab13f97b05052638ab128089cfeb560f8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>RDOM_MAX_NARR</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a688b02c742527982430c7e70d68a8356</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>RARR_MAX_NDIM</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ad3923fd54ea972e213c1b8f1c2b88de1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IWAVE_NDIM</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a82680e7f95a8240216b9f786f5da9f2c</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IWAVE_3NDIM</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a3124a1a3cdf062ace8219c59097eaaba</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IWAVE_NNEI</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a3f7bb4ebfd73dcee5743b6aa58b11a6f</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IX</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a2e6a5e764791eb0062952229677be6c9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IY</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ad2b017eaecdc36ccf2a99d97fc53eff4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IZ</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ae00182be9e84cfaa674e420a191dfbbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>RARR_DUMP_POINTERS</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a09c52b84b35c1877d24564424a22524b</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IARR_MAX_NDIM</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>add356cb58441fc9b688e7ce7da88aa22</anchor>
      <arglist></arglist>
      <docanchor file="utils_8h">IPNT</docanchor>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>RARR_MAX_3NDIM</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a415bd0d4a9d3660e73dd438c8ffeb545</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>_IPNT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a08fb5aa043c3b6a7ca5b34404797e62c</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>_RPNT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>acc5115b3404495dc4d5e802105c21888</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>CHECK_BOUNDS</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a170752d39562f5d87a24111715e76219</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>PS_SEP</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ad504ef520e4f60d9d96e055fe7d3fc44</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>PS_QUO</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ab0c64509682c4bbfee0b8e7568c4a3bc</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>iwave_min</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ace72e2401cfeaeb782f6f22affe8c3ac</anchor>
      <arglist>(a, b)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>iwave_max</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>af01a9ca7c3c687ab6b07bd4a63c2b7ae</anchor>
      <arglist>(a, b)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>iwave_abs</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ab19417f88858f924179d618f490c368c</anchor>
      <arglist>(a)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_SUCCESS</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gad1a31c5316d5f51ca5c1c4fb2b9b1958</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_INTERNAL</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga69388fb6308870dea3811ca1b2295be0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_OTHER</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga5104afad49aabbaf0b0be123e00c3458</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_ALLOC</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gad214f6d8c2fa1f4ae1f7bf1aa7d541dd</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_BADINPUT</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gaa10195595b9896967a77a0d6484ccf9a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_OUTOFBOUNDS</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga8660e9fe6627c3ac74e70f5c389858fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_BADINDEX</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gaae579bde80c55674b925f032f7af1524</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_BADARRINDEX</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gac0a11f4edead40c38600a6c241319198</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_BADDIMINDEX</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga0f98e405f9337a15626f8398994db43d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_FILE</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga7292365000032a045f36c05aa99d378f</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_FILEOPEN</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gaf4eb543352a08c2c3007bf530abc3ddc</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_MPI</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga2deab3bf22bdea474e79d10338c6bd65</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_DOMAINDECOMP</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga152113f2216bd54977647d23b690e5a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_PARSE</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga42c38e4efaeac403ac11bb3deecd8695</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_PARSENONAME</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga9b8d6bcc38f3c50a34014b7798edc804</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_PARSENOVALUE</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga6fdbafc843a71402469e7aa9212a1b9f</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_PARSECONVERT</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga50848035574ef8a19e4365b7d629593b</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_ALREADYALLOC</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gab1c00b7c0f3e787a3566f77167a83d4a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_RANGE</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga47953fce4e19bebbc189f365b76f41cd</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_OVERFLOW</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga2ae9bebfa7ace670d4d2718bd17142b1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_UNDERFLOW</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gacb13f1b2a8324564117de8667601d3a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_NOTIMESTEP</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga7d4b8f3b3865e7c856ca35c9f1d17b15</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_NOTINGRID</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gaee96528a048b56ec975946454661d688</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>SEAMX_BIG_ENDIAN</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a5659a41c483b011dadda086d2be97ce0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>SEAMX_LITTLE_ENDIAN</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>afc980d10564dbac749aed6aee02c22ea</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>int</type>
      <name>IPNT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a7cf9e7a7e450c33146f0bd8c0f9405bb</anchor>
      <arglist>[IARR_MAX_NDIM]</arglist>
    </member>
    <member kind="typedef">
      <type>ireal</type>
      <name>RPNT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>af9855b8dba39f31b3fdc31315a05ebd9</anchor>
      <arglist>[RARR_MAX_NDIM]</arglist>
      <docanchor file="utils_8h">RPNT</docanchor>
    </member>
    <member kind="typedef">
      <type>struct s_SIZEDSTRING</type>
      <name>SIZEDSTRING</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a28d56ae3665e5e758a555a594f6eeaec</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int *</type>
      <name>IASN</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a572d076f1ae0dbfb23b5824ffc5a2b68</anchor>
      <arglist>(IPNT l, const IPNT r)</arglist>
    </member>
    <member kind="function">
      <type>ireal *</type>
      <name>RASN</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ab96a2d28ca3d4698e8c6f02bee8446de</anchor>
      <arglist>(RPNT l, const RPNT r)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>usermalloc_</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a8c8fbe5cf00374dbb8f2d54d3a4b6f55</anchor>
      <arglist>(size_t size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>userfree_</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ac985011cf09be6701748cead22e090e4</anchor>
      <arglist>(void *ptr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>gen_3n1</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a2b2a000cb8ab16b4696656f04c552f06</anchor>
      <arglist>(int ndim, int *n)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>gen_i2pnt</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a371f55c71c30527086f403c103b90b08</anchor>
      <arglist>(int ndim, int i, IPNT p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>gen_pnt2i</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a7d219168d905147b2ff5aa91d77a0dc2</anchor>
      <arglist>(int ndim, const IPNT p, int *i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeRank</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ae5e346fa1026e1d0fe78fd9e79861352</anchor>
      <arglist>(int rk)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>retrieveRank</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>aaa511d5c7db3bf809ac4fff260028339</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeSize</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a80d410064e29a402f6bdc52cd142746c</anchor>
      <arglist>(int sz)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>retrieveSize</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a2db1322fe94532c4345ae9a1d32b6a2e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeComm</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>aa9f2b52b79e17e6c24572b5f3ba3650a</anchor>
      <arglist>(MPI_Comm cm)</arglist>
    </member>
    <member kind="function">
      <type>MPI_Comm</type>
      <name>retrieveComm</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a32fdcd2b8369bb359967534e1a9f1c6f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeOutstream</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>adfe332319422f185352372753c708a77</anchor>
      <arglist>(FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>FILE *</type>
      <name>retrieveOutstream</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a7fa12750d26d20922a65a0d312818d3c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeGlobalRank</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a5a51f09b262d5f1ab1c0708562aae60a</anchor>
      <arglist>(int rk)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>retrieveGlobalRank</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ab30eb6dde1b2054f3262eb7b40722513</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeGlobalSize</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ad27854107b29240b6f2385f0a80fff17</anchor>
      <arglist>(int sz)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>retrieveGlobalSize</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ab904c72e89f67edb9f9090cb62c1454b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeGlobalComm</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a471e62cbdfab6f662995a65d08f07629</anchor>
      <arglist>(MPI_Comm cm)</arglist>
    </member>
    <member kind="function">
      <type>MPI_Comm</type>
      <name>retrieveGlobalComm</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>aa2e4850f45a217137fd12c6077f583f9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeThreadSupp</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>af360bd24b88e158f0043d1773922e46d</anchor>
      <arglist>(int ts)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>retrieveThreadSupp</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ab6a99adc0855c2c3e95c3c3cab2ff9df</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeGroupID</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ae90c10295119f9da0725f437e6066e35</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>retrieveGroupID</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a5dbfe72c845baf02fe7cec9a1bf5b8e4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeNumGroups</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a923508bd60518d811ed07265cf3622f1</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>retrieveNumGroups</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a05eaefe4d9e6082e702e591cf81483ef</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>storeRemComm</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a0a67f6157fd7e3e798733d6652e1a3a1</anchor>
      <arglist>(MPI_Comm cm)</arglist>
    </member>
    <member kind="function">
      <type>MPI_Comm</type>
      <name>retrieveRemComm</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a34d49aa92a1a50c9b62a1f60608c3288</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getMachineEndianness</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a79918554ade9505fd2b287193d5757dc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>swapBytes</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a7fd6c5fb3b8e6994aba80590cdab180f</anchor>
      <arglist>(unsigned char *arr, int arrsize, int atomsize)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>fopen</name>
    <title>IWAVE File Management Functions</title>
    <filename>fopen</filename>
  </compound>
  <compound kind="group">
    <name>create</name>
    <title>PARARRAY: Creation and destruction</title>
    <filename>group__create.html</filename>
    <namespace>RVL</namespace>
    <subgroup>print</subgroup>
    <subgroup>ffaccess</subgroup>
    <subgroup>flaccess</subgroup>
    <subgroup>sfassign</subgroup>
    <subgroup>overloaded</subgroup>
    <member kind="function">
      <type>PARARRAY *</type>
      <name>ps_new</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga22610c013d9baddbff98547d73b959f3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ps_delete</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gae7588685e9fcf8ea7c70bc421c896f43</anchor>
      <arglist>(PARARRAY **p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_setnull</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga403dda40e717a2a6fde0ce23f0bd7972</anchor>
      <arglist>(PARARRAY *parr)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gac8ed76b18b430feb617bd2abe41ceda2</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, char &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga3900c5e5488602f104b1e1d4fad95d83</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, short &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga94ae60450fd62a060e9657df7f64d18a</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, int &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga4b052205a4bae64c413a435afdeb68e6</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, long &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga21407734aa4593653ee79ad686f49fe6</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, unsigned short &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga942bd570653f63468d6ad27db005c54f</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, unsigned int &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga408d0b569973fed89d3532c1b71ebcc3</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, unsigned long &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gac72841c30b2e09f7997527e9896f7168</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, float &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gafc15e3471a2b1ab7e13cbcc8947e65b6</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, double &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga9c3dbb432ad3b6448de1858a6f64e39b</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, bool &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>valparse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gaccf70be59e5711495e5933b4e092d95c</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, T def)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>valparse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gaae0f6b6aad214dc4aea9c12af739f833</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name)</arglist>
    </member>
    <member kind="variable">
      <type>WORD *</type>
      <name>val</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga926dddf828443622fab7989506ae5cad</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>struct s_PSLINK *</type>
      <name>prev</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga4bc28d876e048107fdf59494cf7844fd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>struct s_PSLINK *</type>
      <name>next</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gaf8d009aaf7a077d4a55bab2340c58684</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>print</name>
    <title>Output to stream</title>
    <filename>group__print.html</filename>
    <member kind="function">
      <type>int</type>
      <name>ps_printall</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>ga64ce217021050379834f6f38c5987101</anchor>
      <arglist>(PARARRAY parr, FILE *stream)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_createfp</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>ga7aaefb64870ef48bbbaae5326bd39e73</anchor>
      <arglist>(PARARRAY *parr, FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_createfile_fproto</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>ga1a231cf9b18762847f072bf86d4e3dfe</anchor>
      <arglist>(PARARRAY *parr, FILE **stream, const char *proto, const char *filename)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_createfile</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>gafd1f9e4d4112cd2a5a1a358f536747b7</anchor>
      <arglist>(PARARRAY *parr, const char *filename)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_createargs</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>gaed3ad5b588d5854b59fc2ddc4fd43b62</anchor>
      <arglist>(PARARRAY *parr, int argc, char **argv)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_copy</name>
      <anchorfile>group__print.html</anchorfile>
      <anchor>gad44b69c03bec322c9e87f992cb4f4ea9</anchor>
      <arglist>(PARARRAY **tgt, PARARRAY src)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>ffaccess</name>
    <title>Parameter access - first occurence</title>
    <filename>group__ffaccess.html</filename>
    <member kind="function">
      <type>int</type>
      <name>ps_ffcstring</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>gaf826510b46fad583252705525f124672</anchor>
      <arglist>(PARARRAY parr, const char *key, char **p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffchar</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga0ae7a683ffb52fa1e178d0d177e1cd7e</anchor>
      <arglist>(PARARRAY parr, const char *key, char *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffshort</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga0bfd6c5d1c6c1e14fdc6d47982df927a</anchor>
      <arglist>(PARARRAY parr, const char *key, short *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffint</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga4a5876ce56fbb1b103a524c863d7792d</anchor>
      <arglist>(PARARRAY parr, const char *key, int *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_fflong</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>gae7cbab0c40ca721c6bdc6aa5b872747e</anchor>
      <arglist>(PARARRAY parr, const char *key, long *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffushort</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga615269e94a42c65015a0e87ef2054582</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned short *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffuint</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga1d1b56b5c891de96dd19ea5ba9b3c14d</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned int *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffulong</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga5cc39d8804c55f4ef86e4bca41c314ac</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned long *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_fffloat</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>gacbc0af2c3b9aa11c95665344fb1ac0e7</anchor>
      <arglist>(PARARRAY parr, const char *key, float *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffdouble</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga17e18df3535d2e580af71e6693510560</anchor>
      <arglist>(PARARRAY parr, const char *key, double *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_ffreal</name>
      <anchorfile>group__ffaccess.html</anchorfile>
      <anchor>ga7e403f404da8e23ec947139f83113aac</anchor>
      <arglist>(PARARRAY parr, const char *key, ireal *p)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>flaccess</name>
    <title>Parameter access - last occurence</title>
    <filename>group__flaccess.html</filename>
    <member kind="function">
      <type>int</type>
      <name>ps_flcstring</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>gad600512aa48c23bca168afbd56360759</anchor>
      <arglist>(PARARRAY parr, const char *key, char **p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flchar</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga5e96dc7ee7d414859d6d91c925e2de00</anchor>
      <arglist>(PARARRAY parr, const char *key, char *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flshort</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>gabf116354c182aa9b1afd22208eba244d</anchor>
      <arglist>(PARARRAY parr, const char *key, short *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flint</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga5c39554369dfd58ba1c78f7511a3dea8</anchor>
      <arglist>(PARARRAY parr, const char *key, int *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_fllong</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga666f5e165f499ef37b9af543998dcf16</anchor>
      <arglist>(PARARRAY parr, const char *key, long *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flushort</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga23e235903067cc52b6b02f43cd8fb70a</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned short *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_fluint</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga57e5e102e3cc0687b55c0214f676ec1c</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned int *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flulong</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga07df22608785d6d6fc961db32ed34bbd</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned long *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flfloat</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>gaf47a79ee9eeb11f468154c0ff2c3a0fb</anchor>
      <arglist>(PARARRAY parr, const char *key, float *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_fldouble</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>ga5c692ae40d51f5c1300e05e0b9357a52</anchor>
      <arglist>(PARARRAY parr, const char *key, double *p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_flreal</name>
      <anchorfile>group__flaccess.html</anchorfile>
      <anchor>gac460bbe480477a492b8f49eb48476e6c</anchor>
      <arglist>(PARARRAY parr, const char *key, ireal *p)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>sfassign</name>
    <title>Parameter assignment - first occurence</title>
    <filename>group__sfassign.html</filename>
    <member kind="function">
      <type>int</type>
      <name>ps_sfcstring</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga8e8a1042527967504e4fca64f71d9bef</anchor>
      <arglist>(PARARRAY parr, const char *key, const char *val)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfchar</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gae97418a9c1f5515024d8e29ea39c6424</anchor>
      <arglist>(PARARRAY parr, const char *key, char p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfshort</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga74571071be0305000f3ce6a43c1101e1</anchor>
      <arglist>(PARARRAY parr, const char *key, short p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfint</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga30d103654f68511eca81ced1ffbf0e69</anchor>
      <arglist>(PARARRAY parr, const char *key, int p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sflong</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gab06af338eeaffc655e817a01dfec951f</anchor>
      <arglist>(PARARRAY parr, const char *key, long p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfushort</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga4278e02d6089c43797d40fd4ba1e1a9d</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned short p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfuint</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga3f6da9922f01b68cfac2a3f398c857bf</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned int p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfulong</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga4fd84b864795ef403a8575cfc0ac3a15</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned long p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sffloat</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gaf4f62c7b33f05bd632b6252058a3352a</anchor>
      <arglist>(PARARRAY parr, const char *key, float p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfdouble</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga04b714f7fd7aeaab2afb1c5f45978cdc</anchor>
      <arglist>(PARARRAY parr, const char *key, double p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sfreal</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga105a6eb9a20b1847f5cfcf56edc9dc7a</anchor>
      <arglist>(PARARRAY parr, const char *key, ireal p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slcstring</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga47ed5c5bb39d5579e2c0fe47ef147dc9</anchor>
      <arglist>(PARARRAY parr, const char *key, const char *val)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slchar</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga833aab3c778e10765b133974a7470102</anchor>
      <arglist>(PARARRAY parr, const char *key, char p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slshort</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga66a745782546f06eadb3b627a7e706e5</anchor>
      <arglist>(PARARRAY parr, const char *key, short p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slint</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gad9d8015cab8bc775997486106645a184</anchor>
      <arglist>(PARARRAY parr, const char *key, int p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sllong</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga0e8695938688e7c8ca9b043c6491add9</anchor>
      <arglist>(PARARRAY parr, const char *key, long p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slushort</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gab16848ce5d857bfc249685abab98490d</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned short p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sluint</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga37784810b2472c56419b515d51dcb766</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned int p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slulong</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gac279ebee6e5509e3687a15bc36050486</anchor>
      <arglist>(PARARRAY parr, const char *key, unsigned long p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slfloat</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gac797f58997f7b291c34702fd39cbc9d8</anchor>
      <arglist>(PARARRAY parr, const char *key, float p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_sldouble</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>gabfe9137bdf2fafb0e68c3cb6ebda6673</anchor>
      <arglist>(PARARRAY parr, const char *key, double p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ps_slreal</name>
      <anchorfile>group__sfassign.html</anchorfile>
      <anchor>ga30aa66866932e4e69cfbb35f0feed844</anchor>
      <arglist>(PARARRAY parr, const char *key, ireal p)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>overloaded</name>
    <title>parse functions</title>
    <filename>group__overloaded.html</filename>
  </compound>
  <compound kind="group">
    <name>error</name>
    <title>Error codes</title>
    <filename>group__error.html</filename>
    <member kind="define">
      <type>#define</type>
      <name>E_SUCCESS</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gad1a31c5316d5f51ca5c1c4fb2b9b1958</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_INTERNAL</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga69388fb6308870dea3811ca1b2295be0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_OTHER</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga5104afad49aabbaf0b0be123e00c3458</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_ALLOC</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gad214f6d8c2fa1f4ae1f7bf1aa7d541dd</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_BADINPUT</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gaa10195595b9896967a77a0d6484ccf9a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_OUTOFBOUNDS</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga8660e9fe6627c3ac74e70f5c389858fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_BADINDEX</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gaae579bde80c55674b925f032f7af1524</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_BADARRINDEX</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gac0a11f4edead40c38600a6c241319198</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_BADDIMINDEX</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga0f98e405f9337a15626f8398994db43d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_FILE</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga7292365000032a045f36c05aa99d378f</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_FILEOPEN</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gaf4eb543352a08c2c3007bf530abc3ddc</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_MPI</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga2deab3bf22bdea474e79d10338c6bd65</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_DOMAINDECOMP</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga152113f2216bd54977647d23b690e5a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_PARSE</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga42c38e4efaeac403ac11bb3deecd8695</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_PARSENONAME</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga9b8d6bcc38f3c50a34014b7798edc804</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_PARSENOVALUE</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga6fdbafc843a71402469e7aa9212a1b9f</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_PARSECONVERT</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga50848035574ef8a19e4365b7d629593b</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_ALREADYALLOC</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gab1c00b7c0f3e787a3566f77167a83d4a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_RANGE</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga47953fce4e19bebbc189f365b76f41cd</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_OVERFLOW</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga2ae9bebfa7ace670d4d2718bd17142b1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_UNDERFLOW</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gacb13f1b2a8324564117de8667601d3a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_NOTIMESTEP</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>ga7d4b8f3b3865e7c856ca35c9f1d17b15</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>E_NOTINGRID</name>
      <anchorfile>group__error.html</anchorfile>
      <anchor>gaee96528a048b56ec975946454661d688</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>s_KEYVAL</name>
    <filename>structs__KEYVAL.html</filename>
    <member kind="variable">
      <type>WORD *</type>
      <name>key</name>
      <anchorfile>structs__KEYVAL.html</anchorfile>
      <anchor>acbcf6e9826be25b26bdcff3653dafbda</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>WORD *</type>
      <name>val</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga926dddf828443622fab7989506ae5cad</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>s_PARARRAY</name>
    <filename>structs__PARARRAY.html</filename>
    <member kind="variable">
      <type>PSLINK *</type>
      <name>list</name>
      <anchorfile>structs__PARARRAY.html</anchorfile>
      <anchor>a524aca5186a345b2b66bd5c0f41af251</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>s_PSLINK</name>
    <filename>structs__PSLINK.html</filename>
    <member kind="variable">
      <type>KEYVAL *</type>
      <name>pair</name>
      <anchorfile>structs__PSLINK.html</anchorfile>
      <anchor>aa6a1a221a17213a061740e162be3374d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>struct s_PSLINK *</type>
      <name>prev</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga4bc28d876e048107fdf59494cf7844fd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>struct s_PSLINK *</type>
      <name>next</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gaf8d009aaf7a077d4a55bab2340c58684</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>s_SIZEDSTRING</name>
    <filename>structs__SIZEDSTRING.html</filename>
    <member kind="variable">
      <type>long</type>
      <name>n</name>
      <anchorfile>structs__SIZEDSTRING.html</anchorfile>
      <anchor>a8f352f24a6c04a1c3871c66dbe4a4ffe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>s</name>
      <anchorfile>structs__SIZEDSTRING.html</anchorfile>
      <anchor>a3308df0dacfa9085bdf782b6155aee7a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>s_WORD</name>
    <filename>structs__WORD.html</filename>
    <member kind="variable">
      <type>Users williamsymes Applications RSFSRC trip iwave base include parser h char *</type>
      <name>str</name>
      <anchorfile>structs__WORD.html</anchorfile>
      <anchor>adf453c0c6e4814fe3cfeb976662eb1d2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>RVL</name>
    <filename>namespaceRVL.html</filename>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a3c0789d4f7c051d3ce0ebcb99f0e9143</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, string &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gac8ed76b18b430feb617bd2abe41ceda2</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, char &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga3900c5e5488602f104b1e1d4fad95d83</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, short &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga94ae60450fd62a060e9657df7f64d18a</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, int &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga4b052205a4bae64c413a435afdeb68e6</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, long &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga21407734aa4593653ee79ad686f49fe6</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, unsigned short &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga942bd570653f63468d6ad27db005c54f</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, unsigned int &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga408d0b569973fed89d3532c1b71ebcc3</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, unsigned long &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gac72841c30b2e09f7997527e9896f7168</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, float &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gafc15e3471a2b1ab7e13cbcc8947e65b6</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, double &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>parse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>ga9c3dbb432ad3b6448de1858a6f64e39b</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, bool &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>valparse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gaccf70be59e5711495e5933b4e092d95c</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name, T def)</arglist>
    </member>
    <member kind="function">
      <type>T</type>
      <name>valparse</name>
      <anchorfile>group__create.html</anchorfile>
      <anchor>gaae0f6b6aad214dc4aea9c12af739f833</anchor>
      <arglist>(PARARRAY const &amp;par, std::string name)</arglist>
    </member>
  </compound>
</tagfile>
