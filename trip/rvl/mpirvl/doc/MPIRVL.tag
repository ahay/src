<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="page">
    <name>index</name>
    <title>MPIRVL - MPI wrappers for RVL data management classes</title>
    <filename>index</filename>
  </compound>
  <compound kind="file">
    <name>doc.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/mpirvl/include/</path>
    <filename>doc_8h</filename>
  </compound>
  <compound kind="file">
    <name>mpidatatransfer.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/mpirvl/include/</path>
    <filename>mpidatatransfer_8hh</filename>
    <includes id="mpiutils_8hh" name="mpiutils.hh" local="yes" imported="no">mpiutils.hh</includes>
    <class kind="class">RVL::MPI_Sender</class>
    <class kind="class">RVL::MPI_Receiver</class>
    <class kind="class">RVL::MPI_Broadcaster</class>
    <class kind="class">RVL::MPI_Reducer</class>
    <member kind="function">
      <type>MPI_Datatype</type>
      <name>findMPIDatatype&lt; char &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aad3cf839539bc34c95a11db44354f52d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MPI_Datatype</type>
      <name>findMPIDatatype&lt; short &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a0b79f0de0a773c3b1859f9c5b2960b86</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MPI_Datatype</type>
      <name>findMPIDatatype&lt; int &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>abc7ae4f25d1e925d4661cbe3635b0665</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MPI_Datatype</type>
      <name>findMPIDatatype&lt; long &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a8e8a635b4466d1d46c59eeaddd0a396f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MPI_Datatype</type>
      <name>findMPIDatatype&lt; unsigned char &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a4181a752435916eb4223639bebcd691a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MPI_Datatype</type>
      <name>findMPIDatatype&lt; unsigned short &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a842c2bf21200a71ec4f8bdb50b1a7e71</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MPI_Datatype</type>
      <name>findMPIDatatype&lt; unsigned int &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a1894af11e0e1bc43a2f4968663a82ce2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MPI_Datatype</type>
      <name>findMPIDatatype&lt; unsigned long &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a63e586a18a9ceb8d954a19fa4ebcd88a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MPI_Datatype</type>
      <name>findMPIDatatype&lt; float &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a5dea8eed832b7cb47358a7e0395b039f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MPI_Datatype</type>
      <name>findMPIDatatype&lt; double &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a5d0030f6824679f2287e436e8b1ce3ab</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MPI_Datatype</type>
      <name>findMPIDatatype&lt; long double &gt;</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aa16f3a417ce339bfaca9bf5f82a8b0f4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Build_MPI_Datatype</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a4f5cb9324273932895c13612755237d9</anchor>
      <arglist>(MPI_Datatype *mpitype)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>mpiserialdc.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/mpirvl/include/</path>
    <filename>mpiserialdc_8hh</filename>
    <includes id="mpiserialfo_8hh" name="mpiserialfo.hh" local="yes" imported="no">mpiserialfo.hh</includes>
    <class kind="class">RVL::MPISerialDC</class>
    <class kind="class">RVL::MPISerialDCF</class>
    <member kind="define">
      <type>#define</type>
      <name>FRUITCAKE</name>
      <anchorfile>mpiserialdc_8hh.html</anchorfile>
      <anchor>a85ca02f947edfbad14f2087de0ab4fd9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>mpiserialfo.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/mpirvl/include/</path>
    <filename>mpiserialfo_8hh</filename>
    <class kind="class">RVL::MPISynchRoot</class>
    <class kind="class">RVL::MPISerialFunctionObject</class>
    <class kind="class">RVL::MPISerialFunctionObjectRedn</class>
  </compound>
  <compound kind="file">
    <name>mpiseriallap.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/mpirvl/include/</path>
    <filename>mpiseriallap_8hh</filename>
    <includes id="mpiserialfo_8hh" name="mpiserialfo.hh" local="yes" imported="no">mpiserialfo.hh</includes>
    <class kind="class">RVL::MPISerialLAP</class>
  </compound>
  <compound kind="file">
    <name>mpiutils.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/mpirvl/include/</path>
    <filename>mpiutils_8hh</filename>
  </compound>
  <compound kind="file">
    <name>rkstream.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/mpirvl/include/</path>
    <filename>rkstream_8hh</filename>
    <member kind="function">
      <type>void</type>
      <name>makeRankStream</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a9e906cab4cc825880268be5b1a105aa3</anchor>
      <arglist>(ofstream &amp;outfile, int rk, string prefix=&quot;&quot;)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::MPI_Sender</name>
    <filename>classRVL_1_1MPI__Sender.html</filename>
    <templarg></templarg>
    <member kind="function">
      <type></type>
      <name>MPI_Sender</name>
      <anchorfile>classRVL_1_1MPI__Sender.html</anchorfile>
      <anchor>a2959d87cb5810744801ac01dd8fe8566</anchor>
      <arglist>(int _dest=0, MPI_Comm _comm=MPI_COMM_WORLD)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MPI_Sender</name>
      <anchorfile>classRVL_1_1MPI__Sender.html</anchorfile>
      <anchor>a4bbb5de213df8c74d7909804887d4e6d</anchor>
      <arglist>(MPI_Sender&lt; T &gt; const &amp;s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~MPI_Sender</name>
      <anchorfile>classRVL_1_1MPI__Sender.html</anchorfile>
      <anchor>a7c52d05a5ffe87eac6847b51831569f0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>setDestination</name>
      <anchorfile>classRVL_1_1MPI__Sender.html</anchorfile>
      <anchor>a4382f85cb05f5be40d72da09d2832a81</anchor>
      <arglist>(int d)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1MPI__Sender.html</anchorfile>
      <anchor>a165bb415d9c0023421854d2351dc2ba6</anchor>
      <arglist>(T &amp;x) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::MPI_Receiver</name>
    <filename>classRVL_1_1MPI__Receiver.html</filename>
    <templarg></templarg>
    <member kind="function">
      <type></type>
      <name>MPI_Receiver</name>
      <anchorfile>classRVL_1_1MPI__Receiver.html</anchorfile>
      <anchor>a4301aae3615a8158d3fc717b8e9ed912</anchor>
      <arglist>(MPI_Status *_status, int _src=0, MPI_Comm _comm=MPI_COMM_WORLD)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MPI_Receiver</name>
      <anchorfile>classRVL_1_1MPI__Receiver.html</anchorfile>
      <anchor>afbd5b578450e0f1ff75740b54127f5eb</anchor>
      <arglist>(MPI_Receiver&lt; T &gt; const &amp;s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~MPI_Receiver</name>
      <anchorfile>classRVL_1_1MPI__Receiver.html</anchorfile>
      <anchor>aeefe30a25ed0d75b245e067a0773f897</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>setSource</name>
      <anchorfile>classRVL_1_1MPI__Receiver.html</anchorfile>
      <anchor>a973156512d7f410fbc7a3393b2f5b749</anchor>
      <arglist>(int s)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1MPI__Receiver.html</anchorfile>
      <anchor>a57ca5cb65455d949d254cfcd8580edd1</anchor>
      <arglist>(T &amp;x) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::MPI_Broadcaster</name>
    <filename>classRVL_1_1MPI__Broadcaster.html</filename>
    <templarg>T</templarg>
    <member kind="function">
      <type></type>
      <name>MPI_Broadcaster</name>
      <anchorfile>classRVL_1_1MPI__Broadcaster.html</anchorfile>
      <anchor>aaa05150dd94e20029d74d55fc0f62d39</anchor>
      <arglist>(int _root=0, MPI_Comm _comm=MPI_COMM_WORLD)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MPI_Broadcaster</name>
      <anchorfile>classRVL_1_1MPI__Broadcaster.html</anchorfile>
      <anchor>a8b80fb6ac827fbd79c024eb6e3c4435e</anchor>
      <arglist>(MPI_Broadcaster&lt; T &gt; const &amp;b)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~MPI_Broadcaster</name>
      <anchorfile>classRVL_1_1MPI__Broadcaster.html</anchorfile>
      <anchor>a4c0b4da893546202c3e4d09389bf4c43</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1MPI__Broadcaster.html</anchorfile>
      <anchor>ae07e98734d09acc7ce22b252b9563878</anchor>
      <arglist>(T &amp;x) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::MPI_Reducer</name>
    <filename>classRVL_1_1MPI__Reducer.html</filename>
    <templarg></templarg>
    <member kind="function">
      <type></type>
      <name>MPI_Reducer</name>
      <anchorfile>classRVL_1_1MPI__Reducer.html</anchorfile>
      <anchor>a0d9df657dedec45cf1d6917c6565801d</anchor>
      <arglist>(MPI_Op _op=MPI_SUM, int _root=0, MPI_Comm _comm=MPI_COMM_WORLD)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MPI_Reducer</name>
      <anchorfile>classRVL_1_1MPI__Reducer.html</anchorfile>
      <anchor>a91a79e577a67c4f83da7ce7b66ef7ee0</anchor>
      <arglist>(MPI_Reducer&lt; T &gt; const &amp;b)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~MPI_Reducer</name>
      <anchorfile>classRVL_1_1MPI__Reducer.html</anchorfile>
      <anchor>ada4ca52e0b4371466bbf96a8a650273a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1MPI__Reducer.html</anchorfile>
      <anchor>a38cf66791973796c04c5612f48a0ac97</anchor>
      <arglist>(T &amp;xout, T const &amp;xin) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::MPISerialDC</name>
    <filename>classRVL_1_1MPISerialDC.html</filename>
    <base>RVL::DataContainer</base>
    <member kind="function">
      <type></type>
      <name>MPISerialDC</name>
      <anchorfile>classRVL_1_1MPISerialDC.html</anchorfile>
      <anchor>a96dcd36a488ba0ce403029b48ec70d26</anchor>
      <arglist>(DataContainerFactory const &amp;F)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~MPISerialDC</name>
      <anchorfile>classRVL_1_1MPISerialDC.html</anchorfile>
      <anchor>ae1c44d123ab32a1c44abe88a7a15b00a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1MPISerialDC.html</anchorfile>
      <anchor>ad63e6a5042cbb115909dedbb103885a8</anchor>
      <arglist>(FunctionObject &amp;f, std::vector&lt; DataContainer const * &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1MPISerialDC.html</anchorfile>
      <anchor>a74a28c441323822257154a84a2fc28a5</anchor>
      <arglist>(FunctionObjectConstEval &amp;f, vector&lt; DataContainer const * &gt; &amp;x) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1MPISerialDC.html</anchorfile>
      <anchor>ac0102081f6017bc163b71c078fff7613</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>rk</name>
      <anchorfile>classRVL_1_1MPISerialDC.html</anchorfile>
      <anchor>a94fffedc11d5bdf8e977a3cbf5f3ffbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>DataContainer *</type>
      <name>mydc</name>
      <anchorfile>classRVL_1_1MPISerialDC.html</anchorfile>
      <anchor>a763bbc7130fcd3b0edb03a49090aa24c</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>MPISerialDCF</name>
      <anchorfile>classRVL_1_1MPISerialDC.html</anchorfile>
      <anchor>ac44282e6d289495ca72d73293b61f579</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::MPISerialDCF</name>
    <filename>classRVL_1_1MPISerialDCF.html</filename>
    <base>RVL::DataContainerFactory</base>
    <member kind="function">
      <type></type>
      <name>MPISerialDCF</name>
      <anchorfile>classRVL_1_1MPISerialDCF.html</anchorfile>
      <anchor>af0dd003e3011daf4eed1d4e891eaa5ef</anchor>
      <arglist>(DataContainerFactory const &amp;_f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MPISerialDCF</name>
      <anchorfile>classRVL_1_1MPISerialDCF.html</anchorfile>
      <anchor>a785bb57c954446d50f6d418c6639d3d0</anchor>
      <arglist>(MPISerialDCF const &amp;fact)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~MPISerialDCF</name>
      <anchorfile>classRVL_1_1MPISerialDCF.html</anchorfile>
      <anchor>ab8d06b85322a2301f4e457cd37350a74</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DataContainer *</type>
      <name>build</name>
      <anchorfile>classRVL_1_1MPISerialDCF.html</anchorfile>
      <anchor>a29b962ed6547aed032dceff1d9aa59c2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>compare</name>
      <anchorfile>classRVL_1_1MPISerialDCF.html</anchorfile>
      <anchor>a9ee0802b9a415a6f44af69c8b41bcc29</anchor>
      <arglist>(DataContainerFactory const &amp;dcf) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1MPISerialDCF.html</anchorfile>
      <anchor>a6fcfca7f70351511eeb74c47aa259cc9</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1MPISerialDCF.html</anchorfile>
      <anchor>af90b696358a9ef5b991b36500d60f64c</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::MPISynchRoot</name>
    <filename>classRVL_1_1MPISynchRoot.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>set</name>
      <anchorfile>classRVL_1_1MPISynchRoot.html</anchorfile>
      <anchor>a1f3a9d5c40f6d923e89a92fa5a125c46</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>synch</name>
      <anchorfile>classRVL_1_1MPISynchRoot.html</anchorfile>
      <anchor>ae1939b302680835809111b8b4d7a9fe1</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~MPISynchRoot</name>
      <anchorfile>classRVL_1_1MPISynchRoot.html</anchorfile>
      <anchor>a35ccc8a0baf2e1e99b70bba6542a3961</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::MPISerialFunctionObject</name>
    <filename>classRVL_1_1MPISerialFunctionObject.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::LocalFunctionObject</base>
    <member kind="function">
      <type></type>
      <name>MPISerialFunctionObject</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObject.html</anchorfile>
      <anchor>a33daf2a7d62fdf948df2d4373312561b</anchor>
      <arglist>(FunctionObject &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MPISerialFunctionObject</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObject.html</anchorfile>
      <anchor>a8c58fdddd44b777ac3f8ec321c698ca0</anchor>
      <arglist>(MPISerialFunctionObject&lt; DataType &gt; const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~MPISerialFunctionObject</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObject.html</anchorfile>
      <anchor>a7a92567ad8c17e6dc4df6263cbd49d45</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObject.html</anchorfile>
      <anchor>a4a51547973a45cdfeefba444cd67d549</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;target, vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObject.html</anchorfile>
      <anchor>aac7586e7821d0049b05db6969abf9b94</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::MPISerialFunctionObjectRedn</name>
    <filename>classRVL_1_1MPISerialFunctionObjectRedn.html</filename>
    <templarg>DataType</templarg>
    <templarg>ValType</templarg>
    <base>RVL::FunctionObjectScalarRedn</base>
    <base>RVL::LocalConstEval</base>
    <base>RVL::MPISynchRoot</base>
    <member kind="function">
      <type></type>
      <name>MPISerialFunctionObjectRedn</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObjectRedn.html</anchorfile>
      <anchor>ab4a0b9cd9278509d184c1afb09da6ef1</anchor>
      <arglist>(FunctionObjectScalarRedn&lt; ValType &gt; &amp;_f, int _root=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MPISerialFunctionObjectRedn</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObjectRedn.html</anchorfile>
      <anchor>aa0bb0dccb12721bdab33317ea5361ca0</anchor>
      <arglist>(MPISerialFunctionObjectRedn&lt; DataType, ValType &gt; const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~MPISerialFunctionObjectRedn</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObjectRedn.html</anchorfile>
      <anchor>a5406a2e267b4567e524667e69b1e6254</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObjectRedn.html</anchorfile>
      <anchor>a94873c2e946f5ed16ace434f67ecabf0</anchor>
      <arglist>(vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>synch</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObjectRedn.html</anchorfile>
      <anchor>aa2a61cd225de47b611d151ae52f2dea5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setValue</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObjectRedn.html</anchorfile>
      <anchor>a534b0a43802538626d8cc645f637b1c6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObjectRedn.html</anchorfile>
      <anchor>aa0fcdb787de923a4a83c4c80d99033df</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1MPISerialFunctionObjectRedn.html</anchorfile>
      <anchor>a897fd7040ac97bb87d2d75e925de9586</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::MPISerialLAP</name>
    <filename>classRVL_1_1MPISerialLAP.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearAlgebraPackage</base>
    <member kind="function">
      <type></type>
      <name>MPISerialLAP</name>
      <anchorfile>classRVL_1_1MPISerialLAP.html</anchorfile>
      <anchor>a3a3b32dcb5047c3c7f717d29121ddaa7</anchor>
      <arglist>(Scalar ipscale=ScalarFieldTraits&lt; Scalar &gt;::One(), int _root=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MPISerialLAP</name>
      <anchorfile>classRVL_1_1MPISerialLAP.html</anchorfile>
      <anchor>a9e39207ecf5322e32a7fc2f23b564530</anchor>
      <arglist>(const MPISerialLAP&lt; Scalar &gt; &amp;p)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~MPISerialLAP</name>
      <anchorfile>classRVL_1_1MPISerialLAP.html</anchorfile>
      <anchor>ae5d0ab260fd6ac0a2f29f46eb8399963</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>FunctionObjectScalarRedn&lt; Scalar &gt; &amp;</type>
      <name>inner</name>
      <anchorfile>classRVL_1_1MPISerialLAP.html</anchorfile>
      <anchor>a6082cad2dd73928443622c5cebc6b0d8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>FunctionObject &amp;</type>
      <name>zero</name>
      <anchorfile>classRVL_1_1MPISerialLAP.html</anchorfile>
      <anchor>ad3a6ab0b74ad51ddb7d71dff28489ca7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>LinCombObject&lt; Scalar &gt; &amp;</type>
      <name>linComb</name>
      <anchorfile>classRVL_1_1MPISerialLAP.html</anchorfile>
      <anchor>a391e2e77bc18b6d1354f3e30f43e0773</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>compare</name>
      <anchorfile>classRVL_1_1MPISerialLAP.html</anchorfile>
      <anchor>ad0914da8f87158593f31e28d73eff103</anchor>
      <arglist>(LinearAlgebraPackage&lt; Scalar &gt; const &amp;lap) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setScale</name>
      <anchorfile>classRVL_1_1MPISerialLAP.html</anchorfile>
      <anchor>a908db84fa08b1c5393499539410d16e3</anchor>
      <arglist>(Scalar newscale)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1MPISerialLAP.html</anchorfile>
      <anchor>a21016bcb92d18c6b01dfd5c2b6982e66</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
</tagfile>
