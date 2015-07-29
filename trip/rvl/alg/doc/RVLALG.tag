<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="page">
    <name>index</name>
    <title>Algorithm</title>
    <filename>index</filename>
  </compound>
  <compound kind="file">
    <name>alg.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/alg/include/</path>
    <filename>alg_8hh</filename>
    <class kind="class">RVLAlg::Algorithm</class>
    <class kind="class">RVLAlg::NoAlg</class>
    <class kind="class">RVLAlg::ListAlg</class>
    <class kind="class">RVLAlg::Terminator</class>
    <class kind="class">RVLAlg::LoopAlg</class>
    <class kind="class">RVLAlg::DoLoopAlg</class>
    <class kind="class">RVLAlg::CondListAlg</class>
    <class kind="class">RVLAlg::StateAlg</class>
    <class kind="class">RVLAlg::BranchAlg</class>
    <namespace>RVLAlg</namespace>
    <member kind="define">
      <type>#define</type>
      <name>STOP_LOOP</name>
      <anchorfile>alg_8hh.html</anchorfile>
      <anchor>aa7b70bfc48fefea6fd33dbe037d18240</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>CONTINUE_LOOP</name>
      <anchorfile>alg_8hh.html</anchorfile>
      <anchor>a99e1bb98060f61f8ca7b470e48df20c7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>boolterm.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/alg/include/</path>
    <filename>boolterm_8hh</filename>
    <includes id="alg_8hh" name="alg.hh" local="yes" imported="no">alg.hh</includes>
    <class kind="class">RVLAlg::BoolTerminator</class>
    <class kind="class">RVLAlg::AndTerminator</class>
    <class kind="class">RVLAlg::OrTerminator</class>
    <class kind="class">RVLAlg::NotTerminator</class>
    <class kind="class">RVLAlg::XorTerminator</class>
    <namespace>RVLAlg</namespace>
  </compound>
  <compound kind="file">
    <name>doc.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/alg/include/</path>
    <filename>doc_8h</filename>
  </compound>
  <compound kind="file">
    <name>ioterm.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/alg/include/</path>
    <filename>ioterm_8hh</filename>
    <includes id="alg_8hh" name="alg.hh" local="yes" imported="no">alg.hh</includes>
    <includes id="scalarterm_8hh" name="scalarterm.hh" local="yes" imported="no">scalarterm.hh</includes>
    <class kind="class">RVLAlg::IOTerminator</class>
    <class kind="class">RVLAlg::VecWatchTerminator</class>
    <class kind="class">RVLAlg::IterationTable</class>
    <class kind="class">RVLAlg::SteppedIterationTable</class>
    <class kind="class">RVLAlg::GradientThresholdIterationTable</class>
    <class kind="class">RVLAlg::CountingThresholdIterationTable</class>
    <class kind="class">RVLAlg::VectorCountingThresholdIterationTable</class>
    <class kind="class">RVLAlg::CountingNormTable</class>
    <namespace>RVLAlg</namespace>
  </compound>
  <compound kind="file">
    <name>scalarterm.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/alg/include/</path>
    <filename>scalarterm_8hh</filename>
    <includes id="alg_8hh" name="alg.hh" local="yes" imported="no">alg.hh</includes>
    <class kind="class">RVLAlg::CountTerminator</class>
    <class kind="class">RVLAlg::MaxTerminator</class>
    <class kind="class">RVLAlg::MinTerminator</class>
    <class kind="class">RVLAlg::MinTerminatorFE</class>
    <namespace>RVLAlg</namespace>
  </compound>
  <compound kind="file">
    <name>terminator.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/alg/include/</path>
    <filename>terminator_8hh</filename>
    <includes id="boolterm_8hh" name="boolterm.hh" local="yes" imported="no">boolterm.hh</includes>
    <includes id="scalarterm_8hh" name="scalarterm.hh" local="yes" imported="no">scalarterm.hh</includes>
    <includes id="ioterm_8hh" name="ioterm.hh" local="yes" imported="no">ioterm.hh</includes>
    <includes id="vectorterm_8hh" name="vectorterm.hh" local="yes" imported="no">vectorterm.hh</includes>
  </compound>
  <compound kind="file">
    <name>vectorterm.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/alg/include/</path>
    <filename>vectorterm_8hh</filename>
    <includes id="alg_8hh" name="alg.hh" local="yes" imported="no">alg.hh</includes>
    <class kind="class">RVLAlg::UnaryThresholdTerminator</class>
    <class kind="class">RVLAlg::BinaryThresholdTerminator</class>
    <class kind="class">RVLAlg::TernaryThresholdTerminator</class>
    <class kind="class">RVLAlg::NormThresholdTerminator</class>
    <class kind="class">RVLAlg::Norm2ThresholdTerminator</class>
    <class kind="class">RVLAlg::DiffThresholdTerminator</class>
    <class kind="class">RVLAlg::Diff2ThresholdTerminator</class>
    <class kind="class">RVLAlg::IPThresholdTerminator</class>
    <class kind="class">RVLAlg::AbsIPThresholdTerminator</class>
    <class kind="class">RVLAlg::NormGradientTerminator</class>
    <class kind="class">RVLAlg::DiffBallProjTerminator</class>
    <class kind="class">RVLAlg::BallProjTerminator</class>
    <namespace>RVLAlg</namespace>
  </compound>
  <compound kind="namespace">
    <name>RVLAlg</name>
    <filename>namespaceRVLAlg.html</filename>
    <class kind="class">RVLAlg::Algorithm</class>
    <class kind="class">RVLAlg::NoAlg</class>
    <class kind="class">RVLAlg::ListAlg</class>
    <class kind="class">RVLAlg::Terminator</class>
    <class kind="class">RVLAlg::LoopAlg</class>
    <class kind="class">RVLAlg::DoLoopAlg</class>
    <class kind="class">RVLAlg::CondListAlg</class>
    <class kind="class">RVLAlg::StateAlg</class>
    <class kind="class">RVLAlg::BranchAlg</class>
    <class kind="class">RVLAlg::BoolTerminator</class>
    <class kind="class">RVLAlg::AndTerminator</class>
    <class kind="class">RVLAlg::OrTerminator</class>
    <class kind="class">RVLAlg::NotTerminator</class>
    <class kind="class">RVLAlg::XorTerminator</class>
    <class kind="class">RVLAlg::IOTerminator</class>
    <class kind="class">RVLAlg::VecWatchTerminator</class>
    <class kind="class">RVLAlg::IterationTable</class>
    <class kind="class">RVLAlg::SteppedIterationTable</class>
    <class kind="class">RVLAlg::GradientThresholdIterationTable</class>
    <class kind="class">RVLAlg::CountingThresholdIterationTable</class>
    <class kind="class">RVLAlg::VectorCountingThresholdIterationTable</class>
    <class kind="class">RVLAlg::CountingNormTable</class>
    <class kind="class">RVLAlg::CountTerminator</class>
    <class kind="class">RVLAlg::MaxTerminator</class>
    <class kind="class">RVLAlg::MinTerminator</class>
    <class kind="class">RVLAlg::MinTerminatorFE</class>
    <class kind="class">RVLAlg::UnaryThresholdTerminator</class>
    <class kind="class">RVLAlg::BinaryThresholdTerminator</class>
    <class kind="class">RVLAlg::TernaryThresholdTerminator</class>
    <class kind="class">RVLAlg::NormThresholdTerminator</class>
    <class kind="class">RVLAlg::Norm2ThresholdTerminator</class>
    <class kind="class">RVLAlg::DiffThresholdTerminator</class>
    <class kind="class">RVLAlg::Diff2ThresholdTerminator</class>
    <class kind="class">RVLAlg::IPThresholdTerminator</class>
    <class kind="class">RVLAlg::AbsIPThresholdTerminator</class>
    <class kind="class">RVLAlg::NormGradientTerminator</class>
    <class kind="class">RVLAlg::DiffBallProjTerminator</class>
    <class kind="class">RVLAlg::BallProjTerminator</class>
  </compound>
  <compound kind="class">
    <name>RVLAlg::Algorithm</name>
    <filename>classRVLAlg_1_1Algorithm.html</filename>
    <member kind="function">
      <type></type>
      <name>Algorithm</name>
      <anchorfile>classRVLAlg_1_1Algorithm.html</anchorfile>
      <anchor>a6689419a9052041e71b7d3fb88a43699</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Algorithm</name>
      <anchorfile>classRVLAlg_1_1Algorithm.html</anchorfile>
      <anchor>a267aa869644db6d4a94d1a2435c3b46c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>run</name>
      <anchorfile>classRVLAlg_1_1Algorithm.html</anchorfile>
      <anchor>a1f2607e339d136836946b756b7e21dc0</anchor>
      <arglist>()=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::NoAlg</name>
    <filename>classRVLAlg_1_1NoAlg.html</filename>
    <base>RVLAlg::Algorithm</base>
    <member kind="function">
      <type>void</type>
      <name>run</name>
      <anchorfile>classRVLAlg_1_1NoAlg.html</anchorfile>
      <anchor>a7993c23e933763f64d5ce561faea2b1e</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::ListAlg</name>
    <filename>classRVLAlg_1_1ListAlg.html</filename>
    <base>RVLAlg::Algorithm</base>
    <member kind="function">
      <type></type>
      <name>ListAlg</name>
      <anchorfile>classRVLAlg_1_1ListAlg.html</anchorfile>
      <anchor>aaaac913ca17059807e7751bf1c3bbb9c</anchor>
      <arglist>(Algorithm &amp;first)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ListAlg</name>
      <anchorfile>classRVLAlg_1_1ListAlg.html</anchorfile>
      <anchor>a6de6cf315abb5a1262e95642aa8b31fc</anchor>
      <arglist>(Algorithm &amp;first, Algorithm &amp;next)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>run</name>
      <anchorfile>classRVLAlg_1_1ListAlg.html</anchorfile>
      <anchor>ab7beb752038f55afe4583bda5d9e716e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>islist</name>
      <anchorfile>classRVLAlg_1_1ListAlg.html</anchorfile>
      <anchor>a8e694fdc6d2f68fc5abfa73e20bc0cee</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Algorithm &amp;</type>
      <name>one</name>
      <anchorfile>classRVLAlg_1_1ListAlg.html</anchorfile>
      <anchor>ae92a279d1d0ba8a4b58132c56b7e512d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Algorithm &amp;</type>
      <name>two</name>
      <anchorfile>classRVLAlg_1_1ListAlg.html</anchorfile>
      <anchor>a0f4d6b7279981e9452b50a1d7ba51284</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::Terminator</name>
    <filename>classRVLAlg_1_1Terminator.html</filename>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Terminator</name>
      <anchorfile>classRVLAlg_1_1Terminator.html</anchorfile>
      <anchor>a8497767b414c5bc32a25c7d4c3acbc8b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1Terminator.html</anchorfile>
      <anchor>ab7a71a3442d4c8b6d16f3e56e8e7e985</anchor>
      <arglist>()=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::LoopAlg</name>
    <filename>classRVLAlg_1_1LoopAlg.html</filename>
    <base>RVLAlg::Algorithm</base>
    <member kind="function">
      <type></type>
      <name>LoopAlg</name>
      <anchorfile>classRVLAlg_1_1LoopAlg.html</anchorfile>
      <anchor>a0e2204df84a5f819f13fea9846ce6412</anchor>
      <arglist>(Algorithm &amp;alg, Terminator &amp;stop)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>run</name>
      <anchorfile>classRVLAlg_1_1LoopAlg.html</anchorfile>
      <anchor>af2d7f24aa8cd40a47124f262ca1325d7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Algorithm &amp;</type>
      <name>inside</name>
      <anchorfile>classRVLAlg_1_1LoopAlg.html</anchorfile>
      <anchor>a7579b60272ddce577ab84e2b2fbed810</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Terminator &amp;</type>
      <name>term</name>
      <anchorfile>classRVLAlg_1_1LoopAlg.html</anchorfile>
      <anchor>a8248d5a039e70e89cdee47d97f8718f9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::DoLoopAlg</name>
    <filename>classRVLAlg_1_1DoLoopAlg.html</filename>
    <base>RVLAlg::LoopAlg</base>
    <member kind="function">
      <type></type>
      <name>DoLoopAlg</name>
      <anchorfile>classRVLAlg_1_1DoLoopAlg.html</anchorfile>
      <anchor>a65cb6aac08ecd70d9ba3360a65966f49</anchor>
      <arglist>(Algorithm &amp;alg, Terminator &amp;stop)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>run</name>
      <anchorfile>classRVLAlg_1_1DoLoopAlg.html</anchorfile>
      <anchor>a746d125cc0aeb1951dbfd7a4e347aa4d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::CondListAlg</name>
    <filename>classRVLAlg_1_1CondListAlg.html</filename>
    <base>RVLAlg::ListAlg</base>
    <member kind="function">
      <type></type>
      <name>CondListAlg</name>
      <anchorfile>classRVLAlg_1_1CondListAlg.html</anchorfile>
      <anchor>a1b9200bc668c92a3a763c517e8b64ff8</anchor>
      <arglist>(Algorithm &amp;first, Algorithm &amp;next, Terminator &amp;_stop)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>run</name>
      <anchorfile>classRVLAlg_1_1CondListAlg.html</anchorfile>
      <anchor>a38ad26119440a66a626b9432d7b62b93</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Terminator &amp;</type>
      <name>stop</name>
      <anchorfile>classRVLAlg_1_1CondListAlg.html</anchorfile>
      <anchor>a59f66ed134ced73f3f9148de5ed52ac3</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::StateAlg</name>
    <filename>classRVLAlg_1_1StateAlg.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Algorithm</base>
    <member kind="function" virtualness="pure">
      <type>virtual T &amp;</type>
      <name>getState</name>
      <anchorfile>classRVLAlg_1_1StateAlg.html</anchorfile>
      <anchor>a9040cd2624bcc95915d9193dc8484114</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual const T &amp;</type>
      <name>getState</name>
      <anchorfile>classRVLAlg_1_1StateAlg.html</anchorfile>
      <anchor>a9d548cb80b5f7dc4cc5107280d93ab3c</anchor>
      <arglist>() const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::BranchAlg</name>
    <filename>classRVLAlg_1_1BranchAlg.html</filename>
    <base>RVLAlg::Algorithm</base>
    <member kind="function">
      <type></type>
      <name>BranchAlg</name>
      <anchorfile>classRVLAlg_1_1BranchAlg.html</anchorfile>
      <anchor>ada2750504c555077da6a1e3956d738b7</anchor>
      <arglist>(Terminator &amp;iftest, Algorithm &amp;thenclause, Algorithm &amp;elseclause)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>run</name>
      <anchorfile>classRVLAlg_1_1BranchAlg.html</anchorfile>
      <anchor>ab76bf9affcd437792cdbd38caf7c2280</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Algorithm &amp;</type>
      <name>thencl</name>
      <anchorfile>classRVLAlg_1_1BranchAlg.html</anchorfile>
      <anchor>a42377b241ea88a80c97317d0a0f3113c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Algorithm &amp;</type>
      <name>elsecl</name>
      <anchorfile>classRVLAlg_1_1BranchAlg.html</anchorfile>
      <anchor>a61ab944985c7917bb121a3e2de56786a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Terminator &amp;</type>
      <name>test</name>
      <anchorfile>classRVLAlg_1_1BranchAlg.html</anchorfile>
      <anchor>a09078a0cf58135c99e76437ced4d3109</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::BoolTerminator</name>
    <filename>classRVLAlg_1_1BoolTerminator.html</filename>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>BoolTerminator</name>
      <anchorfile>classRVLAlg_1_1BoolTerminator.html</anchorfile>
      <anchor>a267ef467ea62464e83928475e929b83f</anchor>
      <arglist>(bool _ans=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BoolTerminator</name>
      <anchorfile>classRVLAlg_1_1BoolTerminator.html</anchorfile>
      <anchor>a405f34ac6ce4a2ad39f6ef197ff2ebd5</anchor>
      <arglist>(BoolTerminator &amp;bt)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BoolTerminator</name>
      <anchorfile>classRVLAlg_1_1BoolTerminator.html</anchorfile>
      <anchor>a0439b723cf0ef39d2e32091b8b539853</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setValue</name>
      <anchorfile>classRVLAlg_1_1BoolTerminator.html</anchorfile>
      <anchor>a2d9e096f33b5f5f68a8fd5e11c90a181</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setValue</name>
      <anchorfile>classRVLAlg_1_1BoolTerminator.html</anchorfile>
      <anchor>ae241c2daa44de6e6f155eee7e53a983a</anchor>
      <arglist>(bool _ans)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1BoolTerminator.html</anchorfile>
      <anchor>a4b8dbb0c7eb9ee289d527bc600d3f03d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::AndTerminator</name>
    <filename>classRVLAlg_1_1AndTerminator.html</filename>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>AndTerminator</name>
      <anchorfile>classRVLAlg_1_1AndTerminator.html</anchorfile>
      <anchor>abdc6c20e9b3b65509093720eb5dc6d35</anchor>
      <arglist>(Terminator &amp;a, Terminator &amp;b)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AndTerminator</name>
      <anchorfile>classRVLAlg_1_1AndTerminator.html</anchorfile>
      <anchor>a278fabbbb3a96cc3850fda99688fc2c6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1AndTerminator.html</anchorfile>
      <anchor>af0e1fa1b166bb8477b48317969b5ae4b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Terminator &amp;</type>
      <name>first</name>
      <anchorfile>classRVLAlg_1_1AndTerminator.html</anchorfile>
      <anchor>a748abad65cd958f92488253f80d78046</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Terminator &amp;</type>
      <name>second</name>
      <anchorfile>classRVLAlg_1_1AndTerminator.html</anchorfile>
      <anchor>a85d8b159504d75f8a115d747735731cb</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::OrTerminator</name>
    <filename>classRVLAlg_1_1OrTerminator.html</filename>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>OrTerminator</name>
      <anchorfile>classRVLAlg_1_1OrTerminator.html</anchorfile>
      <anchor>a1a4ef156afe05484ca1a4ec403cdd7c5</anchor>
      <arglist>(Terminator &amp;a, Terminator &amp;b)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~OrTerminator</name>
      <anchorfile>classRVLAlg_1_1OrTerminator.html</anchorfile>
      <anchor>acf0b2c9df93f29474c90d1e6ce57ef13</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1OrTerminator.html</anchorfile>
      <anchor>a69112208d26d4d6e14b9bf637354be06</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Terminator &amp;</type>
      <name>first</name>
      <anchorfile>classRVLAlg_1_1OrTerminator.html</anchorfile>
      <anchor>ac1bb7ba44a63a1652f4bbd8d751f3b1b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Terminator &amp;</type>
      <name>second</name>
      <anchorfile>classRVLAlg_1_1OrTerminator.html</anchorfile>
      <anchor>aba35ad2927ffddf7b7ebe4fa09aadc65</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::NotTerminator</name>
    <filename>classRVLAlg_1_1NotTerminator.html</filename>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>NotTerminator</name>
      <anchorfile>classRVLAlg_1_1NotTerminator.html</anchorfile>
      <anchor>adae7353f8ce5a6b9d7d6ae13d1c9ad02</anchor>
      <arglist>(Terminator &amp;a)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~NotTerminator</name>
      <anchorfile>classRVLAlg_1_1NotTerminator.html</anchorfile>
      <anchor>a2b16fbd60b806fafaa293a9d2f766a84</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1NotTerminator.html</anchorfile>
      <anchor>a634554cab289f0a2587cd4ba9eb9f46f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Terminator &amp;</type>
      <name>first</name>
      <anchorfile>classRVLAlg_1_1NotTerminator.html</anchorfile>
      <anchor>a795570d1e55e01fc16764e0f598093aa</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::XorTerminator</name>
    <filename>classRVLAlg_1_1XorTerminator.html</filename>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>XorTerminator</name>
      <anchorfile>classRVLAlg_1_1XorTerminator.html</anchorfile>
      <anchor>ac2421e71821a3595b014703fd8f577a2</anchor>
      <arglist>(Terminator &amp;a, Terminator &amp;b)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~XorTerminator</name>
      <anchorfile>classRVLAlg_1_1XorTerminator.html</anchorfile>
      <anchor>a7d0ffbfca3b445b82ebc1d995829ce84</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1XorTerminator.html</anchorfile>
      <anchor>a53ddb295b6fa9c41796eb7d26cee22ee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Terminator &amp;</type>
      <name>first</name>
      <anchorfile>classRVLAlg_1_1XorTerminator.html</anchorfile>
      <anchor>aedf53c241429c332a64b69801de3de4c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Terminator &amp;</type>
      <name>second</name>
      <anchorfile>classRVLAlg_1_1XorTerminator.html</anchorfile>
      <anchor>ac2b4fc41e5ebf8b0dceb9ac4f2da5153</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::IOTerminator</name>
    <filename>classRVLAlg_1_1IOTerminator.html</filename>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>IOTerminator</name>
      <anchorfile>classRVLAlg_1_1IOTerminator.html</anchorfile>
      <anchor>a73736ea6cc5c270a92d985276aeee61b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IOTerminator</name>
      <anchorfile>classRVLAlg_1_1IOTerminator.html</anchorfile>
      <anchor>aade5cb5b75b3e74cce7b1b9d2b686150</anchor>
      <arglist>(char t[])</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IOTerminator</name>
      <anchorfile>classRVLAlg_1_1IOTerminator.html</anchorfile>
      <anchor>acec4221b7c0381ac4a0df5a694e4d4e9</anchor>
      <arglist>(istream &amp;in_, ostream &amp;out_)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IOTerminator</name>
      <anchorfile>classRVLAlg_1_1IOTerminator.html</anchorfile>
      <anchor>aeae4cc102cb06ee1a21469f93e816d8f</anchor>
      <arglist>(char t[], istream &amp;in_, ostream &amp;out_)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~IOTerminator</name>
      <anchorfile>classRVLAlg_1_1IOTerminator.html</anchorfile>
      <anchor>ab9709339c463838c02251954edba24cd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1IOTerminator.html</anchorfile>
      <anchor>a0e10616fc369c223bf911ee18e4412d9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>char *</type>
      <name>s</name>
      <anchorfile>classRVLAlg_1_1IOTerminator.html</anchorfile>
      <anchor>a392df189403609990346b89a93ad426a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>istream &amp;</type>
      <name>ins</name>
      <anchorfile>classRVLAlg_1_1IOTerminator.html</anchorfile>
      <anchor>a0510c4a55ca74761dd568459e913aff7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ostream &amp;</type>
      <name>outs</name>
      <anchorfile>classRVLAlg_1_1IOTerminator.html</anchorfile>
      <anchor>ab5cad01eb1520e89a6b2808c1895bf84</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::VecWatchTerminator</name>
    <filename>classRVLAlg_1_1VecWatchTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>VecWatchTerminator</name>
      <anchorfile>classRVLAlg_1_1VecWatchTerminator.html</anchorfile>
      <anchor>a61ff42163c581c58d810401172804082</anchor>
      <arglist>(Vector&lt; Scalar &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>VecWatchTerminator</name>
      <anchorfile>classRVLAlg_1_1VecWatchTerminator.html</anchorfile>
      <anchor>a4214209d3d3235db13a7e62bcda80bb5</anchor>
      <arglist>(Vector&lt; Scalar &gt; &amp;x, ostream &amp;outs)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1VecWatchTerminator.html</anchorfile>
      <anchor>a3c8d15d170b9c4c3d622d62cda2c86f0</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::IterationTable</name>
    <filename>classRVLAlg_1_1IterationTable.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>IterationTable</name>
      <anchorfile>classRVLAlg_1_1IterationTable.html</anchorfile>
      <anchor>a58bb0d647561b94ee90cf9794725d2d1</anchor>
      <arglist>(FunctionalEvaluation&lt; Scalar &gt; &amp;_fx, ostream &amp;_out=cout)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1IterationTable.html</anchorfile>
      <anchor>a18313767c2cbd03cf15c6ea6ace8ca7e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>count</name>
      <anchorfile>classRVLAlg_1_1IterationTable.html</anchorfile>
      <anchor>a66231fb48fc397f6ab1208f567452bcf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>FunctionalEvaluation&lt; Scalar &gt; &amp;</type>
      <name>fx</name>
      <anchorfile>classRVLAlg_1_1IterationTable.html</anchorfile>
      <anchor>aeee11366845efdd6ac55a62f3fa25c54</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ostream &amp;</type>
      <name>out</name>
      <anchorfile>classRVLAlg_1_1IterationTable.html</anchorfile>
      <anchor>a903655933d22548b0fc708bd08a92548</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::SteppedIterationTable</name>
    <filename>classRVLAlg_1_1SteppedIterationTable.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>SteppedIterationTable</name>
      <anchorfile>classRVLAlg_1_1SteppedIterationTable.html</anchorfile>
      <anchor>a5121ce68ef17e136f03b29854da034de</anchor>
      <arglist>(FunctionalEvaluation&lt; Scalar &gt; &amp;_fx, Scalar &amp;_step, ostream &amp;_out=cout)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1SteppedIterationTable.html</anchorfile>
      <anchor>a24527fe2a805c58ba804870e10f86314</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>count</name>
      <anchorfile>classRVLAlg_1_1SteppedIterationTable.html</anchorfile>
      <anchor>a5fc4a09d7afa218e91cbc4bd82561764</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>FunctionalEvaluation&lt; Scalar &gt; &amp;</type>
      <name>fx</name>
      <anchorfile>classRVLAlg_1_1SteppedIterationTable.html</anchorfile>
      <anchor>a9e4c1d4cb0dc172cc71ae974b6c66de5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Scalar &amp;</type>
      <name>step</name>
      <anchorfile>classRVLAlg_1_1SteppedIterationTable.html</anchorfile>
      <anchor>a11ad2e04b978e3ab6143e9c5413e3af8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ostream &amp;</type>
      <name>out</name>
      <anchorfile>classRVLAlg_1_1SteppedIterationTable.html</anchorfile>
      <anchor>a5178652a2a0b0aba33c6644d980866bf</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::GradientThresholdIterationTable</name>
    <filename>classRVLAlg_1_1GradientThresholdIterationTable.html</filename>
    <templarg></templarg>
    <base>RVLAlg::CountTerminator</base>
    <member kind="function">
      <type></type>
      <name>GradientThresholdIterationTable</name>
      <anchorfile>classRVLAlg_1_1GradientThresholdIterationTable.html</anchorfile>
      <anchor>adb9b5c1982b751fd18c538efdffb5ca7</anchor>
      <arglist>(FunctionalEvaluation&lt; Scalar &gt; &amp;_fx, int maxcount=1, atype _atol=numeric_limits&lt; atype &gt;::max(), atype _rtol=ScalarFieldTraits&lt; atype &gt;::One(), ostream &amp;_out=cout, bool _doit=true)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1GradientThresholdIterationTable.html</anchorfile>
      <anchor>ab126627d175929d254445ec61d484def</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>atype</type>
      <name>atol</name>
      <anchorfile>classRVLAlg_1_1GradientThresholdIterationTable.html</anchorfile>
      <anchor>a2f995e39da3fea04a8f40340bb6f9859</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>atype</type>
      <name>rtol</name>
      <anchorfile>classRVLAlg_1_1GradientThresholdIterationTable.html</anchorfile>
      <anchor>ac1403becf4df23ccd013c0157a145a20</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>init</name>
      <anchorfile>classRVLAlg_1_1GradientThresholdIterationTable.html</anchorfile>
      <anchor>ac724bfcf45353529ecdbf32bd2ca7321</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>atype</type>
      <name>fxval0</name>
      <anchorfile>classRVLAlg_1_1GradientThresholdIterationTable.html</anchorfile>
      <anchor>ae8a87f9d09e6e4af0490bbede1eff55f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>atype</type>
      <name>ngfx0</name>
      <anchorfile>classRVLAlg_1_1GradientThresholdIterationTable.html</anchorfile>
      <anchor>a6c73472f69383080c398471e9dbf7680</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>atype</type>
      <name>ingfx</name>
      <anchorfile>classRVLAlg_1_1GradientThresholdIterationTable.html</anchorfile>
      <anchor>a5389e683f557c7765146d91f4a77516a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>doit</name>
      <anchorfile>classRVLAlg_1_1GradientThresholdIterationTable.html</anchorfile>
      <anchor>afeba5c394d742e7fecff6e085e10ea0a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>FunctionalEvaluation&lt; Scalar &gt; &amp;</type>
      <name>fx</name>
      <anchorfile>classRVLAlg_1_1GradientThresholdIterationTable.html</anchorfile>
      <anchor>af354ecd7823a2ce04dd74f9dbaa32350</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ostream &amp;</type>
      <name>out</name>
      <anchorfile>classRVLAlg_1_1GradientThresholdIterationTable.html</anchorfile>
      <anchor>a5017d52a8c4efd4b6997deee8cc2d106</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::CountingThresholdIterationTable</name>
    <filename>classRVLAlg_1_1CountingThresholdIterationTable.html</filename>
    <templarg></templarg>
    <base>RVLAlg::CountTerminator</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>atype</name>
      <anchorfile>classRVLAlg_1_1CountingThresholdIterationTable.html</anchorfile>
      <anchor>a65c7eea3b99aa4021dea1c49a488fd50</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CountingThresholdIterationTable</name>
      <anchorfile>classRVLAlg_1_1CountingThresholdIterationTable.html</anchorfile>
      <anchor>a5eb3be0df85d0a7d70f2d9c7dfd84aaa</anchor>
      <arglist>(int maxcount, const atype &amp;val, atype tol, string &amp;scalarname, ostream &amp;out=cout)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CountingThresholdIterationTable</name>
      <anchorfile>classRVLAlg_1_1CountingThresholdIterationTable.html</anchorfile>
      <anchor>ad62da47203fbecd721e70e50330ef4a9</anchor>
      <arglist>(int maxcount, const atype &amp;val, atype tol, ostream &amp;out=cout)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1CountingThresholdIterationTable.html</anchorfile>
      <anchor>a0fafead714261e4dc13ff1af799c7e14</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>atype</type>
      <name>tol_</name>
      <anchorfile>classRVLAlg_1_1CountingThresholdIterationTable.html</anchorfile>
      <anchor>a9162e217c8d62e0dab3ad6b4af58bd21</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const atype &amp;</type>
      <name>val_</name>
      <anchorfile>classRVLAlg_1_1CountingThresholdIterationTable.html</anchorfile>
      <anchor>a3dc8e9ffbb9f9ac7735f0fce3d7cabd0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ostream &amp;</type>
      <name>out_</name>
      <anchorfile>classRVLAlg_1_1CountingThresholdIterationTable.html</anchorfile>
      <anchor>a84bc9ce5d17c205aaeebee3d09707af1</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::VectorCountingThresholdIterationTable</name>
    <filename>classRVLAlg_1_1VectorCountingThresholdIterationTable.html</filename>
    <templarg></templarg>
    <base>RVLAlg::CountTerminator</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>atype</name>
      <anchorfile>classRVLAlg_1_1VectorCountingThresholdIterationTable.html</anchorfile>
      <anchor>a629e10aadcbd30754ebbf49af4a02f8a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>VectorCountingThresholdIterationTable</name>
      <anchorfile>classRVLAlg_1_1VectorCountingThresholdIterationTable.html</anchorfile>
      <anchor>a7e72b062e6d00387509ea7d0ba393016</anchor>
      <arglist>(int maxcount, vector&lt; string &gt; const &amp;_names, vector&lt; atype * &gt; const &amp;_nums, vector&lt; atype &gt; const &amp;_tols, ostream &amp;_out=cout)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init</name>
      <anchorfile>classRVLAlg_1_1VectorCountingThresholdIterationTable.html</anchorfile>
      <anchor>a3a267ecf1c914fafbc59c46ebd5f9ebe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1VectorCountingThresholdIterationTable.html</anchorfile>
      <anchor>a11c6244249c7180294d24cc3f89af974</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>vector&lt; atype &gt;</type>
      <name>tols</name>
      <anchorfile>classRVLAlg_1_1VectorCountingThresholdIterationTable.html</anchorfile>
      <anchor>ae22046477c06b836f475bd94fe3460a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>vector&lt; atype * &gt;</type>
      <name>nums</name>
      <anchorfile>classRVLAlg_1_1VectorCountingThresholdIterationTable.html</anchorfile>
      <anchor>a224e3369e024ae0dca37f81b14ab9ee1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>vector&lt; string &gt;</type>
      <name>names</name>
      <anchorfile>classRVLAlg_1_1VectorCountingThresholdIterationTable.html</anchorfile>
      <anchor>acf4dd9c1663911a7c52a51817c73afcd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ostream &amp;</type>
      <name>out</name>
      <anchorfile>classRVLAlg_1_1VectorCountingThresholdIterationTable.html</anchorfile>
      <anchor>aa18edc93d563e473775c86859bc5361e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::CountingNormTable</name>
    <filename>classRVLAlg_1_1CountingNormTable.html</filename>
    <templarg></templarg>
    <base>RVLAlg::CountTerminator</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>NormRetType</name>
      <anchorfile>classRVLAlg_1_1CountingNormTable.html</anchorfile>
      <anchor>a155512c5f7854d37893e98f12f98b737</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CountingNormTable</name>
      <anchorfile>classRVLAlg_1_1CountingNormTable.html</anchorfile>
      <anchor>ae071df766e456e05cb1c96f6e7e78460</anchor>
      <arglist>(int maxcount, const Vector&lt; Scalar &gt; &amp;_x, NormRetType _tol, bool _doit=true, ostream &amp;_out=cout)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1CountingNormTable.html</anchorfile>
      <anchor>aa86ba696df92d8b31fd355bfa81255c4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>doit</name>
      <anchorfile>classRVLAlg_1_1CountingNormTable.html</anchorfile>
      <anchor>a87711469d00c394a41cb2271f85625c7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ostream &amp;</type>
      <name>out</name>
      <anchorfile>classRVLAlg_1_1CountingNormTable.html</anchorfile>
      <anchor>a4b2be54a0048303f34103bc719f5f7cb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const Vector&lt; Scalar &gt; &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1CountingNormTable.html</anchorfile>
      <anchor>a7b4ed1e973798a18436c6a238387f408</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>NormRetType</type>
      <name>tol</name>
      <anchorfile>classRVLAlg_1_1CountingNormTable.html</anchorfile>
      <anchor>ae9a1f48ebc460862b0101c74c4e039ea</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::CountTerminator</name>
    <filename>classRVLAlg_1_1CountTerminator.html</filename>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>CountTerminator</name>
      <anchorfile>classRVLAlg_1_1CountTerminator.html</anchorfile>
      <anchor>a89229ea4ef2d0775dad8bf25906288b2</anchor>
      <arglist>(int maxcount)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CountTerminator</name>
      <anchorfile>classRVLAlg_1_1CountTerminator.html</anchorfile>
      <anchor>a23cf99db88b79d9620d6608a937d8f57</anchor>
      <arglist>(int init, int maxcount, int inc)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CountTerminator</name>
      <anchorfile>classRVLAlg_1_1CountTerminator.html</anchorfile>
      <anchor>a0debe8b80fd817586238db109193e09a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1CountTerminator.html</anchorfile>
      <anchor>ad64dce961cbfa14fe1cb97a813525b76</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getCount</name>
      <anchorfile>classRVLAlg_1_1CountTerminator.html</anchorfile>
      <anchor>ab8d8b7849d2991bd8170367a9d78159d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>mc</name>
      <anchorfile>classRVLAlg_1_1CountTerminator.html</anchorfile>
      <anchor>a9b682bf3fe6e6daf99bf7b96054a8d87</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>count</name>
      <anchorfile>classRVLAlg_1_1CountTerminator.html</anchorfile>
      <anchor>a4e66873bc2f5bfa0e42ab99e01e17f76</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>i</name>
      <anchorfile>classRVLAlg_1_1CountTerminator.html</anchorfile>
      <anchor>a597ad03c1ba8c5795cec101ffef4f26b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::MaxTerminator</name>
    <filename>classRVLAlg_1_1MaxTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>MaxTerminator</name>
      <anchorfile>classRVLAlg_1_1MaxTerminator.html</anchorfile>
      <anchor>a0d86eb00fc82cf14a4f317572f693b7e</anchor>
      <arglist>(const Scalar &amp;tx, Scalar maxval)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1MaxTerminator.html</anchorfile>
      <anchor>ac8de37c16b0fecb9b505fc364a405972</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const Scalar &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1MaxTerminator.html</anchorfile>
      <anchor>a15aaffea15ae5be3a4411ffec7e48d5b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Scalar</type>
      <name>mv</name>
      <anchorfile>classRVLAlg_1_1MaxTerminator.html</anchorfile>
      <anchor>a1ac70572ffad0f9274f7865afce42d99</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::MinTerminator</name>
    <filename>classRVLAlg_1_1MinTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>MinTerminator</name>
      <anchorfile>classRVLAlg_1_1MinTerminator.html</anchorfile>
      <anchor>a4e2d0cb5c9917dd40855f5392cc0aabe</anchor>
      <arglist>(const Scalar &amp;tx, Scalar minval)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1MinTerminator.html</anchorfile>
      <anchor>a31fa187f1250a5078b588ea9a0801010</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const Scalar &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1MinTerminator.html</anchorfile>
      <anchor>aab98eb01ce1434653d6e8816e4a0cf14</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Scalar</type>
      <name>mv</name>
      <anchorfile>classRVLAlg_1_1MinTerminator.html</anchorfile>
      <anchor>a16fe363a76aecff15630643b4674b7a2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::MinTerminatorFE</name>
    <filename>classRVLAlg_1_1MinTerminatorFE.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>MinTerminatorFE</name>
      <anchorfile>classRVLAlg_1_1MinTerminatorFE.html</anchorfile>
      <anchor>a72b852dcca4ed7f4f2d7a2edc7e0e4bd</anchor>
      <arglist>(RVL::FunctionalEvaluation&lt; Scalar &gt; &amp;_fx, Scalar minval)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1MinTerminatorFE.html</anchorfile>
      <anchor>af39ea820d6dd531bbbf41dbfe8e1b47b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Scalar</type>
      <name>mv</name>
      <anchorfile>classRVLAlg_1_1MinTerminatorFE.html</anchorfile>
      <anchor>a55033f13fb461b8623a351c89f970172</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::UnaryThresholdTerminator</name>
    <filename>classRVLAlg_1_1UnaryThresholdTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>UnaryThresholdTerminator</name>
      <anchorfile>classRVLAlg_1_1UnaryThresholdTerminator.html</anchorfile>
      <anchor>ac7a5a6e6d86b948ba809bb5dd3308f94</anchor>
      <arglist>(FunctionObjectScalarRedn&lt; Scalar &gt; &amp;tf, Vector&lt; Scalar &gt; &amp;tx, Scalar ttol)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1UnaryThresholdTerminator.html</anchorfile>
      <anchor>a6a189ca5084d0591602b58898f5918d2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>FunctionObjectScalarRedn&lt; Scalar &gt; &amp;</type>
      <name>f</name>
      <anchorfile>classRVLAlg_1_1UnaryThresholdTerminator.html</anchorfile>
      <anchor>abc48944eed788079735b91b551e7da07</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1UnaryThresholdTerminator.html</anchorfile>
      <anchor>a86781de06b6ff1a31f281b8e113e983b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Scalar</type>
      <name>tol</name>
      <anchorfile>classRVLAlg_1_1UnaryThresholdTerminator.html</anchorfile>
      <anchor>a38d30e4a89e9fd8de35fca31ec50f62b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::BinaryThresholdTerminator</name>
    <filename>classRVLAlg_1_1BinaryThresholdTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>BinaryThresholdTerminator</name>
      <anchorfile>classRVLAlg_1_1BinaryThresholdTerminator.html</anchorfile>
      <anchor>aa6db560998626485da9fd177671909ca</anchor>
      <arglist>(FunctionObjectScalarRedn&lt; Scalar &gt; &amp;tf, Vector&lt; Scalar &gt; &amp;tx, Vector&lt; Scalar &gt; &amp;ty, Scalar ttol)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1BinaryThresholdTerminator.html</anchorfile>
      <anchor>a3e3f66e22bb551e650122b819edc5301</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>FunctionObjectScalarRedn&lt; Scalar &gt; &amp;</type>
      <name>f</name>
      <anchorfile>classRVLAlg_1_1BinaryThresholdTerminator.html</anchorfile>
      <anchor>a26b4f459036f9ad7931666f616cf331d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1BinaryThresholdTerminator.html</anchorfile>
      <anchor>a893c8982949f894d9d4e2fcaaa8bd5f6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>y</name>
      <anchorfile>classRVLAlg_1_1BinaryThresholdTerminator.html</anchorfile>
      <anchor>a88b832555e2e953305dd2d2fc5324860</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Scalar</type>
      <name>tol</name>
      <anchorfile>classRVLAlg_1_1BinaryThresholdTerminator.html</anchorfile>
      <anchor>a11d3f7a28ee8792adef7dc120ec1605b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::TernaryThresholdTerminator</name>
    <filename>classRVLAlg_1_1TernaryThresholdTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>TernaryThresholdTerminator</name>
      <anchorfile>classRVLAlg_1_1TernaryThresholdTerminator.html</anchorfile>
      <anchor>a9552cdd40b6b549ab7085ee3000ea669</anchor>
      <arglist>(FunctionObjectScalarRedn&lt; Scalar &gt; &amp;tf, Vector&lt; Scalar &gt; &amp;tx, Vector&lt; Scalar &gt; &amp;ty, Vector&lt; Scalar &gt; &amp;tz, Scalar ttol)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1TernaryThresholdTerminator.html</anchorfile>
      <anchor>a102a7e36969fb7812a7e784f27a37569</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>FunctionObjectScalarRedn&lt; Scalar &gt; &amp;</type>
      <name>f</name>
      <anchorfile>classRVLAlg_1_1TernaryThresholdTerminator.html</anchorfile>
      <anchor>af4f83e3aeaa4b2b58163790bad1c199f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1TernaryThresholdTerminator.html</anchorfile>
      <anchor>a773ea732ee1d7a0f3fa5c2ea4d88499a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>y</name>
      <anchorfile>classRVLAlg_1_1TernaryThresholdTerminator.html</anchorfile>
      <anchor>a2a3d8b2d556ba5d3582ae668d7eaac61</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>z</name>
      <anchorfile>classRVLAlg_1_1TernaryThresholdTerminator.html</anchorfile>
      <anchor>a6cc9a2f0e8807193695042ee38f31524</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Scalar</type>
      <name>tol</name>
      <anchorfile>classRVLAlg_1_1TernaryThresholdTerminator.html</anchorfile>
      <anchor>a2995418aed5a11decec1c5ceb5c76294</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::NormThresholdTerminator</name>
    <filename>classRVLAlg_1_1NormThresholdTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>NormRetType</name>
      <anchorfile>classRVLAlg_1_1NormThresholdTerminator.html</anchorfile>
      <anchor>a8462299180830dd960ffde52147f945d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>NormThresholdTerminator</name>
      <anchorfile>classRVLAlg_1_1NormThresholdTerminator.html</anchorfile>
      <anchor>ae6a088a881ebbbae1647eaab2ca0895c</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;tx, NormRetType ttol)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1NormThresholdTerminator.html</anchorfile>
      <anchor>a8be061819c4c9cd7298bbd7f0bf72b8e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const Vector&lt; Scalar &gt; &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1NormThresholdTerminator.html</anchorfile>
      <anchor>aa380ec01f00dd5560e1c3d0977930869</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>NormRetType</type>
      <name>tol</name>
      <anchorfile>classRVLAlg_1_1NormThresholdTerminator.html</anchorfile>
      <anchor>a4457a4e41539ed5f3a102a25628dde08</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::Norm2ThresholdTerminator</name>
    <filename>classRVLAlg_1_1Norm2ThresholdTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>NormRetType</name>
      <anchorfile>classRVLAlg_1_1Norm2ThresholdTerminator.html</anchorfile>
      <anchor>a650c49af70ce175023f7a18ee34ebd85</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Norm2ThresholdTerminator</name>
      <anchorfile>classRVLAlg_1_1Norm2ThresholdTerminator.html</anchorfile>
      <anchor>aef54f7d45c29f6a2cba85f1fe437c458</anchor>
      <arglist>(Vector&lt; Scalar &gt; &amp;tx, NormRetType ttol)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1Norm2ThresholdTerminator.html</anchorfile>
      <anchor>a65dabddeda0bc4cc7eeaf9f79d637054</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1Norm2ThresholdTerminator.html</anchorfile>
      <anchor>a4858b6b5a9363e5dc05a853fe7490d40</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>NormRetType</type>
      <name>tol</name>
      <anchorfile>classRVLAlg_1_1Norm2ThresholdTerminator.html</anchorfile>
      <anchor>a8b23402e8ddd692cd033424ce1efa132</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::DiffThresholdTerminator</name>
    <filename>classRVLAlg_1_1DiffThresholdTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>NormRetType</name>
      <anchorfile>classRVLAlg_1_1DiffThresholdTerminator.html</anchorfile>
      <anchor>a7ba0e6b183c482593ca1768b28e70f88</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DiffThresholdTerminator</name>
      <anchorfile>classRVLAlg_1_1DiffThresholdTerminator.html</anchorfile>
      <anchor>aea7c28bd5d8cef067ae58014e4b7d1bf</anchor>
      <arglist>(Vector&lt; Scalar &gt; &amp;tx, Vector&lt; Scalar &gt; &amp;ty, NormRetType ttol)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1DiffThresholdTerminator.html</anchorfile>
      <anchor>a1f1045f3e8cf8a9034c7cb702b3a8767</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1DiffThresholdTerminator.html</anchorfile>
      <anchor>a32693b2e180bbdaaa2e10508cccb263d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>y</name>
      <anchorfile>classRVLAlg_1_1DiffThresholdTerminator.html</anchorfile>
      <anchor>a6bf24a98d1a6d4a7d310d346f9680005</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>NormRetType</type>
      <name>tol</name>
      <anchorfile>classRVLAlg_1_1DiffThresholdTerminator.html</anchorfile>
      <anchor>a860f8603356dc3dcea0ad74954c84185</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::Diff2ThresholdTerminator</name>
    <filename>classRVLAlg_1_1Diff2ThresholdTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>NormRetType</name>
      <anchorfile>classRVLAlg_1_1Diff2ThresholdTerminator.html</anchorfile>
      <anchor>ac749fd6c69e4e52cef8d98657cc5b15a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Diff2ThresholdTerminator</name>
      <anchorfile>classRVLAlg_1_1Diff2ThresholdTerminator.html</anchorfile>
      <anchor>a2b17bf5d7420eaa90403ecc060fcedbb</anchor>
      <arglist>(Vector&lt; Scalar &gt; &amp;tx, Vector&lt; Scalar &gt; &amp;ty, NormRetType ttol)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1Diff2ThresholdTerminator.html</anchorfile>
      <anchor>a6114230f1f1f4817e209363d8c2a20f7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1Diff2ThresholdTerminator.html</anchorfile>
      <anchor>a857413532b9a9292bd737609926ab312</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>y</name>
      <anchorfile>classRVLAlg_1_1Diff2ThresholdTerminator.html</anchorfile>
      <anchor>a4187b3fe4e8fe10a3f9dbde9a0c20e15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>NormRetType</type>
      <name>tol</name>
      <anchorfile>classRVLAlg_1_1Diff2ThresholdTerminator.html</anchorfile>
      <anchor>aa7568c0604ecffe84fafdbf94909264e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::IPThresholdTerminator</name>
    <filename>classRVLAlg_1_1IPThresholdTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="function">
      <type></type>
      <name>IPThresholdTerminator</name>
      <anchorfile>classRVLAlg_1_1IPThresholdTerminator.html</anchorfile>
      <anchor>a9ef0746684aec4ba0d80b336ce52d986</anchor>
      <arglist>(Vector&lt; Scalar &gt; &amp;tx, Vector&lt; Scalar &gt; &amp;ty, Scalar ttol)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1IPThresholdTerminator.html</anchorfile>
      <anchor>a140a61ef7a4ae25b91762e0ab49fa065</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1IPThresholdTerminator.html</anchorfile>
      <anchor>ab30865c38630387d5d42133c4638310e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>y</name>
      <anchorfile>classRVLAlg_1_1IPThresholdTerminator.html</anchorfile>
      <anchor>ae3741f46e4ba0b6ee4261accac15e369</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Scalar</type>
      <name>tol</name>
      <anchorfile>classRVLAlg_1_1IPThresholdTerminator.html</anchorfile>
      <anchor>a74ed9b19071b79fb9508dbccb68b98b8</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::AbsIPThresholdTerminator</name>
    <filename>classRVLAlg_1_1AbsIPThresholdTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>NormRetType</name>
      <anchorfile>classRVLAlg_1_1AbsIPThresholdTerminator.html</anchorfile>
      <anchor>a78d16262cf6184d979f6c93f6fad4dc5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>AbsIPThresholdTerminator</name>
      <anchorfile>classRVLAlg_1_1AbsIPThresholdTerminator.html</anchorfile>
      <anchor>af3f5d7de27f80cedf0e61f6073ea73a0</anchor>
      <arglist>(Vector&lt; Scalar &gt; &amp;tx, Vector&lt; Scalar &gt; &amp;ty, NormRetType ttol)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AbsIPThresholdTerminator</name>
      <anchorfile>classRVLAlg_1_1AbsIPThresholdTerminator.html</anchorfile>
      <anchor>af1fa8183b4e4ca9c577e49cf5353ea4c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1AbsIPThresholdTerminator.html</anchorfile>
      <anchor>a2503d241882b3cc8734e3a3b0fcc831a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1AbsIPThresholdTerminator.html</anchorfile>
      <anchor>a716e88acc893f6c3ea1f5b435b68c872</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>y</name>
      <anchorfile>classRVLAlg_1_1AbsIPThresholdTerminator.html</anchorfile>
      <anchor>a73960a8ec7a3f98b5be344b6a9fb01bf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>NormRetType</type>
      <name>tol</name>
      <anchorfile>classRVLAlg_1_1AbsIPThresholdTerminator.html</anchorfile>
      <anchor>ab96f7c7e6a50244c9033b7ab53dbdb69</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::NormGradientTerminator</name>
    <filename>classRVLAlg_1_1NormGradientTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>NormRetType</name>
      <anchorfile>classRVLAlg_1_1NormGradientTerminator.html</anchorfile>
      <anchor>afd488ff65db4668e6fd1388efd439804</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>NormGradientTerminator</name>
      <anchorfile>classRVLAlg_1_1NormGradientTerminator.html</anchorfile>
      <anchor>a9742084bdc6e50bc342a3b63a0f6ff46</anchor>
      <arglist>(Vector&lt; Scalar &gt; &amp;x, Functional&lt; Scalar &gt; &amp;f, NormRetType tol)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1NormGradientTerminator.html</anchorfile>
      <anchor>ab0cbd46b6a2c7ab8a3b89d984078c015</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>NormRetType</type>
      <name>tol_</name>
      <anchorfile>classRVLAlg_1_1NormGradientTerminator.html</anchorfile>
      <anchor>ae4678fe256f2e7e9107e2abb7eec5ccd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>FunctionalEvaluation&lt; Scalar &gt;</type>
      <name>fx_</name>
      <anchorfile>classRVLAlg_1_1NormGradientTerminator.html</anchorfile>
      <anchor>a54871c84c876c06ce67b93a141c0e234</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::DiffBallProjTerminator</name>
    <filename>classRVLAlg_1_1DiffBallProjTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>NormRetType</name>
      <anchorfile>classRVLAlg_1_1DiffBallProjTerminator.html</anchorfile>
      <anchor>a81c6be7c5e35da6384cd3fc24f96c2ec</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DiffBallProjTerminator</name>
      <anchorfile>classRVLAlg_1_1DiffBallProjTerminator.html</anchorfile>
      <anchor>a9fb5199ce72e0418d099c2821d51d2c3</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;tx, Vector&lt; Scalar &gt; &amp;ty, NormRetType _maxstep, ostream &amp;_str=cout)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1DiffBallProjTerminator.html</anchorfile>
      <anchor>ab9941f5216c2cfab585b1c1e6e8ec539</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>static_query</name>
      <anchorfile>classRVLAlg_1_1DiffBallProjTerminator.html</anchorfile>
      <anchor>a1326357d9bb2aa409b6c98fec2bb82bb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; const &amp;</type>
      <name>x</name>
      <anchorfile>classRVLAlg_1_1DiffBallProjTerminator.html</anchorfile>
      <anchor>af3296a9330cee5a50e20d1b6781914a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>y</name>
      <anchorfile>classRVLAlg_1_1DiffBallProjTerminator.html</anchorfile>
      <anchor>aaa4ddf3b3c122b3d3806869662838540</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>NormRetType</type>
      <name>maxstep</name>
      <anchorfile>classRVLAlg_1_1DiffBallProjTerminator.html</anchorfile>
      <anchor>a556cf66dd28ae2d65e0ff527eccb5878</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>res</name>
      <anchorfile>classRVLAlg_1_1DiffBallProjTerminator.html</anchorfile>
      <anchor>a1eced85ce28ee41ad5d2717460602cca</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ostream &amp;</type>
      <name>str</name>
      <anchorfile>classRVLAlg_1_1DiffBallProjTerminator.html</anchorfile>
      <anchor>a377e8e6fc34c40b42b8644726c34189e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVLAlg::BallProjTerminator</name>
    <filename>classRVLAlg_1_1BallProjTerminator.html</filename>
    <templarg></templarg>
    <base>RVLAlg::Terminator</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>NormRetType</name>
      <anchorfile>classRVLAlg_1_1BallProjTerminator.html</anchorfile>
      <anchor>a29d11610077865549d1c5b8cc789ae3e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BallProjTerminator</name>
      <anchorfile>classRVLAlg_1_1BallProjTerminator.html</anchorfile>
      <anchor>a0988095b9d07758502f2099174a0ea4f</anchor>
      <arglist>(Vector&lt; Scalar &gt; &amp;ty, NormRetType _maxstep, ostream &amp;_str=cout)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>query</name>
      <anchorfile>classRVLAlg_1_1BallProjTerminator.html</anchorfile>
      <anchor>a716d10eb7ea90a5cf624f60e1d61d349</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>y</name>
      <anchorfile>classRVLAlg_1_1BallProjTerminator.html</anchorfile>
      <anchor>a92f65ceb59ce1422f36ff259eeb2ebf7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>NormRetType</type>
      <name>maxstep</name>
      <anchorfile>classRVLAlg_1_1BallProjTerminator.html</anchorfile>
      <anchor>a8cabe4d0d8857c8428af4a3951956908</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ostream &amp;</type>
      <name>str</name>
      <anchorfile>classRVLAlg_1_1BallProjTerminator.html</anchorfile>
      <anchor>af929a1773cda0d18a8f08d07d413c0e8</anchor>
      <arglist></arglist>
    </member>
  </compound>
</tagfile>
