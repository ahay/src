# Spec file for making madagascar 1.2 RPMs for Fedora 15
# Usage details at http://ahay.org/wiki/Packaging_madagascar#RPM

%define version 1.2
# If you change version number, remember to also change it in the sed hack in the install section
%define m8rv madagascar-%{version}

Name:      madagascar
Version:   %{version}
Release:   1%{?dist}
License:   GPLv2+
Summary:   Utilities for geophysical data processing and numerical experiments
Group:     Applications/Engineering
URL:       http://ahay.org
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
BuildArch: %{_arch}
Requires:  binutils, gcc, glibc-headers, scons, gcc-c++, gcc-gfortran, numpy, python, swig, libgomp, blas, blas-devel, atlas, atlas-devel, libXaw-devel 

%description
Madagascar is an open-source software package for multidimensional data analysis and reproducible computational experiments. Its mission is to provide a convenient and powerful work environment and technology transfer tool for researchers working with digital image and data processing in geophysics and related fields.

%prep
rm -rf src
svn export https://rsf.svn.sourceforge.net/svnroot/rsf/branches/madagascar-1.2 src
# Patches, only for 1.2
svn export https://rsf.svn.sourceforge.net/svnroot/rsf/trunk/pens/genlib/ -r 7483 genlib_patch
mv -f genlib_patch/{SConstruct,gentext.c} src/pens/genlib/
rm -rf genlib_patch
svn export https://rsf.svn.sourceforge.net/svnroot/rsf/trunk/plot/plplot -r 7484 plsurf_patch
mv -f plsurf_patch/{SConstruct,plsurf.c} src/plot/plplot/
rm -rf plsurf_patch
svn export https://rsf.svn.sourceforge.net/svnroot/rsf/trunk/framework -r 7480 blas_patch
mv -f blas_patch/configure.py src/framework/
rm -rf blas_patch

%build
cd src
./configure --prefix=%{buildroot}/usr DYNLIB=y XLIBS=Xaw,Xt,X11 API=f77,f90,c++,python
make

%install
rm -rf %{buildroot}
cd src
make install

# Hack -- temporary fix. Should parse buildroot to add backslashes instead, or
# add a "rpm" option to configure
sed -i "s/\/home\/makerpm\/rpmbuild\/BUILDROOT\/madagascar-1\.2-1\.fc15\.x86_64//g" %{buildroot}/usr/share/madagascar/etc/env.sh
sed -i "s/\/home\/makerpm\/rpmbuild\/BUILDROOT\/madagascar-1\.2-1\.fc15\.x86_64//g" %{buildroot}/usr/share/madagascar/etc/env.csh
sed -i "s/\/home\/makerpm\/rpmbuild\/BUILDROOT\/madagascar-1\.2-1\.fc15\.x86_64//g" %{buildroot}/usr/share/madagascar/etc/config.py
sed -i "s/\/home\/makerpm\/rpmbuild\/BUILDROOT\/madagascar-1\.2-1\.fc15\.x86_64//g" %{buildroot}/usr/lib/python2.7/site-packages/rsf/prog.py

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root,-)
%{_bindir}/*
# Using _exec_prefix/lib instead of _libdir because on x86_64 _libdir defaults 
# to lib64, and it would be difficult to change hard-coded RSFROOT/lib 
# references throughout the codebase right now
%{_exec_prefix}/lib/*
%{_includedir}/*
%{_mandir}/man1/*
%{_mandir}/man9/*
%dir %{_datadir}/doc/madagascar
%dir %{_datadir}/doc/madagascar/spec
%dir %{_datadir}/doc/madagascar/txt
%dir %{_datadir}/madagascar
%dir %{_datadir}/madagascar/etc
%{_datadir}/doc/madagascar/spec/*
%{_datadir}/doc/madagascar/txt/*
%{_datadir}/madagascar/etc/*
%doc %{_datadir}/doc/madagascar/*.txt
%docdir %{_datadir}/doc/madagascar/html
%{_datadir}/doc/madagascar/html
