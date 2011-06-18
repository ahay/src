# Spec file for making madagascar 1.0 RPMs for Fedora 12
# Usage details at http://ahay.org/wiki/Packaging_madagascar#RPM

%define version 1.0
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
svn export https://rsf.svn.sourceforge.net/svnroot/rsf/branches/madagascar-1.0 src

%build
cd src
./configure --prefix=%{buildroot}/usr DYNLIB=y XLIBS=Xaw,Xt,X11 API=f77,f90,c++,python
make

%install
rm -rf %{buildroot}
cd src
make install
# Fixing m8r differences from the Filesystem Hierarchy Standard:
mkdir -p %{buildroot}/usr/share/%{m8rv}/html
mv %{buildroot}/usr/share/doc/*.html %{buildroot}/usr/share/%{m8rv}/html
mv %{buildroot}/usr/share/{spec,txt} %{buildroot}/usr/share/%{m8rv}
mv %{buildroot}/usr/etc %{buildroot}/etc
# Bad hack -- temporary fix. Should parse buildroot to add backslashes instead, or add a "rpm" option to configure
sed -i "s/\/home\/makerpm\/rpmbuild\/BUILDROOT\/madagascar-1\.0-1\.fc12\.x86_64//g" %{buildroot}/etc/madagascar/env.sh
sed -i "s/\/home\/makerpm\/rpmbuild\/BUILDROOT\/madagascar-1\.0-1\.fc12\.x86_64//g" %{buildroot}/etc/madagascar/env.csh
sed -i "s/\/home\/makerpm\/rpmbuild\/BUILDROOT\/madagascar-1\.0-1\.fc12\.x86_64//g" %{buildroot}/usr/lib/python2.6/site-packages/rsf/prog.py

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
%{_datadir}/%{m8rv}/spec/*
%{_datadir}/%{m8rv}/txt/*
%dir %{_sysconfdir}/madagascar
%{_sysconfdir}/madagascar/*
%doc src/{AUTHORS,COPYING,README}.txt
%docdir %{_datadir}/%{m8rv}/html
%{_datadir}/%{m8rv}/html
