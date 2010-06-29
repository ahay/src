# Spec file for making madagascar RPMs for Fedora 12
#
# To make the RPMs, it is a good idea to first create a user just for this,
# i.e. makerpm. Then, install dependencies and set up. As root:
# yum -y groupinstall "Development Tools"
# yum -y install rpmdevtools mock rpmlint 
# useradd -G mock makerpm
#
# Log in as makerpm. Type: rpmdev-setuptree to create the ~/rpmbuild directory.
# Place this file (m8r-fedora12-x86_64.spec) in ~/rpmbuild/SPECS , then cd to 
# that directory and run:
# rpmbuild -ba madagascar.spec
#
# If you want a log of the entire process, use something like
# rpmbuild -ba madagascar.spec |& tee ~/log_rpm.asc
#
# You may be able to QC the build with
# mock -r fedora-12-x86_64 rebuild path_to_source_RPM
#
# If you want to debug your experiments, you can comment the contents of the
# prep, build or clean sections (for prep and build, you must have executed it
# once uncommented)
#
# A very good reference is
# https://fedoraproject.org/wiki/How_to_create_an_RPM_package

%define version 1.0.0alpha
%define m8rv madagascar-%{version}

Name:      madagascar
Version:   %{version}
Release:   1%{?dist}
License:   GPLv2+
Summary:   Utilities for geophysical data processing and numerical experiments
Group:     Applications/Engineering
URL:       http://reproducibility.org
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
BuildArch: x86_64
Requires:  binutils, gcc, glibc-headers, scons, texlive-latex, gcc-c++, gcc-gfortran, numpy, python, swig, libgomp, openmpi, openmpi-devel, blas, blas-devel, atlas, atlas-devel, units, gifsicle, libtiff-devel, libjpeg-devel, plplot-devel, mesa-libGL-devel, freeglut, freeglut-devel, libXaw-devel, netpbm-devel 

%description
Madagascar is an open-source software package for multidimensional data analysis and reproducible computational experiments. Its mission is to provide a convenient and powerful work environment and technology transfer tool for researchers working with digital image and data processing in geophysics and related fields.

%prep
rm -rf src
svn export -r 6304 https://rsf.svn.sourceforge.net/svnroot/rsf/trunk src

%build
cd src
./configure --prefix=%{buildroot}/usr API=f77,f90,c++,python
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

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root,-)
%{_bindir}/*
# Using this instead of _libdir because that defaults to lib64, and it would
# be difficult to change hard-coded RSFROOT/lib references just for this: 
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
