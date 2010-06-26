# For making m8r RPMs for Fedora 12

Name:       madagascar
Version:    1.0.0
Release:    1%{?dist}
License:    GPLv2+
Summary:    Utilities for geophysical data processing and numerical experiments
Group:      Applications/Engineering
URL:        http://reproducibility.org
Source0:    http://sourceforge.net/projects/rsf/files/%{name}/%{name}-%{version}/%{name}-%{version}.tar.gz
BuildRoot:  %{_tmppath}/%{name}-%{version}-%{release}-root
BuildArch:  x86_64
Requires:   binutils, gcc, glibc-headers, scons, texlive-latex, subversion gcc-c++, gcc-gfortran, numpy, python, swig, octave, octave-devel, libgomp, openmpi, openmpi-devel, blas, blas-devel, atlas, atlas-devel, units, gifsicle, ffmpeg, ffmpeg-devel, libtiff-devel, libjpeg-devel, plplot-devel, mesa-libGL-devel, freeglut, freeglut-devel, libXaw-devel, netpbm-devel 

%description
Utilities for geophysical data processing and numerical experiments

%prep
%setup -q

%build
%configure --prefix=$RPM_BUILD_ROOT
make

%install
rm -rf $RPM_BUILD_ROOT
make install

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root,-)
%doc COPYING README AUTHORS NEWS
%dir %{_datadir}/%{name}/
