<p align="center">
  <img src=http://www.ahay.org//wikilocal/style/Madagascar2.png>
</p>

Madagascar
==========
[![CircleCI](https://img.shields.io/circleci/project/github/ahay/src/master.svg?label=Circle%20CI)](https://circleci.com/gh/ahay/src)
[![TravisCI](https://img.shields.io/travis/com/ahay/src.svg?label=Travis%20CI)](https://travis-ci.com/ahay/src)
###  A package for reproducible geophysical data processing

http://www.ahay.org

## What is Madagascar?

Madagascar is an open-source software package for multidimensional data
analysis and reproducible computational experiments. Its mission is to
provide

* a convenient and powerful environment
* a convenient technology transfer tool 

for researchers working with digital image and data processing in
geophysics and related fields. Technology developed using the
Madagascar project management system is transferred in the form of
recorded processing histories, which become "computational recipes" to
be verified, exchanged, and modified by users of the system.


## Design Principles

* Madagascar is a modern package. Started in 2003 and publicly released in
2006, it was developed entirely from scratch. Being a relatively new package,
it follows modern software engineering practices such as module
encapsulation and test-driven development. A rapid development of a project
of this scope (more than 300 main programs and more than 3,000 tests) would
not be possible without standing on the shoulders of giants and learning
from the 30 years of previous experience in open packages such as SEPlib and
Seismic Unix. We have borrowed and reimplemented functionality and ideas
from these other packages.

* Madagascar is a test-driven package. Test-driven development is not only an
agile software programming practice but also a way of bringing scientific
foundation to geophysical research that involves numerical experiments.
Bringing reproducibility and peer review, the backbone of any real science,
to the field of computational geophysics is the main motivation for
Madagascar development. The package consists of two levels: low-level main
programs (typically developed in the C programming language and working as
data filters) and high-level processing flows (described with the help of
the Python programming language) that combine main programs and completely
document data processing histories for testing and reproducibility.
Experience shows that high-level programming is easily mastered even by
beginning students that have no previous programming experience.

* Madagascar is an open-source package. It is distributed under the standard
GPL open-source license, which places no restriction on the usage and
modification of the code. Moreover, access to modifying the source
repository is not controlled by one organization but shared equally among
different developers. This enables an open collaboration among different
groups spread all over the world, in the true spirit of the open source
movement.

* Madagascar uses a simple, flexible, and universal data format that can
handle very large datasets but is not tied specifically to seismic data or
data of any other particular kind. This "regularly sampled" format is
borrowed from the traditional SEPlib. A universal data format allows us to
share general-purpose data processing tools with scientists from other
disciplines such as petroleum engineers working on large-scale reservoir
simulations.

## Where to get more information about Madagascar

The primary source of information is the web site:

http://www.ahay.org/

Additional information:

users' mailing list ("RSF-user"):

https://lists.sourceforge.net/lists/listinfo/rsf-user

developers' mailing list ("RSF-devel"):

https://lists.sourceforge.net/lists/listinfo/rsf-devel

development blog:

http://ahay.org/blog/

## Compiling, Building, Installing and Testing

See the INSTALL.txt document for build instructions.

## History

While being written from scratch, Madagascar borrows ideas from the
design of SEPlib, a publicly available software package, maintained by
Bob Clapp at the Stanford Exploration Project (SEP). Generations of
SEP students and researchers contributed to SEPlib. Most important
contributions came from Rob Clayton, Jon Claerbout, Dave Hale, Stew
Levin, Rick Ottolini, Joe Dellinger, Steve Cole, Dave Nichols, Martin
Karrenbach, Biondo Biondi, and Bob Clapp.

Madagascar was started, under the name RSF (Regularly Sampled Format)
by Sergey Fomel in 2003. Since then, many other people have
contributed to it. See the AUTHORS.txt file for an incomplete list.
